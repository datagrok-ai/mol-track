from typing import Any, Dict, List
from app.services.search.field_resolver import FieldResolutionError, FieldResolver
from app.services.search.operators import SearchOperators
import app.models as models
from app.services.search.utils import JoinOrderingTool, sanitize_field_name


class QueryBuildError(Exception):
    """Custom exception for query building errors"""

    pass


class QueryBuilder:
    """Builds dynamic SQL queries from search requests"""

    def __init__(self, field_resolver: FieldResolver):
        self.field_resolver = field_resolver
        self.operators = SearchOperators()

        self.parameter_counter = 0

    def build_query(self, request: models.SearchRequest) -> Dict[str, Any]:
        """
        Builds complete SQL query from search request

        Query consists of it's base part and the filter subquery

        SELECT ...
        FROM ...
        WHERE EXISTS (filter_subquery)

        Returns:
            Dict with query components: {
                'sql': str,              # Main query SQL
                'params': Dict[str, Any] # Query parameters
            }
        """
        level = request.level
        schema = self.field_resolver.db_schema
        table_config = self.field_resolver.table_configs[level]

        # Build SQL parts for base query
        base_query_joins = JoinOrderingTool()
        output_info = self.build_base_sql_parts(request.output, level, base_query_joins)

        base_select_clause = output_info["select_clause"]
        group_by = output_info["group_by"]
        table = table_config["table"]
        alias = table_config["alias"]
        has_dynamic = output_info["has_dynamic"]

        if (group_by and f"{table}.id" not in request.output) or (not group_by and has_dynamic):
            group_by = [f"{alias}{alias}.id"] + group_by

        group_by_sql = f"GROUP BY {' ,'.join(group_by)} " if group_by else ""

        # Build FROM clause with primary table
        base_from_clause = f"{schema}.{table} {alias}{alias}"

        query_params = {}
        filter_sql = ""

        # Build sub query from filters
        if request.filter:
            sql_components = self.build_filter_sql_parts(request.filter, level)

            # Create query from the returned sql parts
            if sql_components["sql"]:
                filter_sql = f"WHERE {sql_components['sql']}"
                query_params.update(sql_components["params"])

        # Convert joins to list and remove duplicates
        base_joins = base_query_joins.getJoinSQL()
        # The main query needs to have a different alias compared to the subqueries
        base_joins = base_joins.replace(
            f" {table_config['alias']}.", f" {table_config['alias']}{table_config['alias']}."
        )

        # Build main query
        complete_sql = (
            f"SELECT {base_select_clause} "
            f"FROM {base_from_clause} "
            f"{base_joins} "
            f"{filter_sql} "
            f"{group_by_sql} "
            f"ORDER BY {table_config['alias']}{table_config['alias']}.id "
        )

        return {"sql": complete_sql.strip(), "params": query_params}

    def build_base_sql_parts(
        self, output_list: List[str], search_level: str, all_joins: JoinOrderingTool
    ) -> Dict[str, Any]:
        """
        Generate SQL query parts for base query.

        Returns:
            Dict with SQL components: select clause, where clause, group by
        """
        select_fields, list_of_aliases, conditions, group_by = [], [], [], []
        # has_dynamic is True if any of the output fields are dynamic
        has_dynamic = False

        for field_path in output_list:
            resolved = self.field_resolver.resolve_field(field_path, search_level, all_joins)

            # Create alias for the field
            field_alias = sanitize_field_name(field_path)

            select_fields.append(f"{resolved['sql_expression']} AS {field_alias}")

            if resolved["is_dynamic"]:
                has_dynamic = True
            else:
                list_of_aliases.append(field_alias)

            property_filters = resolved.get("property_filter", None)
            if not resolved["is_dynamic"] and property_filters:
                conditions.append(property_filters)

        # Combine conditions with operator
        operator = " AND "
        combined_sql = operator.join(conditions)

        if has_dynamic:
            group_by = list_of_aliases
        return {
            "select_clause": ", ".join(select_fields),
            "conditions": combined_sql,
            "group_by": group_by,
            "has_dynamic": has_dynamic,
        }

    def build_filter_sql_parts(self, filter_obj: models.Filter, level: str) -> Dict[str, Any]:
        """
        Builds SQL WHERE clause from filter

        Returns:
            Dict with: {'sql': str, 'params': Dict[str, Any]}
        """
        if isinstance(filter_obj, models.AtomicCondition):
            # Handle single atomic condition
            cond_info = self.build_condition(filter_obj, level)
            return {
                "sql": cond_info["sql"],
                "params": cond_info["params"],
            }

        elif isinstance(filter_obj, models.LogicalNode):
            # Handle logical node with multiple conditions
            conditions = []
            all_params = {}

            for condition in filter_obj.conditions:
                if isinstance(condition, models.AtomicCondition):
                    cond_info = self.build_condition(condition, level)
                    # Rename parameters to avoid conflicts
                    renamed_sql = cond_info["sql"]
                    renamed_params = {}
                    for old_name, value in cond_info["params"].items():
                        new_name = f"{old_name}_{self.parameter_counter}"
                        renamed_sql = renamed_sql.replace(f":{old_name}", f":{new_name}")
                        renamed_params[new_name] = value

                    conditions.append(renamed_sql)
                    all_params.update(renamed_params)
                    self.parameter_counter += 1

                elif isinstance(condition, models.LogicalNode):
                    # Recursive filter handling
                    nested_info = self.build_filter_sql_parts(condition, level)

                    conditions.append(f"({nested_info['sql']})")
                    all_params.update(nested_info["params"])

            # Combine conditions with operator
            operator = f" {filter_obj.operator.value} "
            combined_conditions = operator.join(conditions)

            return {"sql": combined_conditions, "params": all_params}

    def build_condition(self, condition: models.AtomicCondition, level: str) -> Dict[str, Any]:
        """
        Builds SQL parts for a single condition

        Returns:
            Dict with: {'sql': str, 'params': Dict[str, Any]}
        """
        try:
            joins = JoinOrderingTool()
            # Resolve field to SQL components
            field_info = self.field_resolver.resolve_field(condition.field, level, joins, True)

            # Validate operator and value
            self.operators.validate_operator_value(condition.operator, condition.value, condition.threshold)
            self.operators.validate_operands(condition.operator, condition.field)

            # Handle dynamic properties with property name filtering
            if field_info["is_dynamic"]:
                return self._build_dynamic_condition(field_info, condition)
            else:
                return self._build_direct_condition(field_info, condition)

        except FieldResolutionError as e:
            raise QueryBuildError(f"Field resolution error: {str(e)}")
        except ValueError as e:
            raise QueryBuildError(f"Condition validation error: {str(e)}")

    def _build_direct_condition(self, field_info: Dict[str, Any], condition: models.AtomicCondition) -> Dict[str, Any]:
        """Build condition for direct field access"""
        field_sql = field_info["sql_expression"]

        sql_expr, params = self.operators.get_sql_expression(
            condition.operator, field_sql, condition.value, condition.threshold
        )

        sql = sql_expr
        if field_info["subquery"]["sql"] != "":
            if field_info["search_level"]["alias"] != field_info["subquery"]["alias"]:
                key = field_info["search_level"]["foreign_key"]
            else:
                key = "id"
            sql = (
                f"EXISTS ( "
                f"{field_info['subquery']['sql']} "
                f"WHERE {field_info['subquery']['alias']}.{key}="
                f"{field_info['search_level']['alias']}{field_info['search_level']['alias']}.id "
                f"AND {sql_expr})  "
            )

        return {"sql": sql, "params": params}

    def _build_dynamic_condition(self, field_info: Dict[str, Any], condition: models.AtomicCondition) -> Dict[str, Any]:
        """Build condition for dynamic property access"""
        # For dynamic properties, we need to filter by property name AND value
        # Handle numeric operators specially to avoid type mismatches
        value_column = field_info["sql_field"]
        if condition.operator in ["<", ">", "<=", ">=", "RANGE"]:
            # Cast the value column to numeric for comparison
            value_column = f"CAST({value_column} AS NUMERIC)"

        value_sql_expr, value_params = self.operators.get_sql_expression(
            condition.operator, value_column, condition.value, condition.threshold
        )

        if field_info["search_level"]["alias"] != field_info["subquery"]["alias"]:
            key = field_info["search_level"]["foreign_key"]
        else:
            key = "id"

        where = (
            f"EXISTS ( "
            f"{field_info['subquery']['sql']} "
            f"WHERE {field_info['subquery']['alias']}.{key}="
            f"{field_info['search_level']['alias']}{field_info['search_level']['alias']}.id "
            f"AND {field_info['property_filter']}"
            f"AND {value_sql_expr})  "
        )
        return {"sql": where, "params": value_params}

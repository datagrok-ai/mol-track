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

        if group_by != [] and f"{table_config['alias']}" not in group_by:
            group_by = [f"{table_config['alias']}.id"] + group_by

        group_by_sql = f"GROUP BY {' ,'.join(group_by)} " if group_by else " "

        # Build FROM clause with primary table
        base_from_clause = f"{schema}.{table_config['table']} {table_config['alias']}"

        # Collect all joins
        filter_query_joins = JoinOrderingTool()

        query_params = {}
        filter_sql = ""

        # Build sub query from filters
        if request.filter:
            sql_components = self.build_filter_sql_parts(request.filter, level, filter_query_joins)

            # Create query from the returned sql parts
            if sql_components["where"] or sql_components["having"]:
                from_clause = f"{schema}.{table_config['table']} {table_config['alias']}{table_config['alias']} "
                where_clause = f"WHERE {sql_components['where']}"
                having_clause = f"HAVING {sql_components['having']}" if sql_components["having"] else ""

                group_by_clause = f"GROUP BY {table_config['alias']}.id" if sql_components["having"] else ""

                query_params.update(sql_components["params"])
                filter_joins = filter_query_joins.getJoinSQL()
                filter_sql = (
                    "WHERE EXISTS ("
                    f"SELECT {table_config['alias']}.id "
                    f"FROM {from_clause} "
                    f"{filter_joins} "
                    f"{where_clause} "
                    "{connect} "
                    f"{group_by_clause} "
                    f"{having_clause}) "
                )
                # The inner query needs to have a different alias compared to the outter query
                filter_sql = filter_sql.replace(
                    f"{table_config['alias']}.", f"{table_config['alias']}{table_config['alias']}."
                )
                # This condition connects inner and outer table
                connect = f"{table_config['alias']}{table_config['alias']}.id={table_config['alias']}.id"
                connect = f" AND {connect}" if sql_components["where"] else f" {connect}"
                filter_sql = filter_sql.format(connect=connect)

        # Convert joins to list and remove duplicates
        base_joins = base_query_joins.getJoinSQL()

        # Build main query
        complete_sql = f"""
        SELECT {base_select_clause}
        FROM {base_from_clause}
        {base_joins}
        {filter_sql}
        {group_by_sql}
        ORDER BY {table_config["alias"]}.id
        """

        return {"sql": complete_sql.strip(), "params": query_params}

    def build_base_sql_parts(
        self, output_list: List[str], search_level: str, all_joins: JoinOrderingTool
    ) -> Dict[str, Any]:
        """
        Generate SQL query parts for base query.

        Returns:
            Dict with SQL components: select clause, where clause, group by
        """
        select_fields = []
        group_by = []
        # has_dynamic is True if any of the output fields are dynamic
        has_dynamic = False
        list_of_aliases = []
        conditions = []

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

        return {"select_clause": ", ".join(select_fields), "conditions": combined_sql, "group_by": group_by}

    def build_filter_sql_parts(
        self, filter_obj: models.Filter, level: str, all_joins: JoinOrderingTool
    ) -> Dict[str, Any]:
        """
        Builds SQL parts for nested filters recursively

        Returns:
            Dict with: {'where': str, 'having': str, 'params': Dict[str, Any]}
        """
        if isinstance(filter_obj, models.AtomicCondition):
            # Handle single atomic condition
            cond_info = self.build_condition(filter_obj, level, all_joins)
            return {
                "where": cond_info["where"] if cond_info["where"] else "",
                "having": cond_info["having"] if cond_info["having"] else "",
                "params": cond_info["params"],
            }

        elif isinstance(filter_obj, models.LogicalNode):
            # Handle logical node with multiple conditions
            where = []
            having = []
            all_params = {}

            for condition in filter_obj.conditions:
                if isinstance(condition, models.AtomicCondition):
                    cond_info = self.build_condition(condition, level, all_joins)
                    # Rename parameters to avoid conflicts
                    replaced_params = cond_info["where"] if cond_info["where"] else cond_info["having"]
                    renamed_params = {}
                    for old_name, value in cond_info["params"].items():
                        new_name = f"{old_name}_{self.parameter_counter}"
                        replaced_params = replaced_params.replace(f":{old_name}", f":{new_name}")
                        renamed_params[new_name] = value

                    if cond_info["where"]:
                        where.append(replaced_params)
                    else:
                        having.append(replaced_params)

                    all_params.update(renamed_params)
                    self.parameter_counter += 1

                elif isinstance(condition, models.LogicalNode):
                    # Recursive filter handling
                    nested_info = self.build_filter_sql_parts(condition, level, all_joins)
                    if nested_info["where"]:
                        where.append(f"({nested_info['where']})")
                    if nested_info["having"]:
                        having.append(f"({nested_info['having']})")
                    all_params.update(nested_info["params"])

            # Combine conditions with operator
            operator = f" {filter_obj.operator.value} "
            combined_where = operator.join(where)
            combined_having = operator.join(having)

            return {"where": combined_where, "having": combined_having, "params": all_params}

    def build_condition(
        self, condition: models.AtomicCondition, level: str, all_joins: JoinOrderingTool
    ) -> Dict[str, Any]:
        """
        Builds SQL parts for a single condition

        Returns:
            Dict with: {'where': str, 'having': str, 'params': Dict[str, Any]}
        """
        try:
            # Resolve field to SQL components
            field_info = self.field_resolver.resolve_field(condition.field, level, all_joins)

            # Validate operator and value
            self.operators.validate_operator_value(condition.operator, condition.value, condition.threshold)

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

        # Handle molecular searches (special case for structure field)
        if condition.field.endswith(".structure") and condition.operator in ["IS SIMILAR", "CONTAINS", "IS CONTAINED"]:
            return self._build_molecular_condition(field_sql, condition)

        # Standard condition building
        sql_expr, params = self.operators.get_sql_expression(
            condition.operator, field_sql, condition.value, condition.threshold
        )

        return {"where": sql_expr, "having": [], "params": params}

    def _build_dynamic_condition(self, field_info: Dict[str, Any], condition: models.AtomicCondition) -> Dict[str, Any]:
        """Build condition for dynamic property access"""
        # For dynamic properties, we need to filter by property name AND value
        # Handle numeric operators specially to avoid type mismatches
        value_column = field_info["sql_expression"]
        if condition.operator in ["<", ">", "<=", ">=", "RANGE"]:
            # Cast the value column to numeric for comparison
            value_column = f"CAST({value_column} AS NUMERIC)"

        value_sql_expr, value_params = self.operators.get_sql_expression(
            condition.operator, value_column, condition.value, condition.threshold
        )

        having = f"{value_sql_expr} "

        return {"where": [], "having": having, "params": value_params}

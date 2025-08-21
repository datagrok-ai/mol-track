"""
Field resolver for advanced search functionality
Handles resolution of field paths like 'compounds.details.chembl' to SQL components
"""

from typing import Any, Dict, get_args
from sqlalchemy.orm import Session
from app.services.search.utils.helper_functions import create_alias, singularize, get_table_columns
from app.services.search.utils.join_tools import JoinOrderingTool, JoinResolver
from app.models import Level


class FieldResolutionError(Exception):
    """Custom exception for field resolution errors"""

    pass


class FieldResolver:
    """Resolves field paths to SQL expressions and joins"""

    def __init__(self, db_schema: str, db: Session):
        self.db_schema = db_schema
        self._generate_table_config(db)
        self.join_resolver = JoinResolver(db_schema, self.table_configs)

    def _generate_table_config(self, db):
        tables = get_args(Level)
        self.table_configs = {}
        for table in tables:
            alias = create_alias(table)
            singular_name = singularize(table)
            self.table_configs[table] = {
                "table": table,
                "alias": alias,
                "details_table": f"{singular_name}_details",
                "details_alias": f"{alias}d",
                "details_fk": f"{singular_name}_id",
                "direct_fields": {column: f"{alias}.{column}" for column in get_table_columns(table, db)},
            }
        self.table_configs["compounds"]["direct_fields"]["structure"] = "c.canonical_smiles"

    def validate_field_path(self, field_path: str, level: Level = None) -> bool:
        """
        Validate that a field path is properly formatted
        Argument level is passed only when the output validation is being done
        """
        parts = field_path.split(".")

        table_name = parts[0]
        if table_name not in self.table_configs or (level is not None and table_name != level):
            return False

        if len(parts) == 2:
            # Direct field access
            field_name = parts[1]
            return field_name in self.table_configs[table_name]["direct_fields"]
        elif len(parts) == 3:
            # Dynamic property access
            return parts[1] == "details"

        return False

    def resolve_field(
        self, field_path: str, search_level: Level, all_joins: JoinOrderingTool, subquery: bool = False
    ) -> Dict[str, Any]:
        """
        Resolves a field path like 'compounds.details.chembl' to SQL components

        Args:
            field_path: Field path in format 'table.field' or 'table.details.property'
            search_level: Level  # The primary search level (compounds, batches, assay_results)
            all_joins: Join organization object
            subquery: True if the field is from filter. Indicates that EXISTS subquery needs to be created

        Returns:
            Dict with SQL components: {
                'sql_expression': str,     # SQL expression for the field
                'is_dynamic': bool,        # Whether this is a dynamic property
                'property_info': Dict,     # Property metadata if dynamic
                'table_alias': str,        # Table alias used
                'value_column': str        # Column containing the value
                'subquery' : {
                    'sql': str,            # SQL for subquery
                    'alias': str           # Alias
                }
                'search_level': {           # Search level information
                        'foreign_key': str,  # Foreign key
                        'alias' : str       # Alias
                    }
            }
        """
        search_level_info = {
            "search_level": {
                "foreign_key": self.table_configs[search_level]["details_fk"],
                "alias": self.table_configs[search_level]["alias"],
            }
        }
        parts = field_path.split(".")
        table_name = parts[0]
        field_or_details = parts[1]

        if table_name not in self.table_configs:
            raise FieldResolutionError(f"Unknown table: {table_name}")

        table_config = self.table_configs[table_name]

        cross_from = ""
        # Add cross-level join if needed
        if table_name != search_level:
            cross_joins, cross_tables, cross_from = self.join_resolver.resolve_join_components(
                search_level, table_name, subquery, field_or_details == "details"
            )
            if cross_joins:
                all_joins.add(cross_joins, cross_tables)

        # Handle direct field access
        if field_or_details != "details":
            if field_or_details in table_config["direct_fields"]:
                return (
                    self._resolve_direct_property(
                        table_config, field_or_details, search_level, all_joins, subquery, cross_from
                    )
                    | search_level_info
                )
            else:
                raise FieldResolutionError(f"Unknown direct field: {field_path}")

        property_name = parts[2]

        return (
            self._resolve_dynamic_property(table_config, property_name, all_joins, subquery, cross_from)
            | search_level_info
        )

    def _resolve_dynamic_property(
        self, table_config: Dict, property_name: str, joins: JoinOrderingTool, subquery: bool, cross_from: str
    ) -> Dict[str, Any]:
        alias = table_config["details_alias"]
        table = table_config["details_table"]
        base_table = table_config["table"]

        if subquery and cross_from == "":
            cross_from = table

        if table != cross_from:
            # Add details table join
            field = f"{table_config['alias']}.id"
            if subquery:
                if joins.joinCount() and not joins.checkLastJoin(base_table):
                    field = f"{joins.getLastTableAlias()}.{table_config['details_fk']}"
                elif not joins.joinCount() and base_table != cross_from:
                    field = f"{self.table_configs[cross_from]['alias']}.{table_config['details_fk']}"
            details_join = (
                f"LEFT JOIN {self.db_schema}.{table} {alias} ON {alias}.{table_config['details_fk']} = {field}"
            )
            joins.add([details_join], [table])

        # Add property join
        property_alias = f"p_{alias}"
        property_join = (
            f"LEFT JOIN {self.db_schema}.properties {property_alias} ON {property_alias}.id = {alias}.property_id"
        )
        joins.add([property_join], ["properties"])

        subquery_sql = ""
        subquery_alias = ""
        if subquery:
            joins_sql = joins.getJoinSQL()
            subquery_alias = self.table_configs[cross_from]["alias"] if table != cross_from else alias
            subquery_sql = f"SELECT 1 FROM {self.db_schema}.{cross_from} {subquery_alias} {joins_sql} "

        sql_expression = self.get_details_sql(table_config["table"], property_alias, alias)
        return {
            "sql_expression": sql_expression,
            "sql_field": sql_expression,
            "is_dynamic": True,
            "property_name": f"{property_alias}.name",
            "property_alias": f"{property_alias}_name",
            "table_alias": alias,
            "property_filter": f"LOWER({property_alias}_name) = LOWER('{property_name}')",
            "subquery": {
                "sql": subquery_sql,
                "alias": subquery_alias if subquery else "",
                "property_filter": f"LOWER({property_alias}.name) = LOWER('{property_name}')",
            },
        }

    def _resolve_direct_property(
        self,
        table_config: Dict,
        property_name: str,
        search_level: str,
        joins: JoinOrderingTool,
        subquery: bool,
        cross_from: str,
    ) -> Dict[str, Any]:
        subquery_sql = ""
        if subquery and search_level != table_config["table"]:
            joins_sql = joins.getJoinSQL()
            subquery_sql = (
                "SELECT 1 "
                f"FROM {self.db_schema}.{self.table_configs[cross_from]['table']} "
                f"{self.table_configs[cross_from]['alias']} "
                f"{joins_sql} "
            )

        search_level_alias = self.table_configs[search_level]["alias"]
        return {
            "sql_expression": table_config["direct_fields"][property_name].replace(
                f"{search_level_alias}.", f"{search_level_alias}{search_level_alias}."
            ),
            "sql_field": "",
            "is_dynamic": False,
            "table_alias": table_config["alias"],
            "property_filter": None,
            "subquery": {
                "sql": subquery_sql,
                "alias": self.table_configs[cross_from]["alias"] if subquery and cross_from != "" else "",
            },
        }

    def get_details_sql(self, table: str, property_alias, alias) -> str:
        """
        Get SQL for details table based on the main table
        """

        assay_parts = f"WHEN 'bool' THEN {alias}.value_bool::text " if table == "assay_results" else ""
        details_parts = (
            (f"WHEN 'datetime' THEN {alias}.value_datetime::text WHEN 'uuid' THEN {alias}.value_uuid::text ")
            if not table == "assay_results"
            else ""
        )
        sql_expression = (
            f"CASE {property_alias}.value_type "
            f"WHEN 'int' THEN {alias}.value_num::text "
            f"WHEN 'double' THEN {alias}.value_num::text "
            f"WHEN 'string' THEN {alias}.value_string "
            f"{assay_parts}"
            f"{details_parts}"
            f"END"
        )
        return sql_expression

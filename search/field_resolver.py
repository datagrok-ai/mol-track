"""
Field resolver for advanced search functionality
Handles resolution of field paths like 'compounds.details.chembl' to SQL components
"""

from typing import Dict, Any, List, Optional
import models
import enums
from .utils import sanitize_field_name


class FieldResolutionError(Exception):
    """Custom exception for field resolution errors"""
    pass


class FieldResolver:
    """Resolves field paths to SQL expressions and joins"""
    
    def __init__(self, db_schema: str):
        self.db_schema = db_schema
        self.table_configs = {
            "compounds": {
                "table": "compounds",
                "alias": "c",
                "details_table": "compound_details",
                "details_alias": "cd",
                "details_fk": "compound_id",
                "direct_fields": {
                    "id": "c.id",
                    "canonical_smiles": "c.canonical_smiles",
                    "inchi": "c.inchi", 
                    "inchikey": "c.inchikey",
                    "formula": "c.formula",
                    "molregno": "c.molregno",
                    "hash_mol": "c.hash_mol",
                    "hash_tautomer": "c.hash_tautomer",
                    "hash_canonical_smiles": "c.hash_canonical_smiles",
                    "hash_no_stereo_smiles": "c.hash_no_stereo_smiles",
                    "hash_no_stereo_tautomer": "c.hash_no_stereo_tautomer",
                    "created_at": "c.created_at",
                    "updated_at": "c.updated_at",
                    "is_archived": "c.is_archived",
                    "structure": "c.canonical_smiles"  # Alias for molecular searches
                }
            },
            "batches": {
                "table": "batches",
                "alias": "b",
                "details_table": "batch_details",
                "details_alias": "bd",
                "details_fk": "batch_id",
                "direct_fields": {
                    "id": "b.id",
                    "compound_id": "b.compound_id",
                    "batch_regno": "b.batch_regno",
                    "notes": "b.notes",
                    "created_at": "b.created_at",
                    "updated_at": "b.updated_at"
                }
            },
            "assay_results": {
                "table": "assay_results",
                "alias": "ar",
                "details_table": None,  # assay_results stores values directly
                "details_alias": None,
                "details_fk": None,
                "direct_fields": {
                    "id": "ar.id",
                    "batch_id": "ar.batch_id",
                    "assay_run_id": "ar.assay_run_id",
                    "property_id": "ar.property_id",
                    "value_qualifier": "ar.value_qualifier",
                    "value_num": "ar.value_num",
                    "value_string": "ar.value_string",
                    "value_bool": "ar.value_bool"
                }
            }
        }
        
        # Cross-level relationship definitions
        self.relationships = {
            ("compounds", "batches"): {
                "join": "INNER JOIN {schema}.batches b ON b.compound_id = c.id",
                "condition": "b.compound_id = c.id"
            },
            ("compounds", "assay_results"): {
                "join": "INNER JOIN {schema}.batches b ON b.compound_id = c.id INNER JOIN {schema}.assay_results ar ON ar.batch_id = b.id",
                "condition": "ar.batch_id = b.id AND b.compound_id = c.id"
            },
            ("batches", "assay_results"): {
                "join": "INNER JOIN {schema}.assay_results ar ON ar.batch_id = b.id",
                "condition": "ar.batch_id = b.id"
            },
            ("batches", "compounds"): {
                "join": "INNER JOIN {schema}.compounds c ON c.id = b.compound_id",
                "condition": "c.id = b.compound_id"
            },
            ("assay_results", "batches"): {
                "join": "INNER JOIN {schema}.batches b ON b.id = ar.batch_id",
                "condition": "b.id = ar.batch_id"
            },
            ("assay_results", "compounds"): {
                "join": "INNER JOIN {schema}.batches b ON b.id = ar.batch_id INNER JOIN {schema}.compounds c ON c.id = b.compound_id",
                "condition": "c.id = b.compound_id AND b.id = ar.batch_id"
            }
        }
    
    def resolve_field(self, field_path: str, search_level: str) -> Dict[str, Any]:
        """
        Resolves a field path like 'compounds.details.chembl' to SQL components
        
        Args:
            field_path: Field path in format 'table.field' or 'table.details.property'
            search_level: The primary search level (compounds, batches, assay_results)
            
        Returns:
            Dict with SQL components: {
                'sql_expression': str,     # SQL expression for the field
                'joins': List[str],        # Required JOIN clauses
                'is_dynamic': bool,        # Whether this is a dynamic property
                'property_info': Dict,     # Property metadata if dynamic
                'table_alias': str,        # Table alias used
                'value_column': str        # Column containing the value
            }
        """
        parts = field_path.split(".")
        if len(parts) < 2:
            raise FieldResolutionError(f"Invalid field path: {field_path}")
        
        table_name = parts[0]
        field_or_details = parts[1]
        
        if table_name not in self.table_configs:
            raise FieldResolutionError(f"Unknown table: {table_name}")
        
        table_config = self.table_configs[table_name]
        joins = []
        
        # Add cross-level join if needed
        if table_name != search_level:
            cross_join = self._get_cross_level_join(search_level, table_name)
            if cross_join:
                joins.extend(cross_join)
        
        # Handle direct field access
        if field_or_details != "details":
            if field_or_details in table_config["direct_fields"]:
                return {
                    'sql_expression': table_config["direct_fields"][field_or_details],
                    'joins': joins,
                    'is_dynamic': False,
                    'property_info': None,
                    'table_alias': table_config["alias"],
                    'value_column': table_config["direct_fields"][field_or_details]
                }
            else:
                raise FieldResolutionError(f"Unknown direct field: {field_path}")
        
        # Handle dynamic property access (table.details.property_name)
        if len(parts) != 3:
            raise FieldResolutionError(f"Dynamic property path must be 'table.details.property': {field_path}")
        
        property_name = parts[2]
        
        # Handle assay_results special case (no details table)
        if table_name == "assay_results":
            return self._resolve_assay_result_property(property_name, joins)
        
        # Handle regular details table lookup
        if not table_config["details_table"]:
            raise FieldResolutionError(f"Table {table_name} does not support dynamic properties")
        
        return self._resolve_dynamic_property(table_config, property_name, joins)
    
    def _get_cross_level_join(self, from_level: str, to_level: str) -> List[str]:
        """Get JOIN clauses for cross-level relationships"""
        relationship_key = (from_level, to_level)
        if relationship_key in self.relationships:
            join_sql = self.relationships[relationship_key]["join"].format(schema=self.db_schema)
            return [join_sql]
        else:
            raise FieldResolutionError(f"No relationship defined from {from_level} to {to_level}")
    
    def _resolve_assay_result_property(self, property_name: str, joins: List[str]) -> Dict[str, Any]:
        """Handle assay_results property resolution (values stored directly)"""
        # Add property join to get property info
        property_join = f"INNER JOIN {self.db_schema}.properties p_ar ON p_ar.id = ar.property_id"
        joins.append(property_join)

        # For output, use CASE for display
        sql_expression = (
            "CASE p_ar.value_type "
            "WHEN 'int' THEN ar.value_num::text "
            "WHEN 'double' THEN ar.value_num::text "
            "WHEN 'string' THEN ar.value_string "
            "WHEN 'bool' THEN ar.value_bool::text "
            "END"
        )

        # For filtering, return the correct column name and value_type for the property
        # The query builder will use this to select the right column
        return {
            'sql_expression': sql_expression,
            'joins': joins,
            'is_dynamic': True,
            'property_info': {'name': property_name, 'table': 'assay_results'},
            'table_alias': 'ar',
            'value_column': None,  # Will be set in the query builder
            'property_value_type': f"(SELECT value_type FROM {self.db_schema}.properties WHERE name = '{property_name}' LIMIT 1)",
            'property_filter': f"p_ar.name = '{property_name}'"
        }
    
    def _resolve_dynamic_property(self, table_config: Dict, property_name: str, joins: List[str]) -> Dict[str, Any]:
        """Resolve dynamic property from details table"""
        details_alias = table_config["details_alias"]
        details_table = table_config["details_table"]
        
        # Add details table join
        details_join = (f"INNER JOIN {self.db_schema}.{details_table} {details_alias} "
                       f"ON {details_alias}.{table_config['details_fk']} = {table_config['alias']}.id")
        joins.append(details_join)
        
        # Add property join
        property_alias = f"p_{details_alias}"
        property_join = f"INNER JOIN {self.db_schema}.properties {property_alias} ON {property_alias}.id = {details_alias}.property_id"
        joins.append(property_join)
        
        return {
            'sql_expression': f"CASE {property_alias}.value_type "
                            f"WHEN 'int' THEN {details_alias}.value_num::text "
                            f"WHEN 'double' THEN {details_alias}.value_num::text "
                            f"WHEN 'string' THEN {details_alias}.value_string "
                            f"WHEN 'datetime' THEN {details_alias}.value_datetime::text "
                            f"WHEN 'uuid' THEN {details_alias}.value_uuid::text "
                            f"END",
            'joins': joins,
            'is_dynamic': True,
            'property_info': {'name': property_name, 'table': details_table},
            'table_alias': details_alias,
            'value_column': 'dynamic',
            'property_filter': f"{property_alias}.name = '{property_name}'"
        }
    
    def get_output_fields(self, output_list: List[str], search_level: str) -> Dict[str, Any]:
        """
        Resolve all output fields for SELECT clause
        
        Returns:
            Dict with SELECT components and required JOINs
        """
        select_fields = []
        all_joins = []
        
        for field_path in output_list:
            resolved = self.resolve_field(field_path, search_level)
            
            # Create alias for the field
            field_alias = sanitize_field_name(field_path)
            select_fields.append(f"{resolved['sql_expression']} AS {field_alias}")
            
            # Collect joins
            all_joins.extend(resolved['joins'])
        
        # Remove duplicate joins
        unique_joins = list(dict.fromkeys(all_joins))
        
        print(f"DEBUG: select_fields = {select_fields}")
        return {
            'select_clause': ", ".join(select_fields),
            'joins': unique_joins
        }
    
    def validate_field_path(self, field_path: str) -> bool:
        """Validate that a field path is properly formatted"""
        parts = field_path.split(".")
        
        if len(parts) < 2:
            return False
        
        table_name = parts[0]
        if table_name not in self.table_configs:
            return False
        
        if len(parts) == 2:
            # Direct field access
            field_name = parts[1]
            return field_name in self.table_configs[table_name]["direct_fields"] or field_name == "details"
        elif len(parts) == 3:
            # Dynamic property access
            return parts[1] == "details"
        
        return False 
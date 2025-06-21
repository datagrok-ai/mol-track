"""
Query builder for advanced search functionality
Builds dynamic SQL queries from search requests
"""

from typing import Dict, Any, List, Tuple, Union
from sqlalchemy import text
import models
from .field_resolver import FieldResolver, FieldResolutionError
from .operators import SearchOperators


class QueryBuildError(Exception):
    """Custom exception for query building errors"""
    pass


class QueryBuilder:
    """Builds dynamic SQL queries from search requests"""
    
    def __init__(self, field_resolver: FieldResolver):
        self.field_resolver = field_resolver
        self.operators = SearchOperators()
    
    def build_query(self, request: models.SearchRequest) -> Dict[str, Any]:
        """
        Builds complete SQL query from search request
        
        Returns:
            Dict with query components: {
                'sql': str,              # Main query SQL
                'params': Dict[str, Any] # Query parameters
            }
        """
        level = request.level
        schema = self.field_resolver.db_schema
        table_config = self.field_resolver.table_configs[level]
        
        # Build SELECT clause
        output_info = self.field_resolver.get_output_fields(request.output, level)
        select_clause = output_info['select_clause']
        
        # Build FROM clause with primary table
        from_clause = f"{schema}.{table_config['table']} {table_config['alias']}"
        
        # Collect all joins
        all_joins = set(output_info['joins'])
        where_clause = ""
        query_params = {}
        
        # Build WHERE clause from filters
        if request.filter:
            where_info = self.build_filter(request.filter, level)
            where_clause = f"WHERE {where_info['sql']}"
            query_params.update(where_info['params'])
            all_joins.update(where_info['joins'])
        
        # Convert joins to list and remove duplicates
        join_clauses = list(all_joins)
        joins_sql = " ".join(join_clauses) if join_clauses else ""
        
        # Build main query
        main_sql = f"""
        SELECT {select_clause}
        FROM {from_clause}
        {joins_sql}
        {where_clause}
        ORDER BY {table_config['alias']}.id
        """
        
        return {
            'sql': main_sql.strip(),
            'params': query_params
        }
    
    def build_filter(self, filter_obj: models.Filter, level: str) -> Dict[str, Any]:
        """
        Builds WHERE clause for nested filters
        
        Returns:
            Dict with: {'sql': str, 'params': Dict[str, Any], 'joins': Set[str]}
        """
        if isinstance(filter_obj, models.AtomicCondition):
            # Handle single atomic condition
            cond_info = self.build_condition(filter_obj, level)
            return {
                'sql': cond_info['sql'],
                'params': cond_info['params'],
                'joins': cond_info['joins']
            }
        
        elif isinstance(filter_obj, models.LogicalNode):
            # Handle logical node with multiple conditions
            conditions = []
            all_params = {}
            all_joins = set()
            param_counter = 0
            
            for condition in filter_obj.conditions:
                if isinstance(condition, models.AtomicCondition):
                    cond_info = self.build_condition(condition, level)
                    # Rename parameters to avoid conflicts
                    renamed_sql = cond_info['sql']
                    renamed_params = {}
                    for old_name, value in cond_info['params'].items():
                        new_name = f"{old_name}_{param_counter}"
                        renamed_sql = renamed_sql.replace(f":{old_name}", f":{new_name}")
                        renamed_params[new_name] = value
                    
                    conditions.append(renamed_sql)
                    all_params.update(renamed_params)
                    all_joins.update(cond_info['joins'])
                    param_counter += 1
                
                elif isinstance(condition, models.LogicalNode):
                    # Recursive filter handling
                    nested_info = self.build_filter(condition, level)
                    conditions.append(f"({nested_info['sql']})")
                    all_params.update(nested_info['params'])
                    all_joins.update(nested_info['joins'])
            
            # Combine conditions with operator
            operator = f" {filter_obj.operator.value} "
            combined_sql = operator.join(conditions)
            
            return {
                'sql': combined_sql,
                'params': all_params,
                'joins': all_joins
            }
    
    def build_condition(self, condition: models.AtomicCondition, level: str) -> Dict[str, Any]:
        """
        Builds WHERE clause for a single condition
        
        Returns:
            Dict with: {'sql': str, 'params': Dict[str, Any], 'joins': Set[str]}
        """
        try:
            # Resolve field to SQL components
            field_info = self.field_resolver.resolve_field(condition.field, level)
            
            # Validate operator and value
            self.operators.validate_operator_value(condition.operator, condition.value, condition.threshold)
            
            # Handle dynamic properties with property name filtering
            if field_info['is_dynamic']:
                return self._build_dynamic_condition(field_info, condition)
            else:
                return self._build_direct_condition(field_info, condition)
                
        except FieldResolutionError as e:
            raise QueryBuildError(f"Field resolution error: {str(e)}")
        except ValueError as e:
            raise QueryBuildError(f"Condition validation error: {str(e)}")
    
    def _build_direct_condition(self, field_info: Dict[str, Any], condition: models.AtomicCondition) -> Dict[str, Any]:
        """Build condition for direct field access"""
        field_sql = field_info['sql_expression']
        
        # Handle molecular searches (special case for structure field)
        if condition.field.endswith('.structure') and condition.operator in ['IS SIMILAR', 'CONTAINS', 'IS CONTAINED']:
            return self._build_molecular_condition(field_sql, condition)
        
        # Standard condition building
        sql_expr, params = self.operators.get_sql_expression(
            condition.operator, 
            field_sql, 
            condition.value, 
            condition.threshold
        )
        
        return {
            'sql': sql_expr,
            'params': params,
            'joins': set(field_info['joins'])
        }
    
    def _build_dynamic_condition(self, field_info: Dict[str, Any], condition: models.AtomicCondition) -> Dict[str, Any]:
        """Build condition for dynamic property access"""
        property_filter = field_info['property_filter']
        
        # For dynamic properties, we need to filter by property name AND value
        # Handle numeric operators specially to avoid type mismatches
        if condition.operator in ['<', '>', '<=', '>=', 'RANGE']:
            # For numeric operators, cast the text value back to numeric
            value_column = (
                "CASE p_ar.value_type "
                "WHEN 'int' THEN ar.value_num::text "
                "WHEN 'double' THEN ar.value_num::text "
                "WHEN 'string' THEN ar.value_string "
                "WHEN 'bool' THEN ar.value_bool::text "
                "END"
            )
            # Cast the value column to numeric for comparison
            value_column = f"CAST({value_column} AS NUMERIC)"
        else:
            # For non-numeric operators, use text comparison
            value_column = (
                "CASE p_ar.value_type "
                "WHEN 'int' THEN ar.value_num::text "
                "WHEN 'double' THEN ar.value_num::text "
                "WHEN 'string' THEN ar.value_string "
                "WHEN 'bool' THEN ar.value_bool::text "
                "END"
            )

        value_sql_expr, value_params = self.operators.get_sql_expression(
            condition.operator,
            value_column,
            condition.value,
            condition.threshold
        )

        combined_sql = f"({property_filter} AND {value_sql_expr})"

        return {
            'sql': combined_sql,
            'params': value_params,
            'joins': set(field_info['joins'])
        }
    
    def _build_molecular_condition(self, field_sql: str, condition: models.AtomicCondition) -> Dict[str, Any]:
        """Build condition for molecular searches using RDKit"""
        if condition.operator == 'IS SIMILAR':
            if condition.threshold is None:
                raise QueryBuildError("Molecular similarity search requires a threshold value")
            
            # Use RDKit Tanimoto similarity
            sql = f"tanimoto_sml(mol_from_smiles({field_sql}), mol_from_smiles(:param1)) >= :param2"
            params = {"param1": condition.value, "param2": condition.threshold}
        
        elif condition.operator == 'CONTAINS':
            # Substructure search
            sql = f"mol_from_smiles({field_sql}) @> mol_from_smiles(:param)"
            params = {"param": condition.value}
        
        elif condition.operator == 'IS CONTAINED':
            # Superstructure search  
            sql = f"mol_from_smiles({field_sql}) <@ mol_from_smiles(:param)"
            params = {"param": condition.value}
        
        else:
            raise QueryBuildError(f"Unsupported molecular operator: {condition.operator}")
        
        return {
            'sql': sql,
            'params': params,
            'joins': set()
        }
    
    def validate_cross_level_constraints(self, request: models.SearchRequest) -> List[str]:
        """
        Validate that conditions within same filter group don't mix levels inappropriately
        
        Returns:
            List of validation errors
        """
        errors = []
        
        if not request.filter:
            return errors
        
        def check_filter_levels(filter_obj: models.Filter, path: str = "root") -> None:
            if isinstance(filter_obj, models.AtomicCondition):
                # Single atomic condition - no cross-level issues
                return
            elif isinstance(filter_obj, models.LogicalNode):
                condition_levels = set()
                
                for i, condition in enumerate(filter_obj.conditions):
                    if isinstance(condition, models.AtomicCondition):
                        # Extract table level from field
                        field_level = condition.field.split('.')[0]
                        condition_levels.add(field_level)
                    elif isinstance(condition, models.LogicalNode):
                        # Recursively check nested filters
                        check_filter_levels(condition, f"{path}.conditions[{i}]")
                
                # Check if we have mixed levels in same filter group
                if len(condition_levels) > 1:
                    # This is actually allowed - it creates cross-level joins
                    # But we might want to warn about potential performance impact
                    pass
        
        try:
            check_filter_levels(request.filter)
        except Exception as e:
            errors.append(f"Filter validation error: {str(e)}")
        
        return errors 
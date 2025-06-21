"""
Main search engine for advanced search functionality
Orchestrates field resolution, query building, and execution
"""

from typing import Dict, Any, List
from sqlalchemy.orm import Session
from sqlalchemy import text
import models
from .field_resolver import FieldResolver, FieldResolutionError
from .query_builder import QueryBuilder, QueryBuildError


class SearchEngineError(Exception):
    """Custom exception for search engine errors"""
    pass


class SearchEngine:
    """Main search orchestration"""
    
    def __init__(self, db: Session, db_schema: str = None):
        self.db = db
        self.db_schema = db_schema or models.DB_SCHEMA
        self.field_resolver = FieldResolver(self.db_schema)
        self.query_builder = QueryBuilder(self.field_resolver)
    
    def search(self, request: models.SearchRequest) -> models.SearchResponse:
        """
        Executes search request
        
        Args:
            request: SearchRequest object with search parameters
            
        Returns:
            SearchResponse with results and metadata
        """
        try:
            # Validate the request
            validation_errors = self.validate_request(request)
            if validation_errors:
                raise SearchEngineError(f"Request validation failed: {'; '.join(validation_errors)}")
            
            # Build the SQL query
            query_info = self.query_builder.build_query(request)
            
            # Execute main query
            results = self._execute_main_query(query_info['sql'], query_info['params'])
            
            # Extract column names from output fields
            columns = [field.replace(".", "_") for field in request.output]
            
            return models.SearchResponse(
                data=results,
                total_count=len(results),
                level=request.level,
                columns=columns
            )
            
        except (FieldResolutionError, QueryBuildError) as e:
            raise SearchEngineError(f"Search execution error: {str(e)}")
        except Exception as e:
            raise SearchEngineError(f"Unexpected error during search: {str(e)}")
    
    def validate_request(self, request: models.SearchRequest) -> List[str]:
        """
        Validates search request constraints
        
        Returns:
            List of validation error messages
        """
        errors = []
        
        # Validate output fields
        for field_path in request.output:
            if not self.field_resolver.validate_field_path(field_path):
                errors.append(f"Invalid output field: {field_path}")
        
        # Validate filter conditions if present
        if request.filter:
            filter_errors = self._validate_filter(request.filter, request.level)
            errors.extend(filter_errors)
        
        # Validate cross-level constraints
        cross_level_errors = self.query_builder.validate_cross_level_constraints(request)
        errors.extend(cross_level_errors)
        
        return errors
    
    def _validate_filter(self, filter_obj: models.Filter, level: str, path: str = "filter") -> List[str]:
        """Recursively validate filter conditions"""
        errors = []
        
        if isinstance(filter_obj, models.AtomicCondition):
            # Validate single atomic condition
            if not self.field_resolver.validate_field_path(filter_obj.field):
                errors.append(f"Invalid field at {path}: {filter_obj.field}")
            
        elif isinstance(filter_obj, models.LogicalNode):
            # Validate logical node with multiple conditions
            for i, condition in enumerate(filter_obj.conditions):
                condition_path = f"{path}.conditions[{i}]"
                
                if isinstance(condition, models.AtomicCondition):
                    # Validate field path
                    if not self.field_resolver.validate_field_path(condition.field):
                        errors.append(f"Invalid field at {condition_path}: {condition.field}")
                    
                elif isinstance(condition, models.LogicalNode):
                    # Recursively validate nested filters
                    nested_errors = self._validate_filter(condition, level, condition_path)
                    errors.extend(nested_errors)
        
        return errors
    
    def _execute_main_query(self, sql: str, params: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Execute main query and return results as list of dictionaries"""
        try:
            # SQLAlchemy expects parameters to be passed as keyword arguments
            result = self.db.execute(text(sql), params)
            
            # Get column names from result
            if result.returns_rows:
                columns = list(result.keys())
                rows = result.fetchall()
                
                # Convert rows to dictionaries
                return [dict(zip(columns, row)) for row in rows]
            else:
                return []
                
        except Exception as e:
            raise SearchEngineError(f"Main query execution failed: {str(e)}")
    
    def explain_query(self, request: models.SearchRequest) -> Dict[str, Any]:
        """
        Generate query explanation without executing
        Useful for debugging and optimization
        
        Returns:
            Dict with query details and execution plan
        """
        try:
            # Validate the request
            validation_errors = self.validate_request(request)
            
            # Build the SQL query
            query_info = self.query_builder.build_query(request)
            
            # Get execution plan
            explain_sql = f"EXPLAIN ANALYZE {query_info['sql']}"
            
            execution_plan = None
            if self.db:
                try:
                    result = self.db.execute(text(explain_sql), query_info['params'])
                    execution_plan = [row[0] for row in result.fetchall()]
                except Exception as e:
                    execution_plan = f"Could not get execution plan: {str(e)}"
            
            return {
                'validation_errors': validation_errors,
                'main_sql': query_info['sql'],
                'parameters': query_info['params'],
                'execution_plan': execution_plan,
                'is_valid': len(validation_errors) == 0
            }
            
        except Exception as e:
            return {
                'error': str(e),
                'is_valid': False
            }
    
    def get_available_fields(self, level: str) -> Dict[str, List[str]]:
        """
        Get list of available fields for a given search level
        
        Returns:
            Dict with direct fields and example dynamic properties
        """
        if level not in self.field_resolver.table_configs:
            raise SearchEngineError(f"Unknown search level: {level}")
        
        table_config = self.field_resolver.table_configs[level]
        
        # Direct fields
        direct_fields = list(table_config['direct_fields'].keys())
        
        # Example dynamic properties (would need to query database for actual list)
        dynamic_examples = [
            f"{level}.details.example_property",
            f"{level}.details.another_property"
        ]
        
        # Cross-level fields available
        cross_level_fields = []
        for other_level in self.field_resolver.table_configs.keys():
            if other_level != level:
                other_config = self.field_resolver.table_configs[other_level]
                for field_name in other_config['direct_fields'].keys():
                    cross_level_fields.append(f"{other_level}.{field_name}")
                cross_level_fields.append(f"{other_level}.details.*")
        
        return {
            'direct_fields': direct_fields,
            'dynamic_properties': dynamic_examples,
            'cross_level_fields': cross_level_fields
        }
    
    def get_supported_operators(self) -> Dict[str, Any]:
        """Get list of supported operators with descriptions"""
        from .operators import SearchOperators
        
        operators_info = {}
        for op_name, op_def in SearchOperators.OPERATORS.items():
            operators_info[op_name] = {
                'description': op_def['description'],
                'type': op_def['type'].value,
                'requires_threshold': op_def.get('requires_threshold', False)
            }
        
        return operators_info 
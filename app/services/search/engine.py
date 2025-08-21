from typing import Any, Dict, List
from sqlalchemy.orm import Session
import app.models as models
from app.services.search.field_resolver import FieldResolver, FieldResolutionError
from app.services.search.query_builder import QueryBuilder, QueryBuildError
from sqlalchemy import text
from app.setup.database import DB_SCHEMA
from app.services.search.operators import SearchOperators
from app.services.search.utils.helper_functions import create_alias_mapping, prepare_search_output
from app.models import Level


class SearchEngineError(Exception):
    """Custom exception for search engine errors"""

    pass


class SearchEngine:
    """Search orchestration engine"""

    def __init__(self, db: Session):
        self.db = db
        self.db_schema = DB_SCHEMA
        self.field_resolver = FieldResolver(self.db_schema, db)
        self.query_builder = QueryBuilder(self.field_resolver)
        self.results = None

    def search(self, request: models.SearchRequest):
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

            # Prepare output fields
            self.prepare_output_fields(request)

            # Build the SQL query
            query_info = self.query_builder.build_query(request, self.output_aliases)

            # Execute main query
            results, headers = self._execute_main_query(query_info["sql"], query_info["params"])

            return prepare_search_output(results, headers, request.output_format)

        except (FieldResolutionError, QueryBuildError) as e:
            raise SearchEngineError(f"Search execution error: {str(e)}")
        except Exception as e:
            raise SearchEngineError(f"Error: {str(e)}")

    def validate_request(self, request: models.SearchRequest) -> List[str]:
        """
        Validates search request constraints

        Returns:
            List of validation error messages
        """
        errors = []

        # Validate output fields
        for field_path in request.output:
            if not self.field_resolver.validate_field_path(field_path, request.level):
                errors.append(f"Invalid output field: {field_path}")

        # Validate output fields
        for aggregation in request.aggregations:
            field_path = aggregation.field
            if not self.field_resolver.validate_field_path(field_path):
                errors.append(f"Invalid output field: {field_path}")

        # Validate filter conditions if present
        if request.filter:
            filter_errors = self._validate_filter(request.filter, request.level)
            errors.extend(filter_errors)

        return errors

    def prepare_output_fields(self, request: models.SearchRequest):
        """
        Prepares output fields by resolving aliases and ensuring valid field paths
        """
        if not request.output:
            raise SearchEngineError("Output fields cannot be empty")

        columns = [field.lower() for field in request.output]
        self.output_aliases = create_alias_mapping(columns, request.aggregations)

    def _validate_filter(self, filter_obj: models.Filter, level: Level, path: str = "filter") -> List[str]:
        """Recursively validate filter conditions"""
        errors = []

        if isinstance(filter_obj, models.AtomicCondition):
            # Validate single atomic condition
            if not self.field_resolver.validate_field_path(filter_obj.field):
                errors.append(f"Invalid field at {path}: {filter_obj.field}")
            SearchOperators.validate_operands(filter_obj.operator, filter_obj.field)
            SearchOperators.validate_operator_value(filter_obj.operator, filter_obj.value, filter_obj.threshold)

        elif isinstance(filter_obj, models.LogicalNode):
            # Validate logical node with multiple conditions
            for i, condition in enumerate(filter_obj.conditions):
                condition_path = f"{path}.conditions[{i}]"

                if isinstance(condition, models.AtomicCondition):
                    # Validate field path
                    if not self.field_resolver.validate_field_path(condition.field):
                        errors.append(f"Invalid field at {condition_path}: {condition.field}")
                    SearchOperators.validate_operands(condition.operator, condition.field)
                    SearchOperators.validate_operator_value(condition.operator, condition.value, condition.threshold)

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
                headers = [self.output_aliases[item][0] for item in list(result.keys())]
                rows = result.fetchall()
                return rows, headers
            else:
                return [], []

        except Exception as e:
            raise SearchEngineError(f"Main query execution failed: {str(e)}")

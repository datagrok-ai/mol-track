"""
Advanced Search Module for MolTrack

This module provides comprehensive search functionality across compounds, batches, and assay results
with support for dynamic filtering, cross-level queries, and molecular structure searches.
"""

from .engine import SearchEngine, SearchEngineError
from .field_resolver import FieldResolver, FieldResolutionError
from .query_builder import QueryBuilder, QueryBuildError
from .operators import SearchOperators, OperatorType
from .validators import SearchValidator
from .utils import (
    sanitize_field_name, 
    build_field_alias, 
    validate_smiles_string,
    format_sql_query,
    estimate_query_complexity
)

__all__ = [
    "SearchEngine", 
    "SearchEngineError",
    "FieldResolver", 
    "FieldResolutionError",
    "QueryBuilder",
    "QueryBuildError", 
    "SearchOperators",
    "OperatorType",
    "SearchValidator",
    "sanitize_field_name",
    "build_field_alias",
    "validate_smiles_string", 
    "format_sql_query",
    "estimate_query_complexity"
] 
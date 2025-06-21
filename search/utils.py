"""
Utility functions for search functionality
"""

from typing import Dict, Any, List, Set
import re


def sanitize_field_name(field_name: str) -> str:
    """Sanitize field name for use in SQL aliases"""
    # Replace dots and special characters with underscores
    sanitized = re.sub(r'[^a-zA-Z0-9_]', '_', field_name)
    # Ensure it doesn't start with a number
    if sanitized and sanitized[0].isdigit():
        sanitized = f"field_{sanitized}"
    return sanitized


def extract_unique_joins(joins_list: List[str]) -> List[str]:
    """Extract unique JOIN clauses, preserving order"""
    seen = set()
    unique_joins = []
    
    for join in joins_list:
        # Normalize whitespace for comparison
        normalized = ' '.join(join.split())
        if normalized not in seen:
            seen.add(normalized)
            unique_joins.append(join)
    
    return unique_joins


def build_field_alias(field_path: str) -> str:
    """Build a safe SQL alias from a field path"""
    return sanitize_field_name(field_path.replace(".", "_"))


def validate_smiles_string(smiles: str) -> bool:
    """Validate a SMILES string using RDKit"""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except Exception:
        return False


def format_sql_query(sql: str) -> str:
    """Format SQL query for better readability"""
    # Basic SQL formatting - could be enhanced with a proper SQL formatter
    formatted = sql.strip()
    
    # Add line breaks after major clauses
    clauses = ['SELECT', 'FROM', 'WHERE', 'JOIN', 'ORDER BY', 'LIMIT', 'OFFSET']
    for clause in clauses:
        formatted = formatted.replace(f' {clause} ', f'\n{clause} ')
    
    return formatted


def extract_table_references(field_paths: List[str]) -> Set[str]:
    """Extract unique table references from field paths"""
    tables = set()
    for field_path in field_paths:
        parts = field_path.split('.')
        if parts:
            tables.add(parts[0])
    return tables


def merge_query_parameters(param_lists: List[List[Any]]) -> List[Any]:
    """Merge multiple parameter lists into one"""
    merged = []
    for param_list in param_lists:
        merged.extend(param_list)
    return merged


def estimate_query_complexity(request) -> Dict[str, Any]:
    """Estimate query complexity for performance warnings"""
    complexity = {
        'score': 0,
        'factors': []
    }
    
    # Count output fields
    output_count = len(request.output)
    if output_count > 10:
        complexity['score'] += 2
        complexity['factors'].append(f"Many output fields ({output_count})")
    
    # Count filter conditions
    if request.filter:
        condition_count = count_filter_conditions(request.filter)
        if condition_count > 5:
            complexity['score'] += condition_count // 5
            complexity['factors'].append(f"Many filter conditions ({condition_count})")
    
    # Check for cross-level queries
    table_refs = extract_table_references(request.output)
    if request.filter:
        filter_refs = extract_filter_table_references(request.filter)
        table_refs.update(filter_refs)
    
    if len(table_refs) > 1:
        complexity['score'] += len(table_refs)
        complexity['factors'].append(f"Cross-level query ({len(table_refs)} tables)")
    
    return complexity


def count_filter_conditions(filter_obj) -> int:
    """Recursively count filter conditions"""
    count = 0
    
    for condition in filter_obj.conditions:
        if hasattr(condition, 'field'):  # AtomicCondition
            count += 1
        elif hasattr(condition, 'conditions'):  # Nested LogicalNode
            count += count_filter_conditions(condition)
    
    return count


def extract_filter_table_references(filter_obj) -> Set[str]:
    """Extract table references from filter conditions"""
    tables = set()
    
    for condition in filter_obj.conditions:
        if hasattr(condition, 'field'):  # AtomicCondition
            parts = condition.field.split('.')
            if parts:
                tables.add(parts[0])
        elif hasattr(condition, 'conditions'):  # Nested LogicalNode
            nested_refs = extract_filter_table_references(condition)
            tables.update(nested_refs)
    
    return tables 
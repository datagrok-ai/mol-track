"""
Search operators and their SQL translations
"""
from typing import Dict, Any, List, Tuple
from enum import Enum


class OperatorType(Enum):
    """Types of operators for different data types"""
    STRING = "string"
    NUMERIC = "numeric"
    DATETIME = "datetime"
    BOOLEAN = "boolean"
    MOLECULAR = "molecular"


class SearchOperators:
    """Maps search operators to SQL expressions and validation"""
    
    #TODO: CHECK OPERATORS FOR VALIDITY AND CONSISTENCY
    # Operator definitions with their SQL translations and parameter handling
    OPERATORS = {
        # String operators
        "=": {
            "sql": "= :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Exact match"
        },
        "!=": {
            "sql": "!= :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Not equal"
        },
        "LIKE": {
            "sql": "LIKE :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Pattern match with wildcards"
        },
        "CONTAINS": {
            "sql": "ILIKE :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Case-insensitive contains",
            "value_transform": lambda x: f"%{x}%"
        },
        "STARTS WITH": {
            "sql": "ILIKE :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Starts with pattern",
            "value_transform": lambda x: f"{x}%"
        },
        "ENDS WITH": {
            "sql": "ILIKE :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Ends with pattern",
            "value_transform": lambda x: f"%{x}"
        },
        "IN": {
            "sql": "IN :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Match any value in list",
            "value_transform": lambda x: tuple(x) if isinstance(x, list) else x
        },
        
        # Numeric operators
        "<": {
            "sql": "< :param",
            "type": OperatorType.NUMERIC,
            "params": 1,
            "description": "Less than"
        },
        ">": {
            "sql": "> :param",
            "type": OperatorType.NUMERIC,
            "params": 1,
            "description": "Greater than"
        },
        "<=": {
            "sql": "<= :param",
            "type": OperatorType.NUMERIC,
            "params": 1,
            "description": "Less than or equal"
        },
        ">=": {
            "sql": ">= :param",
            "type": OperatorType.NUMERIC,
            "params": 1,
            "description": "Greater than or equal"
        },
        "RANGE": {
            "sql": "BETWEEN :param1 AND :param2",
            "type": OperatorType.NUMERIC,
            "params": 2,
            "description": "Between two values (inclusive)",
            "value_transform": lambda x: (x[0], x[1]) if isinstance(x, (list, tuple)) and len(x) == 2 else x
        },
        
        # Datetime operators
        "BEFORE": {
            "sql": "< :param",
            "type": OperatorType.DATETIME,
            "params": 1,
            "description": "Before date/time"
        },
        "AFTER": {
            "sql": "> :param",
            "type": OperatorType.DATETIME,
            "params": 1,
            "description": "After date/time"
        },
        "ON": {
            "sql": "DATE(:param) = DATE(:param)",
            "type": OperatorType.DATETIME,
            "params": 1,
            "description": "On specific date"
        },
        
        # Molecular operators (RDKit)
        "IS SIMILAR": {
            "sql": "mol_from_smiles(:param) %% mol_from_smiles(:param)",
            "type": OperatorType.MOLECULAR,
            "params": 1,
            "description": "Molecular similarity using RDKit",
            "requires_threshold": True
        },
        "CONTAINS": {
            "sql": "mol_from_smiles(:param) @> mol_from_smiles(:param)",
            "type": OperatorType.MOLECULAR,
            "params": 1,
            "description": "Molecular substructure search"
        },
        "IS CONTAINED": {
            "sql": "mol_from_smiles(:param) <@ mol_from_smiles(:param)",
            "type": OperatorType.MOLECULAR,
            "params": 1,
            "description": "Molecular superstructure search"
        }
    }

    @classmethod
    def get_operator(cls, operator: str) -> Dict[str, Any]:
        """Get operator definition"""
        if operator not in cls.OPERATORS:
            raise ValueError(f"Unsupported operator: {operator}")
        return cls.OPERATORS[operator]
    
    @classmethod
    def validate_operator_value(cls, operator: str, value: Any, threshold: float = None) -> bool:
        """Validate that a value is appropriate for the given operator"""
        op_def = cls.get_operator(operator)
        
        # Check if threshold is required
        if op_def.get("requires_threshold", False) and threshold is None:
            raise ValueError(f"Operator '{operator}' requires a threshold value")
        
        # Check value type based on operator
        if operator in ["IN"] and not isinstance(value, (list, tuple)):
            raise ValueError(f"Operator '{operator}' requires a list or tuple value")
        
        if operator == "RANGE":
            if not isinstance(value, (list, tuple)) or len(value) != 2:
                raise ValueError(f"Operator '{operator}' requires a list/tuple with exactly 2 values")
            if value[0] >= value[1]:
                raise ValueError(f"RANGE operator requires first value to be less than second value")
        
        return True
    
    @classmethod
    def get_sql_expression(cls, operator: str, field: str, value: Any, threshold: float = None) -> Tuple[str, Dict[str, Any]]:
        """
        Get SQL expression and parameters for an operator
        
        Returns:
            Tuple of (sql_expression, parameters_dict)
        """
        #TODO: Make it use the dict defined at the top.....
        op_def = cls.get_operator(operator)
        
        # Transform value if needed
        if "value_transform" in op_def:
            value = op_def["value_transform"](value)
        
        # Handle molecular similarity with threshold
        if operator == "IS SIMILAR" and threshold is not None:
            sql_expr = f"mol_from_smiles({field}) %% mol_from_smiles(:param1) AND tanimoto_sml(mol_from_smiles({field}), mol_from_smiles(:param2)) >= :param3"
            return sql_expr, {"param1": value, "param2": value, "param3": threshold}
        
        # Handle special cases
        if operator == "ON":
            sql_expr = f"DATE({field}) = DATE(:param)"
            return sql_expr, {"param": value}
        elif operator == "IN":
            # For IN clauses, we need to create multiple parameters
            param_names = [f"param{i+1}" for i in range(len(value))]
            placeholders = "(" + ",".join([f":{name}" for name in param_names]) + ")"
            sql_expr = f"{field} IN {placeholders}"
            params = {name: val for name, val in zip(param_names, value)}
            return sql_expr, params
        elif operator == "RANGE":
            sql_expr = f"{field} BETWEEN :param1 AND :param2"
            return sql_expr, {"param1": value[0], "param2": value[1]}
        else:
            # Standard operator
            sql_expr = f"{field} {op_def['sql']}"
            return sql_expr, {"param": value}

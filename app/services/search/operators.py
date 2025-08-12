"""
Search operators and their SQL translations
"""

from typing import Dict, Any, Tuple
from datetime import datetime
from app.utils.enums import OperatorType


class SearchOperators:
    """Maps search operators to SQL expressions and validation"""

    # TODO: CHECK OPERATORS FOR VALIDITY AND CONSISTENCY
    # Operator definitions with their SQL translations and parameter handling
    OPERATORS = {
        # String operators
        "=": {"sql": "= :param", "type": OperatorType.STRING, "params": 1, "description": "Exact match"},
        "!=": {"sql": "!= :param", "type": OperatorType.STRING, "params": 1, "description": "Not equal"},
        "LIKE": {
            "sql": "LIKE :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Pattern match with wildcards",
        },
        "CONTAINS": {
            "sql": "ILIKE :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Case-insensitive contains",
            "value_transform": lambda x: f"%{x}%",
        },
        "STARTS WITH": {
            "sql": "ILIKE :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Starts with pattern",
            "value_transform": lambda x: f"{x}%",
        },
        "ENDS WITH": {
            "sql": "ILIKE :param",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Ends with pattern",
            "value_transform": lambda x: f"%{x}",
        },
        "IN": {
            "sql": "{field} IN {params}",
            "type": OperatorType.STRING,
            "params": 1,
            "description": "Match any value in list",
            "value_transform": lambda x: tuple(x) if isinstance(x, list) else x,
        },
        # Numeric operators
        "<": {"sql": "< :param", "type": OperatorType.NUMERIC, "params": 1, "description": "Less than"},
        ">": {"sql": "> :param", "type": OperatorType.NUMERIC, "params": 1, "description": "Greater than"},
        "<=": {"sql": "<= :param", "type": OperatorType.NUMERIC, "params": 1, "description": "Less than or equal"},
        ">=": {"sql": ">= :param", "type": OperatorType.NUMERIC, "params": 1, "description": "Greater than or equal"},
        "RANGE": {
            "sql": "{field} BETWEEN :param1 AND :param2",
            "type": OperatorType.NUMERIC,
            "params": 2,
            "description": "Between two values (inclusive)",
            "value_transform": lambda x: (x[0], x[1]) if isinstance(x, (list, tuple)) and len(x) == 2 else x,
        },
        # Datetime operators
        "BEFORE": {"sql": "< :param", "type": OperatorType.DATETIME, "params": 1, "description": "Before date/time"},
        "AFTER": {"sql": "> :param", "type": OperatorType.DATETIME, "params": 1, "description": "After date/time"},
        "ON": {
            "sql": "DATE({field}) = DATE(:param)",
            "type": OperatorType.DATETIME,
            "params": 1,
            "description": "On specific date",
        },
        # Molecular operators (RDKit)
        "IS SIMILAR": {
            "sql": "public.tanimoto_sml(public.morganbv_fp(public.mol_from_smiles({field}::cstring)), "
            "public.morganbv_fp(public.mol_from_smiles(:param2))) >= :param3",
            "type": OperatorType.MOLECULAR,
            "params": 1,
            "description": "Molecular similarity using RDKit",
            "requires_threshold": True,
        },
        "HAS SUBSTRUCTURE": {
            "sql": "public.mol_from_smiles({field}::cstring) OPERATOR(public.@>) :param\\:\\:public.qmol",
            "type": OperatorType.MOLECULAR,
            "params": 1,
            "description": "Molecular substructure search",
        },
        "IS SUBSTRUCTURE OF": {
            "sql": "public.mol_from_smiles(:param) OPERATOR(public.@>) public.mol_from_smiles({field}::cstring)",
            "type": OperatorType.MOLECULAR,
            "params": 1,
            "description": "Molecular superstructure search",
        },
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
        if operator == "IN" and not isinstance(value, (list, tuple)):
            raise ValueError(f"Operator '{operator}' requires a list or tuple value")

        if operator == "RANGE":
            if not isinstance(value, (list, tuple)) or len(value) != 2:
                raise ValueError(f"Operator '{operator}' requires a list/tuple with exactly 2 values")
            if value[0] >= value[1]:
                raise ValueError("RANGE operator requires first value to be less than second value")

        if operator in ["ON", "BEFORE", "AFTER"]:

            def check_date(date_str: str):
                formats = ["%Y-%m-%d", "%Y-%m-%d %H:%M"]  # Supported formats: date-only and date-time
                for fmt in formats:
                    try:
                        # Try to parse the date string using the YYYY-MM-DD format
                        datetime.strptime(date_str, fmt)
                        return True
                    except ValueError:
                        continue
                return False

            if not check_date(value):
                raise ValueError("DATE value must be in format YYYY-MM-DD or YYYY-MM-DD hh:mm")
        return True

    @classmethod
    def validate_operands(cls, operator: str, field: str):
        """Validate that operand is appropriate for operation"""
        if operator in ["IS SIMILAR", "HAS SUBSTRUCTURE", "IS SUBSTRUCTURE OF"]:
            if not field.endswith(".structure"):
                raise ValueError("Molecular operators can only be applied to compounds.structure")
        else:
            if field.endswith(".structure"):
                raise ValueError("Only molecular operators can be applied to compounds.structure")

    @classmethod
    def get_sql_expression(
        cls, operator: str, field: str, value: Any, threshold: float = None
    ) -> Tuple[str, Dict[str, Any]]:
        """
        Get SQL expression and parameters for an operator

        Returns:
            Tuple of (sql_expression, parameters_dict)
        """
        op_def = cls.get_operator(operator)

        # Transform value if needed
        if "value_transform" in op_def:
            value = op_def["value_transform"](value)

        # Handle special cases
        match operator:
            case "IN":
                # For IN clauses, we need to create multiple parameters
                param_names = [f"param{i + 1}" for i in range(len(value))]
                placeholders = "(" + ",".join([f":{name}" for name in param_names]) + ")"
                sql_expr = op_def["sql"].format(field=field, params=placeholders)
                params = {name: val for name, val in zip(param_names, value)}
                return sql_expr, params
            case "RANGE":
                sql_expr = op_def["sql"].format(field=field)
                return sql_expr, {"param1": value[0], "param2": value[1]}
            case "IS SIMILAR":
                sql_expr = op_def["sql"].format(field=field)
                return sql_expr, {"param1": value, "param2": value, "param3": threshold}
            case "IS SUBSTRUCTURE OF" | "HAS SUBSTRUCTURE" | "ON":
                sql_expr = op_def["sql"].format(field=field)
                return sql_expr, {"param": value}
            case _:
                # Standard operator
                sql_expr = f"{field} {op_def['sql']}"
                return sql_expr, {"param": value}

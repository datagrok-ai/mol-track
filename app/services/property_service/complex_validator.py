from datetime import datetime
from typing import Dict, Mapping, List, Any
from cel import evaluate, Context
import re


class ComplexValidationError(ValueError):
    """Raised when a record does not satisfy a validation expression."""


class ComplexValidator:
    """
    Validates entire records (multi-field objects) using CEL expressions.
    Rules can use placeholders like ${field}, ${field}.length, is null, etc.
    """

    @classmethod
    def validate_record(cls, record: Mapping[str, Any], rules: List[str] = None, validate_rule=False) -> None:
        """
        Validate a record against a list of CEL rules.
        Raises RecordValidationError if any rule fails.
        """
        if not rules:
            return
        safe_ctx = cls._sanitize_context(record)
        ctx = cls._build_context(safe_ctx)

        for raw_expr in rules:
            evaluate_rule = True
            variables = cls._extract_variables(raw_expr)
            for var in variables:
                if var not in record.keys():
                    evaluate_rule = False
                    break
            if not evaluate_rule:
                continue
            expr = cls._preprocess(raw_expr)
            try:
                result = evaluate(expr, ctx)
            except Exception as e:
                raise ComplexValidationError(
                    f"Error while evaluating rule '{raw_expr}' (translated to '{expr}'): {e}"
                ) from e

            if not validate_rule and not bool(result):
                raise ComplexValidationError(f"Record does not satisfy rule: {raw_expr}")

    @classmethod
    def validate_rule(cls, expr: str, properties: Dict[str, str]) -> bool:
        """
        Validate a single CEL rule against a dictionary of properties.
        Creates a mock context where each property is set to a non-null value.
        Returns True if the rule compiles, False otherwise.
        """

        variables = cls._extract_variables(expr)
        for var in variables:
            if var not in properties.keys():
                raise ComplexValidationError(f"Unknown property '{var}' in rule: {expr}")

        mock_record = cls._create_mock_record(variables, properties)

        cls.validate_record(mock_record, [expr], validate_rule=True)

    # ------------------------
    # Internal helpers
    # ------------------------
    @staticmethod
    def _sanitize_context(record: Mapping[str, Any]) -> Mapping[str, Any]:
        """Convert record values into JSON-serializable primitives."""
        safe_ctx = {}
        for k, v in record.items():
            k = k.replace(" ", "_")
            if isinstance(v, (str, int, float, bool)) or v is None:
                safe_ctx[k] = v
            else:
                safe_ctx[k] = str(v)
        return safe_ctx

    @staticmethod
    def _build_context(record: Mapping[str, Any]) -> Context:
        """
        Build a CEL Context with custom functions.
        """
        ctx = Context(record)

        # size(value) → length of string or list
        ctx.add_function("size", lambda x: len(x) if x is not None else 0)

        # matches(value, pattern) → regex match
        ctx.add_function("matches", lambda val, pat: re.match(pat, val) is not None)

        # today() → current date as ISO string
        ctx.add_function("today", lambda: datetime.today().date().isoformat())

        # date(str) → convert string to date object
        ctx.add_function("date", lambda s: datetime.fromisoformat(s).date())

        return ctx

    @staticmethod
    def _preprocess(expr: str) -> str:
        """Translate DSL-style expressions into CEL syntax."""
        # ${field} → field
        expr = re.sub(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}", r"\1", expr)

        # ${field}.length → size(field)
        expr = re.sub(r"([a-zA-Z_][a-zA-Z0-9_]*)\.length", r"size(\1)", expr)

        # is null → == null
        expr = re.sub(r"([a-zA-Z_][a-zA-Z0-9_]*)\s+is\s+null", r"\1 == null", expr)

        # is not null → != null
        expr = re.sub(r"([a-zA-Z_][a-zA-Z0-9_]*)\s+is\s+not\s+null", r"\1 != null", expr)

        # replaces white spaces inside ${...} with underscores
        expr = re.sub(r"\$\{([^}]+)\}", lambda m: re.sub(r"[^a-zA-Z0-9_]", "_", m.group(1)), expr)

        return expr.strip()

    @staticmethod
    def _extract_variables(expr: str) -> list[str]:
        """
        Extract all variables from an expression.
        """

        pattern = re.compile(r"\$\{([^}]+)\}")
        vars = pattern.findall(expr)
        vars_cel = {re.sub(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}", r"\1", v) for v in vars}
        return list(vars_cel)

    @staticmethod
    def _create_mock_record(variables, properties) -> dict[str, Any]:
        """
        Create a mock record with sample values for testing.
        """

        mock_values = {
            "string": "mock",
            "int": 2,
            "double": 2.0,
            "bool": True,
            "datetime": datetime.now().isoformat(),
            "uuid": "123e4567-e89b-12d3-a456-426614174000",
        }

        record = {}
        for var in variables:
            value = mock_values.get(properties.get(var, "string"), "mock")
            record[var] = value

        return record

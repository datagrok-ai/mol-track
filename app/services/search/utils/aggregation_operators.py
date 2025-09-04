from app.utils.enums import AggregationNumericOp, AggregationStringOp


class AggregationOperators:
    """Maps aggregation operators to SQL expressions and validation"""

    OPERATORS = {
        AggregationNumericOp.FIRST.value: "MIN({column})",
        AggregationNumericOp.TOTAL_COUNT.value: "COUNT(*)",
        AggregationNumericOp.VALUE_COUNT.value: "COUNT({column})",
        AggregationNumericOp.UNIQUE_COUNT.value: "COUNT(DISTINCT {column})",
        AggregationNumericOp.MISSING_VALUE_COUNT.value: "SUM(({column} IS NULL)::int)",
        AggregationNumericOp.MIN.value: "MIN({column})",
        AggregationNumericOp.MAX.value: "MAX({column})",
        AggregationNumericOp.SUM.value: "SUM({column})",
        AggregationNumericOp.MED.value: "PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY {column})",
        AggregationNumericOp.AVG.value: "AVG({column})",
        AggregationNumericOp.STDEV.value: "STDDEV_SAMP({column})",
        AggregationNumericOp.VARIANCE.value: "VAR_SAMP({column})",
        AggregationNumericOp.SKEW.value: None,
        AggregationNumericOp.KURT.value: None,
        AggregationNumericOp.Q1.value: "PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY {column})",
        AggregationNumericOp.Q2.value: "PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY {column})",
        AggregationNumericOp.Q3.value: "PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY {column})",
        AggregationStringOp.CONCAT_ALL.value: "STRING_AGG({column}, ', ')",
        AggregationStringOp.CONCAT_UNIQUE.value: "STRING_AGG(DISTINCT {column}, ', ')",
        AggregationStringOp.LONGEST.value: "(SPLIT_PART(STRING_AGG(CASE WHEN {condition} THEN {column} END, ',' ORDER BY length({column}) DESC), ',', 1))",
        AggregationStringOp.SHORTEST.value: "(SPLIT_PART(STRING_AGG(CASE WHEN {condition} THEN {column} END, ',' ORDER BY length({column}) ASC), ',', 1))",
        AggregationStringOp.MOST_FREQUENT.value: "MODE() WITHIN GROUP (ORDER BY {column}) ",
        AggregationStringOp.CONCAT_COUNTS.value: None,
    }

    @classmethod
    def get_sql_expression(cls, operator: str, column: str, condition: str) -> str:
        if not operator:
            return f"MAX({column})"
        match operator:
            case (
                AggregationNumericOp.SKEW.value
                | AggregationNumericOp.KURT.value
                | AggregationStringOp.CONCAT_COUNTS.value
            ):
                raise ValueError(f"Operator {operator} is not supported for SQL aggregation yet.")
            case (
                AggregationStringOp.LONGEST.value | AggregationStringOp.SHORTEST.value | AggregationNumericOp.SKEW.value
            ):
                return cls.OPERATORS[operator].format(column=column, condition=condition)
            case _:
                return cls.OPERATORS[operator].format(column=column)

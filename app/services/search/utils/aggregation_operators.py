from app.utils.enums import AggregationNumericOp, AggregationStringOp


class AggregationOperators:
    """Maps aggregation operators to SQL expressions and validation"""

    OPERATORS = {
        AggregationNumericOp.TOTAL_COUNT.value: "COUNT(*)",
        AggregationNumericOp.VALUE_COUNT.value: "COUNT({column})",
        AggregationNumericOp.UNIQUE_COUNT.value: "COUNT(DISTINCT {column})",
        AggregationNumericOp.MISSING_VALUE_COUNT.value: "SUM(({column} IS NULL)::int)",
        AggregationNumericOp.MIN.value: "MIN(CAST({column} AS NUMERIC))",
        AggregationNumericOp.MAX.value: "MAX(CAST({column} AS NUMERIC))",
        AggregationNumericOp.SUM.value: "SUM(CAST({column} AS NUMERIC))",
        AggregationNumericOp.MED.value: "PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY CAST({column} AS NUMERIC))",
        AggregationNumericOp.AVG.value: "AVG(CAST({column} AS NUMERIC))",
        AggregationNumericOp.STDEV.value: "STDDEV_SAMP(CAST({column} AS NUMERIC))",
        AggregationNumericOp.VARIANCE.value: "VAR_SAMP(CAST({column} as NUMERIC))",
        AggregationNumericOp.STDDEV_POP.value: "STDDEV_POP(CAST({column} AS NUMERIC))",
        AggregationNumericOp.VAR_POP.value: "VAR_POP(CAST({column} as NUMERIC))",
        AggregationNumericOp.Q1.value: "PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY CAST({column} AS NUMERIC))",
        AggregationNumericOp.Q2.value: "PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY CAST({column} AS NUMERIC))",
        AggregationNumericOp.Q3.value: "PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY CAST({column} AS NUMERIC))",
        AggregationNumericOp.ARRAY_AGG.value: "ARRAY_AGG({column})",
        AggregationStringOp.CONCAT_ALL.value: "STRING_AGG({column}, ', ')",
        AggregationStringOp.CONCAT_UNIQUE.value: "STRING_AGG(DISTINCT {column}, ', ')",
        AggregationStringOp.MOST_FREQUENT.value: "MODE() WITHIN GROUP (ORDER BY {column}) ",
    }

    @classmethod
    def get_sql_expression(cls, operator: str | None, column: str, condition: str) -> str:
        if not operator:
            return f"MAX({column})"

        try:
            return cls.OPERATORS[operator].format(column=column)
        except KeyError:
            raise ValueError(f"Unsupported operator: {operator}")

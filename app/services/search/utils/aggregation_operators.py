from app.utils.enums import AggregationNumericOp, AggregationStringOp


class AggregationOperators:
    """Maps aggregation operators to SQL expressions and validation"""

    OPERATORS = {
        AggregationNumericOp.FIRST.value: "FIRST_VALUE({column}) OVER (ORDER BY some_order)",
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
        AggregationNumericOp.SKEW.value: "(AVG(POWER({column} - AVG({column}), 3)) / POWER(STDDEV_SAMP({column}), 3))",
        AggregationNumericOp.KURT.value: "(AVG(POWER({column} - AVG({column}), 4)) / POWER(STDDEV_SAMP({column}), 4))",
        AggregationNumericOp.Q1.value: "PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY {column})",
        AggregationNumericOp.Q2.value: "PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY {column})",
        AggregationNumericOp.Q3.value: "PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY {column})",
        AggregationStringOp.CONCAT_ALL.value: "STRING_AGG({column}, ', ')",
        AggregationStringOp.CONCAT_UNIQUE.value: "STRING_AGG(DISTINCT {column}, ', ')",
        AggregationStringOp.LONGEST.value: "",
        AggregationStringOp.SHORTEST.value: "",
        AggregationStringOp.MOST_FREQUENT.value: "",
        AggregationStringOp.CONCAT_COUNTS.value: "",
    }

    @classmethod
    def get_sql_expression(cls, operator: str, column: str) -> str:
        if not operator:
            operator = AggregationNumericOp.MAX.value
        return cls.OPERATORS[operator].format(column=column)

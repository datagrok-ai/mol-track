import enum


class CaseInsensitiveEnum(str, enum.Enum):
    @classmethod
    def _missing_(cls, value):
        if isinstance(value, str):
            for member in cls:
                if member.value.lower() == value.lower():
                    return member
        return None


class ValueType(str, enum.Enum):
    int = "int"
    double = "double"
    bool = "bool"
    datetime = "datetime"
    string = "string"
    uuid = "uuid"


class PropertyClass(CaseInsensitiveEnum):
    DECLARED = "DECLARED"
    CALCULATED = "CALCULATED"
    MEASURED = "MEASURED"
    PREDICTED = "PREDICTED"


class EntityType(CaseInsensitiveEnum):
    BATCH = "BATCH"
    COMPOUND = "COMPOUND"
    ASSAY = "ASSAY"
    ASSAY_RUN = "ASSAY_RUN"
    ASSAY_RESULT = "ASSAY_RESULT"
    SYSTEM = "SYSTEM"


class SearchEntityType(str, enum.Enum):
    BATCHES = "batches"
    COMPOUNDS = "compounds"
    ASSAYS = "assays"
    ASSAY_RUNS = "assay_runs"
    ASSAY_RESULTS = "assay_results"


class EntityTypeReduced(CaseInsensitiveEnum):
    BATCH = "BATCH"
    COMPOUND = "COMPOUND"


class AdditionsRole(CaseInsensitiveEnum):
    SALT = "SALT"
    SOLVATE = "SOLVATE"


class SynonymLevel(CaseInsensitiveEnum):
    BATCH = "BATCH"
    COMPOUND = "COMPOUND"


class ErrorHandlingOptions(str, enum.Enum):
    reject_all = "reject_all"
    reject_row = "reject_row"


class OutputFormat(str, enum.Enum):
    json = "json"
    csv = "csv"
    sdf = "sdf"


class SearchOutputFormat(str, enum.Enum):
    json = "json"
    csv = "csv"
    parquet = "parquet"


class CompoundMatchingRule(CaseInsensitiveEnum):
    ALL_LAYERS = "ALL_LAYERS"
    STEREO_INSENSITIVE_LAYERS = "STEREO_INSENSITIVE_LAYERS"
    TAUTOMER_INSENSITIVE_LAYERS = "TAUTOMER_INSENSITIVE_LAYERS"


class LogicOp(str, enum.Enum):
    """Logical operators for combining conditions"""

    AND = "AND"
    OR = "OR"

    @classmethod
    def _missing_(cls, value):
        return lowercase_enum_values(cls, value)


class CompareOp(str, enum.Enum):
    """Comparison operators for atomic conditions"""

    # String operators
    EQUALS = "="
    NOT_EQUALS = "!="
    IN = "IN"
    STARTS_WITH = "STARTS WITH"
    ENDS_WITH = "ENDS WITH"
    LIKE = "LIKE"
    CONTAINS = "CONTAINS"

    # Numeric operators
    LESS_THAN = "<"
    GREATER_THAN = ">"
    LESS_THAN_OR_EQUAL = "<="
    GREATER_THAN_OR_EQUAL = ">="
    RANGE = "RANGE"

    # Datetime operators
    BEFORE = "BEFORE"
    AFTER = "AFTER"
    ON = "ON"

    # Molecular operators (RDKit)
    IS_SIMILAR = "IS SIMILAR"
    IS_SUBSTRUCTURE_OF = "IS SUBSTRUCTURE OF"
    HAS_SUBSTRUCTURE = "HAS SUBSTRUCTURE"

    @classmethod
    def _missing_(cls, value):
        return lowercase_enum_values(cls, value)


class OperatorType(enum.Enum):
    """Types of operators for different data types"""

    STRING = "string"
    NUMERIC = "numeric"
    DATETIME = "datetime"
    BOOLEAN = "boolean"
    MOLECULAR = "molecular"


class SettingName(str, enum.Enum):
    COMPOUND_MATCHING_RULE = "COMPOUND_MATCHING_RULE"
    COMPOUND_SEQUENCE_START = "COMPOUND_SEQUENCE_START"
    BATCH_SEQUENCE_START = "BATCH_SEQUENCE_START"
    CORPORATE_COMPOUND_ID_PATTERN = "CORPORATE_COMPOUND_ID_PATTERN"
    CORPORATE_BATCH_ID_PATTERN = "CORPORATE_BATCH_ID_PATTERN"


class AggregationNumericOp(str, enum.Enum):
    FIRST = "FIRST"
    TOTAL_COUNT = "COUNT"
    VALUE_COUNT = "VALUES"
    UNIQUE_COUNT = "UNIQUE"
    MISSING_VALUE_COUNT = "NULLS"
    MIN = "MIN"
    MAX = "MAX"
    SUM = "SUM"
    MED = "MED"
    AVG = "AVG"
    STDEV = "STDEV"
    VARIANCE = "VARIANCE"
    SKEW = "SKEW"
    KURT = "KURT"
    Q1 = "Q1"
    Q2 = "Q2"
    Q3 = "Q3"

    @classmethod
    def _missing_(cls, value):
        return lowercase_enum_values(cls, value)


class AggregationStringOp(str, enum.Enum):
    CONCAT_ALL = "CONCAT ALL"
    CONCAT_UNIQUE = "CONCAT UNIQUE"
    LONGEST = "LONGEST"
    SHORTEST = "SHORTEST"
    MOST_FREQUENT = "MOST FREQUENT"
    CONCAT_COUNTS = "CONCAT COUNTS"

    @classmethod
    def _missing_(cls, value):
        return lowercase_enum_values(cls, value)


def lowercase_enum_values(enum_class, value):
    if isinstance(value, str):
        v_norm = value.strip().upper()
        for member in enum_class:
            if member.value == v_norm:
                return member
    return None

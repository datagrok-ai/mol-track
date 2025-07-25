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


class ScopeClass(CaseInsensitiveEnum):
    BATCH = "BATCH"
    COMPOUND = "COMPOUND"
    ASSAY = "ASSAY"
    ASSAY_RUN = "ASSAY_RUN"
    ASSAY_RESULT = "ASSAY_RESULT"
    SYSTEM = "SYSTEM"


class ScopeClassReduced(CaseInsensitiveEnum):
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


class CompoundMatchingRule(CaseInsensitiveEnum):
    ALL_LAYERS = "ALL_LAYERS"
    STEREO_INSENSITIVE_LAYERS = "STEREO_INSENSITIVE_LAYERS"
    TAUTOMER_INSENSITIVE_LAYERS = "TAUTOMER_INSENSITIVE_LAYERS"


class LogicOp(str, enum.Enum):
    """Logical operators for combining conditions"""

    AND = "AND"
    OR = "OR"


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


class OperatorType(enum.Enum):
    """Types of operators for different data types"""

    STRING = "string"
    NUMERIC = "numeric"
    DATETIME = "datetime"
    BOOLEAN = "boolean"
    MOLECULAR = "molecular"

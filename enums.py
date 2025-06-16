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
    CALCULATED = "CALCULATED"
    MEASURED = "MEASURED"
    PREDICTED = "PREDICTED"
    ASSERTED = "ASSERTED"


class ScopeClass(CaseInsensitiveEnum):
    BATCH = "BATCH"
    COMPOUND = "COMPOUND"
    ASSAY = "ASSAY"
    ASSAY_TYPES = "ASSAY_TYPES"
    ASSAY_RESULTS = "ASSAY_RESULTS"
    SYSTEM = "SYSTEM"


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

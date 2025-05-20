import enum


class ValueType(str, enum.Enum):
    int = "int"
    double = "double"
    bool = "bool"
    datetime = "datetime"
    string = "string"


class PropertyClass(str, enum.Enum):
    CALCULATED = "CALCULATED"
    MEASURED = "MEASURED"
    PREDICTED = "PREDICTED"

class ScopeClass(str, enum.Enum):
    BATCH  = "BATCH"
    COMPOUND = "COMPOUND"
    ASSAY = "ASSAY"
    SYSTEM = "SYSTEM" 
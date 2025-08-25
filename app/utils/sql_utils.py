from typing import List, Dict, Any
from psycopg2.extensions import adapt
from sqlmodel import SQLModel
from app.utils import enums

column_types = {
    "value_datetime": "timestamptz",
    "value_num": "double precision",
    "value_uuid": "uuid",
    "value_string": "text",
    "value_bool": "boolean",
    "value_qualifier": "smallint",
}


def values_sql(data: List[Dict[str, Any]], columns: List[str]) -> str:
    def escape_val(val, col_name):
        if val is None:
            sql_type = column_types.get(col_name)
            if sql_type:
                return f"NULL::{sql_type}"
            return "NULL"
        return adapt(val).getquoted().decode()

    rows = []
    for row in data:
        values = [escape_val(row.get(col), col) for col in columns]
        rows.append(f"({', '.join(values)})")
    return ",\n".join(rows)


def prepare_sql_parts(records: List[Dict[str, Any]]):
    cols = list(records[0].keys())
    key, *cols_without_key = cols
    values_to_sql = values_sql(records, cols)
    return cols_without_key, values_to_sql


def generate_sql(*sql_parts: str, terminate_with_select: bool = True) -> str:
    filtered_parts = [part.strip() for part in sql_parts if part and part.strip()]
    if not filtered_parts:
        return ""
    combined_sql = ",\n".join(filtered_parts)
    combined_sql += "\nSELECT 1;" if terminate_with_select else ";"
    return combined_sql


def chunked(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i : i + size]


def normalize_type(type_name: str) -> str:
    t = type_name.lower()

    for vt in enums.ValueType:
        if vt.value in t:
            return vt.value
    return enums.ValueType.string.value


def get_table_fields(table_name: str) -> list[dict[str, Any]] | None:
    mapping = {
        "batches": enums.EntityType.BATCH,
        "compounds": enums.EntityType.COMPOUND,
        "assays": enums.EntityType.ASSAY,
        "assay_runs": enums.EntityType.ASSAY_RUN,
        "assay_results": enums.EntityType.ASSAY_RESULT,
    }

    table = SQLModel.metadata.tables.get(f"moltrack.{table_name}")
    return [
        {
            "entity_type": mapping.get(table_name),
            "name": col.name,
            "value_type": normalize_type(col.type.__class__.__name__),
        }
        for col in table.columns
    ]


def get_direct_fields():
    return [field for table_name in enums.SearchEntityType for field in get_table_fields(table_name.value)]

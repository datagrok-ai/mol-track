from datetime import datetime
import uuid
from typing import List, Dict, Any


class SQLService:
    def values_sql(self, data: List[Dict[str, Any]], columns: List[str]) -> str:
        def quote(val: Any) -> str:
            if val is None:
                return "NULL"
            if isinstance(val, str):
                escaped = val.replace("'", "''")
                return f"'{escaped}'"
            if isinstance(val, datetime):
                return f"'{val.isoformat()}'::timestamp with time zone"
            if isinstance(val, uuid.UUID):
                return f"'{val}'::uuid"
            return str(val)

        rows = []
        for row in data:
            values = [quote(row.get(col)) for col in columns]
            rows.append(f"({', '.join(values)})")
        return ",\n".join(rows)

    def prepare_sql_parts(self, records: List[Dict[str, Any]]):
        cols = list(records[0].keys())
        key, *cols_without_key = cols
        values_sql = self.values_sql(records, cols)
        return cols_without_key, values_sql

    def generate_sql(self, *sql_parts: str, terminate_with_select: bool = True) -> str:
        filtered_parts = [part.strip() for part in sql_parts if part and part.strip()]
        if not filtered_parts:
            return ""
        combined_sql = ",\n".join(filtered_parts)
        combined_sql += "\nSELECT 1;" if terminate_with_select else ";"
        return combined_sql

    def chunked(self, lst, size):
        for i in range(0, len(lst), size):
            yield lst[i : i + size]

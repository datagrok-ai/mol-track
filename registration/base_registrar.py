import csv
import io
import json
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
import uuid
from datetime import datetime

from fastapi.responses import StreamingResponse
from sqlalchemy import select, text
from fastapi import HTTPException

import enums


class BaseRegistrar(ABC):
    def __init__(self, db, mapping: Optional[str], error_handling: str):
        """
        Base class for processing and registering data to a database.
        :param db: SQLAlchemy database session.
        :param mapping: Optional JSON string defining field mappings.
        :param error_handling: Strategy for handling errors during processing.
        """
        self.db = db
        self.error_handling = error_handling
        self.user_mapping = self._load_mapping(mapping)
        self.output_records: List[Dict[str, Any]] = []
        self.sql_statements = []

    # === Input processing methods ===

    def _load_mapping(self, mapping: Optional[str]) -> Dict[str, str]:
        if not mapping:
            return {}
        try:
            return json.loads(mapping)
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid JSON for mapping")

    def process_csv(self, csv_content: str) -> List[Dict[str, Any]]:
        rows = list(csv.DictReader(io.StringIO(csv_content), skipinitialspace=True))
        if not rows:
            raise HTTPException(status_code=400, detail="CSV file is empty or invalid")
        self.normalized_mapping = self.user_mapping or {k: k for k in rows[0]}
        return rows

    def _group_data(self, row: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        grouped = {}
        for src_key, mapped_key in self.normalized_mapping.items():
            value = row.get(src_key)
            table, field = mapped_key.split(".", 1) if "." in mapped_key else ("compounds", mapped_key)
            grouped.setdefault(table, {})[field] = value
        return grouped

    # === Reference loading methods ===

    def _load_reference_map(self, model, key: str = "id"):
        result = self.db.execute(select(model)).scalars().all()
        return {getattr(row, key): row for row in result}

    def model_to_dict(self, obj):
        return {c.key: getattr(obj, c.key) for c in obj.__table__.columns}

    # === SQL construction methods ===

    def _values_sql(self, data: List[Dict[str, Any]], columns: List[str]) -> str:
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

    def register_all(self, rows: List[Dict[str, Any]]):
        self.build_sql(rows)
        if self.sql_statements:
            for sql in self.sql_statements:
                self.db.execute(text(sql))
                self.db.commit()

    @abstractmethod
    def build_sql(self, rows: List[Dict[str, Any]]):
        pass

    @abstractmethod
    def generate_sql(self) -> Optional[str]:
        pass

    # === Output formatting methods ===

    def _add_output_row(self, compound_data, grouped, status, error_msg=None):
        output = {
            **compound_data,
            **{f"property_{k}": v for k, v in grouped.get("compounds_details", {}).items()},
            "registration_status": status,
            "registration_error_message": error_msg,
        }
        if hasattr(self, "get_additional_output_info"):
            output.update(self.get_additional_output_info(grouped))
        self.output_records.append(output)

    def result(self, output_format: str = enums.OutputFormat.json) -> Dict[str, str]:
        def get_csv() -> str:
            output = io.StringIO()
            fieldnames = list({key for rec in self.output_records for key in rec})
            writer = csv.DictWriter(output, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for row in self.output_records:
                writer.writerow(row)
            csv_data = output.getvalue()
            output.close()
            return csv_data

        if output_format == enums.OutputFormat.csv:
            csv_data = get_csv()
            return StreamingResponse(
                io.StringIO(csv_data),
                media_type="text/csv",
                headers={"Content-Disposition": "attachment; filename=compounds_result.csv"},
            )
        return {"status": "Success", "data": self.output_records}

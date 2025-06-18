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
import utils
import main
import models


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
        self.property_records_map = self._load_reference_map(models.Property, "name")
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

    # === Some methods that can be shared across all registrars
    def _build_details_records(
        self, properties: Dict[str, Any], entity_ids: Dict[str, Any], include_user_fields: bool = True
    ) -> List[Dict[str, Any]]:
        records = []
        for prop_name, value in properties.items():
            prop = self.property_records_map.get(prop_name)

            if prop is None:
                raise HTTPException(status_code=400, detail=f"Unknown property: {prop_name}")

            value_type = getattr(prop, "value_type", None)
            if value_type not in utils.value_type_to_field or value_type not in utils.value_type_cast_map:
                raise HTTPException(
                    status_code=400, detail=f"Unsupported or unknown value type for property: {prop_name}"
                )

            field_name = utils.value_type_to_field[value_type]
            cast_fn = utils.value_type_cast_map[value_type]
            detail = {
                **entity_ids,
                "property_id": getattr(prop, "id"),
                "value_datetime": datetime.now(),
                "value_num": 0,
                "value_string": None,
            }

            if include_user_fields:
                detail["created_by"] = main.admin_user_id
                detail["updated_by"] = main.admin_user_id

            try:
                detail[field_name] = cast_fn(value)
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Error casting value for property {prop_name}: {e}")

            records.append(detail)
        return records

    def _prepare_sql_parts(self, records: List[Dict[str, Any]]):
        cols = list(records[0].keys())
        key, *cols_without_key = cols
        values_sql = self._values_sql(records, cols)
        return cols_without_key, values_sql

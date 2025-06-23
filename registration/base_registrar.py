import csv
import io
import json
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional

from fastapi.responses import StreamingResponse
from sqlalchemy import select, text
from fastapi import HTTPException

import enums
from services import property_service
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
        self.property_records_map = self._load_reference_map(models.Property, "name")
        self.property_service = property_service.PropertyService(self.property_records_map)

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

    # === SQL construction and registration methods ===

    @abstractmethod
    def build_sql(self, rows: List[Dict[str, Any]]):
        pass

    @abstractmethod
    def generate_sql(self) -> Optional[str]:
        pass

    def register_all(self, rows: List[Dict[str, Any]]):
        self.build_sql(rows)
        if self.sql_statements:
            for sql in self.sql_statements:
                try:
                    self.db.execute(text(sql))
                    self.db.commit()
                except Exception:
                    self.db.rollback()

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

    # === Error handling methods ===

    def handle_row_error(self, row, exception, global_idx, all_rows):
        self._add_output_row(row, {}, "failed", str(exception))
        if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
            for remaining_row in all_rows[global_idx + 1 :]:
                self._add_output_row(remaining_row, {}, "not_processed")
            raise HTTPException(status_code=400, detail=self.result())

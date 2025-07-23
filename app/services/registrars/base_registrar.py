import csv
import io
import json
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional

from fastapi.responses import StreamingResponse
from sqlalchemy import select, text
from fastapi import HTTPException

from app.utils import enums
from app.utils.logging_utils import logger
from app.services import property_service
from app import models


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
        self.addition_records_map = self._load_reference_map(models.Addition, "name")
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

        if self.user_mapping:
            self.normalized_mapping = self.user_mapping
        else:
            self.normalized_mapping = {}
            for col in rows[0].keys():
                assigned = self._assign_column(col)
                self.normalized_mapping[col] = assigned
        return rows

    def _assign_column(self, col: str) -> str:
        if col in self.property_records_map:
            scope = getattr(self.property_records_map[col], "scope", None)
            prefix = {
                enums.ScopeClass.COMPOUND: "compound_details",
                enums.ScopeClass.BATCH: "batch_details",
                enums.ScopeClass.ASSAY_RUN: "assay_run_details",
                enums.ScopeClass.ASSAY_RESULT: "assay_results",
            }.get(scope)
            return f"{prefix}.{col}" if prefix else col

        if col in self.addition_records_map:
            return f"batch_additions.{col}"
        return col

    def _group_data(self, row: Dict[str, Any], entity_name: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
        grouped = {}
        for src_key, mapped_key in self.normalized_mapping.items():
            value = row.get(src_key)
            table, field = (
                mapped_key.split(".", 1)
                if "." in mapped_key
                else (entity_name if entity_name else "compound", mapped_key)
            )
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
                except Exception as e:
                    logger.error(f"An exception occurred: {e}")
                    self.db.rollback()

    # === Output formatting methods ===

    def _add_output_row(self, compound_data, grouped, status, error_msg=None):
        output = {
            **compound_data,
            **{f"property_{k}": v for k, v in grouped.get("compound_details", {}).items()},
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

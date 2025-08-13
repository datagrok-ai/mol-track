import csv
import io
import json
from abc import ABC, abstractmethod
from typing import Iterator, List, Dict, Any, Optional

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
        self._property_records_map = None
        self._addition_records_map = None

        self.property_service = property_service.PropertyService(self.property_records_map)
        self.user_mapping = self._load_mapping(mapping)
        self.output_rows = []

    @property
    def property_records_map(self):
        if self._property_records_map is None:
            self._property_records_map = self._load_reference_map(models.Property, "name")
        return self._property_records_map

    @property
    def addition_records_map(self):
        if self._addition_records_map is None:
            self._addition_records_map = self._load_reference_map(models.Addition, "name")
        return self._addition_records_map

    # === Input processing methods ===

    def _load_mapping(self, mapping: Optional[str]) -> Dict[str, str]:
        if not mapping:
            return {}
        try:
            return json.loads(mapping)
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid JSON for mapping")

    def process_csv(self, file_stream: io.TextIOBase, chunk_size=5000) -> Iterator[List[Dict[str, Any]]]:
        reader = csv.DictReader(file_stream, skipinitialspace=True)

        try:
            first_row = next(reader)
        except StopIteration:
            raise HTTPException(status_code=400, detail="CSV file is empty or invalid")

        if self.user_mapping:
            self.normalized_mapping = self.user_mapping
        else:
            self.normalized_mapping = {}
            for col in first_row.keys():
                assigned = self._assign_column(col)
                self.normalized_mapping[col] = assigned

        chunk = [first_row]

        for row in reader:
            chunk.append(row)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []

        if chunk:
            yield chunk

    def _assign_column(self, col: str) -> str:
        if col in self.property_records_map:
            entity_type = getattr(self.property_records_map[col], "entity_type", None)
            prefix = {
                enums.EntityType.COMPOUND: "compound_details",
                enums.EntityType.BATCH: "batch_details",
                enums.EntityType.ASSAY_RUN: "assay_run_details",
                enums.EntityType.ASSAY_RESULT: "assay_result_details",
            }.get(entity_type)
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
    def build_sql(self, rows: List[Dict[str, Any]]) -> str:
        pass

    @abstractmethod
    def generate_sql(self) -> Optional[str]:
        pass

    def register_all(self, rows: List[Dict[str, Any]]):
        batch_sql = self.build_sql(rows)

        if batch_sql:
            try:
                self.db.execute(text(batch_sql))
                self.db.commit()
            except Exception as e:
                logger.error(f"An exception occurred: {e}")
                self.db.rollback()

    # === Output formatting methods ===

    def _add_output_row(self, row, status, error_msg=None):
        row["registration_status"] = status
        row["registration_error_message"] = error_msg or ""
        self.output_rows.append(row)

    def cleanup_chunk(self):
        self.output_rows.clear()

    def cleanup(self):
        self.cleanup_chunk()
        self.user_mapping.clear()
        self.normalized_mapping.clear()
        self._property_records_map = None
        self._addition_records_map = None
        self.property_service = None

    # === Error handling methods ===

    def handle_row_error(self, row, exception, global_idx, all_rows):
        self._add_output_row(row, "failed", str(exception))
        if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
            for remaining_row in all_rows[global_idx + 1 :]:
                self._add_output_row(remaining_row, "not_processed")
            raise HTTPException(status_code=400, detail="Registration failed. See the output file for details.")

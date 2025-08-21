from datetime import datetime
from typing import Any, Dict, List, Optional
from fastapi import HTTPException
from pytest import Session
from app.services.registrars.compound_registrar import CompoundRegistrar
from app.utils.admin_utils import admin
from app import models
from app.utils import enums, sql_utils
from sqlalchemy.sql import text


class BatchRegistrar(CompoundRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self._additions_map = None

        self.batches_to_insert = []
        self.batch_details = []
        self.batch_additions = []

        self.entity_type = enums.EntityType.BATCH

    @property
    def additions_map(self):
        if self._additions_map is None:
            self._additions_map = self._load_reference_map(models.Addition, "name")
        return self._additions_map

    def _next_batch_regno(self) -> int:
        return self.db.execute(text("SELECT nextval('moltrack.batch_regno_seq');")).scalar()

    def _build_batch_record(self, inchikey: str) -> Dict[str, Any]:
        return {
            "inchikey": inchikey,
            "notes": None,
            "created_by": admin.admin_user_id,
            "updated_by": admin.admin_user_id,
            "created_at": datetime.now(),
            "batch_regno": self._next_batch_regno(),
        }

    def _build_batch_addition_record(self, batch_additions: Dict[str, Any], batch_regno: int) -> List[Dict[str, Any]]:
        records = []
        for name, value in batch_additions.items():
            if value in ("", "none", None):
                continue

            addition = self.additions_map.get(name)
            if addition is None:
                raise HTTPException(status_code=400, detail=f"Unknown addition: {name}")
            records.append(
                {
                    "batch_regno": batch_regno,
                    "addition_id": getattr(addition, "id"),
                    "addition_equivalent": float(value),
                    "created_by": admin.admin_user_id,
                    "updated_by": admin.admin_user_id,
                }
            )
        return records

    def get_additional_records(self, row, grouped, molregno):
        batch_record = self._build_batch_record(molregno)
        batch_regno = batch_record["batch_regno"]

        self.inject_corporate_property(row, grouped, batch_regno, enums.EntityType.BATCH)
        inserted, updated = self.property_service.build_details_records(
            models.BatchDetail,
            grouped.get("batch_details", {}),
            {"batch_regno": batch_regno},
            enums.EntityType.BATCH,
        )
        additions = self._build_batch_addition_record(grouped.get("batch_additions", {}), batch_regno)

        self.batches_to_insert.append(batch_record)
        self.batch_details.extend(inserted)
        self.batch_additions.extend(additions)

    def get_additional_cte(self):
        if not self.batches_to_insert:
            return ""
        return self._build_batch_ctes(self.batches_to_insert, self.batch_details, self.batch_additions)

    def _build_batch_ctes(self, batches, details, additions) -> str:
        batch_cte = self._build_inserted_batches_cte(batches)
        if details:
            batch_cte += self._build_batch_details_cte(details)
        if additions:
            batch_cte += self._build_batch_additions_cte(additions)

        return batch_cte

    def _build_inserted_batches_cte(self, batches) -> str:
        cols_without_key, values_sql = sql_utils.prepare_sql_parts(batches)
        return f"""
            inserted_batches AS (
                INSERT INTO moltrack.batches (compound_id, {", ".join(cols_without_key)})
                SELECT ic.id, {", ".join([f"b.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS b (molregno, {", ".join(cols_without_key)})
                JOIN available_compounds ic ON b.molregno = ic.molregno
                ON CONFLICT (batch_regno) DO NOTHING
                RETURNING id, batch_regno
            )"""

    def _build_batch_details_cte(self, details) -> str:
        cols_without_key, values_sql = sql_utils.prepare_sql_parts(details)
        return f""",
            inserted_batch_details AS (
                INSERT INTO moltrack.batch_details (batch_id, {", ".join(cols_without_key)})
                SELECT ib.id, {", ".join([f"bd.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS bd(batch_regno, {", ".join(cols_without_key)})
                JOIN inserted_batches ib ON bd.batch_regno = ib.batch_regno
            )"""

    def _build_batch_additions_cte(self, additions) -> str:
        cols_without_key, values_sql = sql_utils.prepare_sql_parts(additions)
        return f""",
            inserted_batch_additions AS (
                INSERT INTO moltrack.batch_additions (batch_id, {", ".join(cols_without_key)})
                SELECT ib.id, {", ".join([f"ba.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS ba(batch_regno, {", ".join(cols_without_key)})
                JOIN inserted_batches ib ON ba.batch_regno = ib.batch_regno
            )"""

    def cleanup_chunk(self):
        super().cleanup_chunk()
        self.batches_to_insert.clear()
        self.batch_details.clear()
        self.batch_additions.clear()

    def cleanup(self):
        super().cleanup()
        self.cleanup_chunk()
        self._additions_map = None

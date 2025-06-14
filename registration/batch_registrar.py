from datetime import datetime
import random
from typing import Any, Dict, List, Optional
from fastapi import HTTPException
from pytest import Session
from registration.compound_registrar import CompoundRegistrar
import main
import models


class BatchRegistrar(CompoundRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.batch_records_map = self._load_reference_map(models.Batch, "batch_regno")
        self.additions_map = self._load_reference_map(models.Addition, "name")

        self.batches_to_insert = []
        self.batch_synonyms = []
        self.batch_details = []
        self.batch_additions = []

    def _build_batch_record(self, inchikey: str) -> Dict[str, Any]:
        batch_regno = random.randint(0, 1000)
        if self.batch_records_map.get(batch_regno) is not None:
            raise HTTPException(status_code=400, detail=f"Batch with batch_regno {batch_regno} already exists")

        return {
            "inchikey": inchikey,
            "notes": None,
            "created_by": main.admin_user_id,
            "updated_by": main.admin_user_id,
            "created_at": datetime.now(),
            "batch_regno": batch_regno,
        }

    def _build_batch_addition_record(self, batch_additions: Dict[str, Any], batch_regno: int) -> List[Dict[str, Any]]:
        records = []
        for name, value in batch_additions.items():
            addition = self.additions_map.get(name)
            if addition is None:
                raise HTTPException(status_code=400, detail=f"Unknown addition: {name}")
            records.append(
                {
                    "batch_regno": batch_regno,
                    "addition_id": getattr(addition, "id"),
                    "addition_equivalent": float(value) if value not in (None, "") else 1,
                    "created_by": main.admin_user_id,
                    "updated_by": main.admin_user_id,
                }
            )
        return records

    def get_additional_records(self, grouped, inchikey):
        batch_record = self._build_batch_record(inchikey)
        self.batches_to_insert.append(batch_record)

        self.batch_synonyms.extend(
            self._build_synonym_records(grouped.get("batches_synonyms", {}), batch_record["batch_regno"], "batch_regno")
        )
        self.batch_details.extend(
            self._build_details_records(
                grouped.get("batches_properties", {}), batch_record["batch_regno"], "batch_regno"
            )
        )
        self.batch_additions.extend(
            self._build_batch_addition_record(grouped.get("batches_additions", {}), batch_record["batch_regno"])
        )

    def get_additional_cte(self):
        if not self.batches_to_insert:
            return ""
        return self._build_batch_ctes(
            self.batches_to_insert, self.batch_synonyms, self.batch_details, self.batch_additions
        )

    def _build_batch_ctes(self, batches, synonyms, details, additions) -> str:
        batch_cte = self._build_inserted_batches_cte(batches)
        if synonyms:
            batch_cte += self._build_batch_synonyms_cte(synonyms)
        if details:
            batch_cte += self._build_batch_details_cte(details)
        if additions:
            batch_cte += self._build_batch_additions_cte(additions)

        return batch_cte

    def _build_inserted_batches_cte(self, batches) -> str:
        cols_without_key, values_sql = self._prepare_sql_parts(batches)
        return f"""
            inserted_batches AS (
                INSERT INTO moltrack.batches (compound_id, {", ".join(cols_without_key)})
                SELECT ic.id, {", ".join([f"b.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS b (inchikey, {", ".join(cols_without_key)})
                JOIN inserted_compounds ic ON b.inchikey = ic.inchikey
                RETURNING id, batch_regno
            )"""

    def _build_batch_synonyms_cte(self, synonyms) -> str:
        cols_without_key, values_sql = self._prepare_sql_parts(synonyms)
        return f""",
            inserted_batch_synonyms AS (
                INSERT INTO moltrack.batch_synonyms (batch_id, {", ".join(cols_without_key)})
                SELECT ib.id, {", ".join([f"bs.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS bs(batch_regno, {", ".join(cols_without_key)})
                JOIN inserted_batches ib ON bs.batch_regno = ib.batch_regno
            )"""

    def _build_batch_details_cte(self, details) -> str:
        cols_without_key, values_sql = self._prepare_sql_parts(details)
        return f""",
            inserted_batch_details AS (
                INSERT INTO moltrack.batch_details (batch_id, {", ".join(cols_without_key)})
                SELECT ib.id, {", ".join([f"bd.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS bd(batch_regno, {", ".join(cols_without_key)})
                JOIN inserted_batches ib ON bd.batch_regno = ib.batch_regno
            )"""

    def _build_batch_additions_cte(self, additions) -> str:
        cols_without_key, values_sql = self._prepare_sql_parts(additions)
        return f""",
            inserted_batch_additions AS (
                INSERT INTO moltrack.batch_additions (batch_id, {", ".join(cols_without_key)})
                SELECT ib.id, {", ".join([f"ba.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS ba(batch_regno, {", ".join(cols_without_key)})
                JOIN inserted_batches ib ON ba.batch_regno = ib.batch_regno
            )"""

    def get_additional_output_info(self, grouped) -> dict:
        if not self.batches_to_insert:
            return {}
        last_batch = self.batches_to_insert[-1]
        subset = {k: last_batch[k] for k in ("notes", "batch_regno")}
        return {
            **subset,
            **{f"batch_synonym_{k}": v for k, v in grouped.get("batches_synonyms", {}).items()},
            **{f"batch_property_{k}": v for k, v in grouped.get("batches_properties", {}).items()},
            **{f"batch_addition_{k}": v for k, v in grouped.get("batches_additions", {}).items()},
        }

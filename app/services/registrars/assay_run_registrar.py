from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy.orm import Session

from app import main
from app import models
from app.utils import enums, sql_utils
from app.services.registrars.base_registrar import BaseRegistrar


class AssayRunRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.assay_records_map = self._load_reference_map(models.Assay, "name")
        self.assay_runs_to_insert = []
        self.output_records: List[Dict[str, Any]] = []

    def _build_assay_run_record(self, assay_data: Dict[str, Any], assay_details: Dict[str, Any]) -> Dict[str, Any]:
        assay_name = assay_data.get("name")
        existing_assay = self.assay_records_map.get(assay_name)
        if not existing_assay:
            raise HTTPException(status_code=400, detail=f"Assay {assay_name} not found.")

        # TODO: Clarify and implement the rules for creating the 'name' attribute in assay_run entries
        return {
            "assay_id": getattr(existing_assay, "id"),
            "name": assay_name + assay_details.get("Assay Run Date"),
            "created_by": main.admin_user_id,
            "updated_by": main.admin_user_id,
        }

    def build_sql(self, rows: List[Dict[str, Any]], batch_size: int = 5000):
        global_idx = 0
        for batch in sql_utils.chunked(rows, batch_size):
            self.assay_runs_to_insert = []
            details = []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row, "assay")
                    assay_data = grouped.get("assay", {})
                    assay_run = self._build_assay_run_record(assay_data, grouped.get("assay_run_details"))
                    self.assay_runs_to_insert.append(assay_run)

                    inserted, updated = self.property_service.build_details_records(
                        models.AssayRunDetail,
                        grouped.get("assay_run_details", {}),
                        {"rn": idx + 1},
                        enums.ScopeClass.ASSAY_RUN,
                        False,
                    )
                    details.extend(inserted)
                    self._add_output_row(assay_run, grouped, "success")
                except Exception as e:
                    self.handle_row_error(row, e, global_idx, rows)
                global_idx += 1

            if self.assay_runs_to_insert:
                batch_sql = self.generate_sql(self.assay_runs_to_insert, details)
                self.sql_statements.append(batch_sql)

    def generate_sql(self, assay_runs, details) -> str:
        assay_runs_sql = self._generate_assay_run_sql(assay_runs)
        details_sql = self._generate_details_sql(details)
        return sql_utils.generate_sql(assay_runs_sql, details_sql)

    # TODO: Think of a more robust key than row number for joining
    def _generate_assay_run_sql(self, assay_runs) -> str:
        cols = list(assay_runs[0].keys())
        values_sql = sql_utils.values_sql(assay_runs, cols)
        return f"""
            WITH inserted_assay_runs AS (
                INSERT INTO moltrack.assay_runs ({", ".join(cols)})
                VALUES {values_sql}
                RETURNING id
            ),
            numbered_assay_runs AS (
                SELECT id, ROW_NUMBER() OVER (ORDER BY id) AS rn
                FROM inserted_assay_runs
            )"""

    def _generate_details_sql(self, details) -> str:
        if not details:
            return ""

        cols_without_key, values_sql = sql_utils.prepare_sql_parts(details)
        return f"""
            inserted_assay_run_details AS (
                INSERT INTO moltrack.assay_run_details (assay_run_id, {", ".join(cols_without_key)})
                SELECT nr.id, {", ".join([f"d.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS d(rn, {", ".join(cols_without_key)})
                JOIN numbered_assay_runs nr ON d.rn = nr.rn
            )"""

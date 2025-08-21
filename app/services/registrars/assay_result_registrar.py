from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy.orm import Session
from sqlalchemy import and_, func, or_, select

from app import models
from app.utils.admin_utils import admin
from app.utils import enums, sql_utils
from app.services.registrars.base_registrar import BaseRegistrar


class AssayResultsRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[Dict[str, str]], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.assay_results_to_insert = []
        self.entity_type = enums.EntityType.ASSAY_RESULT

    def _check_single_result(self, results: list, error_context: str):
        if len(results) == 0:
            raise HTTPException(status_code=400, detail=f"No {error_context} found matching the provided details")
        if len(results) > 1:
            raise HTTPException(status_code=400, detail=f"Multiple {error_context} found matching the provided details")
        return results[0]

    def _lookup_by_details(
        self,
        details: Dict[str, Any],
        details_model,
        parent_id_field: str,
        entity_type: str,
        parent_id_value: Optional[int],
    ):
        """
        Generic lookup for detail-based models, e.g. batch_details or assay_run_details.

        :param details: dict of property_name -> value
        :param details_model: model (BatchDetail or AssayRunDetail)
        :param parent_id_field: str, e.g. 'batch_id' or 'assay_run_id'
        :param parent_id_value: int or None, the parent entity id to filter details by;
                                if None, don't filter by parent id (optional)
        :return: subquery returning matching parent ids
        """
        if not details:
            raise HTTPException(status_code=400, detail="No details provided for lookup")

        property_values = []
        for prop_name, value in details.items():
            if value in (None, "", []):
                continue  # Skip this property if value is empty (we are not inserting empty values)
            prop_info = self.property_service.get_property_info(prop_name, entity_type)

            property_values.append(
                {
                    "property_id": getattr(prop_info["property"], "id"),
                    "value_column_name": prop_info["field_name"],
                    "value": value,
                }
            )

        # Construct OR conditions, where each condition ensures that a row matches a specific property and value
        num_details = len(property_values)
        or_conditions = []
        for pv in property_values:
            col = getattr(details_model, pv["value_column_name"])
            conditions = [
                details_model.property_id == pv["property_id"],
                col == pv["value"],
            ]
            if parent_id_value is not None:
                conditions.append(getattr(details_model, parent_id_field) == parent_id_value)
            cond = and_(*conditions)
            or_conditions.append(cond)

        subq = (
            self.db.query(getattr(details_model, parent_id_field))
            .filter(or_(*or_conditions))  # Match any of the provided pairs
            .group_by(getattr(details_model, parent_id_field))  # Group by parent id
            .having(
                func.count(getattr(details_model, parent_id_field)) == num_details
            )  # Keep only those with all matches
            .subquery()
        )
        return subq

    def _lookup_batch_by_details(self, batch_details: Dict[str, Any]) -> models.Batch:
        subq = self._lookup_by_details(batch_details, models.BatchDetail, "batch_id", enums.EntityType.BATCH, None)

        batch_matches = self.db.query(models.Batch).filter(models.Batch.id.in_(select(subq))).all()
        return self._check_single_result(batch_matches, "batches")

    def _lookup_assay_run_by_details(
        self, assay_filter: Dict[str, Any], assay_run_details: Dict[str, Any]
    ) -> models.AssayRun:
        assay_query = self.db.query(models.Assay)
        for col_name, val in assay_filter.items():
            col = getattr(models.Assay, col_name, None)
            if col is None:
                continue
            assay_query = assay_query.filter(col == val)

        assays = assay_query.all()
        assay = self._check_single_result(assays, "assays")
        subq = self._lookup_by_details(
            assay_run_details, models.AssayRunDetail, "assay_run_id", enums.EntityType.ASSAY_RUN, None
        )
        assay_runs = (
            self.db.query(models.AssayRun)
            .filter(models.AssayRun.assay_id == assay.id)
            .filter(models.AssayRun.id.in_(select(subq)))
            .all()
        )
        return self._check_single_result(assay_runs, "assay runs")

    def _build_assay_result_record(self, batch_id: int, assay_run_id: int) -> Dict[str, Any]:
        return {
            "batch_id": batch_id,
            "assay_run_id": assay_run_id,
            "created_by": admin.admin_user_id,
            "updated_by": admin.admin_user_id,
        }

    # TODO: Identify the specific data row(s) in assay_results.csv causing failures
    def build_sql(self, rows: List[Dict[str, Any]]) -> str:
        self.assay_results_to_insert = []
        details = []

        for idx, row in enumerate(rows):

            def process_row(row):
                grouped = self._group_data(row, "assay")
                batch_record = self._lookup_batch_by_details(grouped.get("batch_details"))
                assay_run_record = self._lookup_assay_run_by_details(
                    grouped.get("assay"), grouped.get("assay_run_details")
                )

                batch_id = getattr(batch_record, "id")
                assay_run_id = getattr(assay_run_record, "id")
                assay_result = self._build_assay_result_record(batch_id, assay_run_id)

                inserted, updated = self.property_service.build_details_records(
                    models.AssayResultDetail,
                    grouped.get("assay_result_details", {}),
                    {"rn": idx + 1},
                    enums.EntityType.ASSAY_RESULT,
                    False,
                )

                self.assay_results_to_insert.append(assay_result)
                details.extend(inserted)

            self._process_row(row, process_row)

        if self.assay_results_to_insert:
            batch_sql = self.generate_sql(self.assay_results_to_insert, details)
            details.clear()
            return batch_sql

    def generate_sql(self, assay_results, details) -> str:
        assay_results_sql = self._generate_assay_result_sql(assay_results)
        details_sql = self._generate_details_sql(details)
        return sql_utils.generate_sql(assay_results_sql, details_sql)

    def _generate_assay_result_sql(self, assay_results) -> str:
        cols = list(assay_results[0].keys())
        values_sql = sql_utils.values_sql(assay_results, cols)
        return f"""
            WITH inserted_assay_results AS (
                INSERT INTO moltrack.assay_results ({", ".join(cols)})
                VALUES {values_sql}
                RETURNING id
            ),
            numbered_assay_results AS (
                SELECT id, ROW_NUMBER() OVER (ORDER BY id) AS rn
                FROM inserted_assay_results
            )"""

    def _generate_details_sql(self, details) -> str:
        if not details:
            return ""

        cols_without_key, values_sql = sql_utils.prepare_sql_parts(details)
        return f"""
            inserted_assay_result_details AS (
                INSERT INTO moltrack.assay_result_details (assay_result_id, {", ".join(cols_without_key)})
                SELECT nr.id, {", ".join([f"d.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS d(rn, {", ".join(cols_without_key)})
                JOIN numbered_assay_results nr ON d.rn = nr.rn
            )"""

    def cleanup_chunk(self):
        super().cleanup_chunk()
        self.assay_results_to_insert.clear()

    def cleanup(self):
        super().cleanup()
        self.cleanup_chunk()

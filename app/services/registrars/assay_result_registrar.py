from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy.orm import Session
from sqlalchemy import and_, func, or_

from app import models
from app.utils import enums, sql_utils
from app.services.registrars.base_registrar import BaseRegistrar


class AssayResultsRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[Dict[str, str]], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.output_records: List[Dict[str, Any]] = []
        self.assay_type_records_map = self._load_reference_map(models.Assay, "name")
        self.assay_results_to_insert = []

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
        scope: str,
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
            prop_info = self.property_service.get_property_info(prop_name, scope)

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
        subq = self._lookup_by_details(batch_details, models.BatchDetail, "batch_id", enums.ScopeClass.BATCH, None)

        batch_matches = self.db.query(models.Batch).filter(models.Batch.id.in_(subq)).all()
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
            assay_run_details, models.AssayRunDetail, "assay_run_id", enums.ScopeClass.ASSAY_RUN, None
        )
        assay_runs = (
            self.db.query(models.AssayRun)
            .filter(models.AssayRun.assay_id == assay.id)
            .filter(models.AssayRun.id.in_(subq))
            .all()
        )
        return self._check_single_result(assay_runs, "assay runs")

    # TODO: Identify the specific data row(s) in assay_results.csv causing failures
    def build_sql(self, rows: List[Dict[str, Any]], batch_size: int = 10):
        global_idx = 0
        for batch in sql_utils.chunked(rows, batch_size):
            self.assay_results_to_insert = []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row, "assay")
                    batch_record = self._lookup_batch_by_details(grouped.get("batch_details"))
                    assay_run_record = self._lookup_assay_run_by_details(
                        grouped.get("assay"), grouped.get("assay_run_details")
                    )

                    inserted, updated = self.property_service.build_details_records(
                        models.AssayResult,
                        grouped.get("assay_results", {}),
                        {"batch_id": getattr(batch_record, "id"), "assay_run_id": getattr(assay_run_record, "id")},
                        enums.ScopeClass.ASSAY_RESULT,
                        False,
                    )
                    self.assay_results_to_insert.extend(inserted)
                    self._add_output_row(row, grouped, "success")
                except Exception as e:
                    self.handle_row_error(row, e, global_idx, rows)
                global_idx += 1

            if self.assay_results_to_insert:
                batch_sql = self.generate_sql(self.assay_results_to_insert)
                self.sql_statements.append(batch_sql)

    def generate_sql(self, assay_results) -> str:
        details_sql = self._generate_details_sql(assay_results)
        return sql_utils.generate_sql(details_sql, terminate_with_select=False)

    def _generate_details_sql(self, details) -> str:
        if not details:
            return ""

        cols_without_key, values_sql = sql_utils.prepare_sql_parts(details)
        return f"""
            INSERT INTO moltrack.assay_results (batch_id, {", ".join(cols_without_key)})
            VALUES {values_sql}
        """

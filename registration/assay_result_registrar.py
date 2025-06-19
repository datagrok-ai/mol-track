from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy.orm import Session
from sqlalchemy import and_, func, or_

import models
from registration.base_registrar import BaseRegistrar


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
        parent_id_value: Optional[int],
    ):
        if not details:
            raise HTTPException(status_code=400, detail="No details provided for lookup")

        property_values = []
        for prop_name, value in details.items():
            prop_info = self.property_service.get_property_info(prop_name)

            property_values.append(
                {
                    "property_id": getattr(prop_info["property"], "id"),
                    "value_column_name": prop_info["field_name"],
                    "value": value,
                }
            )

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
            .filter(or_(*or_conditions))
            .group_by(getattr(details_model, parent_id_field))
            .having(func.count(getattr(details_model, parent_id_field)) == num_details)
            .subquery()
        )
        return subq

    def _lookup_batch_by_details(self, batch_details: Dict[str, Any]) -> models.Batch:
        subq = self._lookup_by_details(batch_details, models.BatchDetail, "batch_id", None)

        batch_matches = self.db.query(models.Batch).filter(models.Batch.id.in_(subq)).all()
        return self._check_single_result(batch_matches, "batches")

    def _lookup_assay_run_by_details(
        self, assay_filter: Dict[str, Any], assay_run_details: Dict[str, Any]
    ) -> models.AssayRun:
        assay_query = self.db.query(models.Assay)
        for col_name, val in assay_filter.items():
            col = getattr(models.Assay, col_name, None)
            if col is None:
                raise HTTPException(status_code=400, detail=f"Invalid assay filter column: {col_name}")
            assay_query = assay_query.filter(col == val)

        assays = assay_query.all()
        assay = self._check_single_result(assays, "assays")
        subq = self._lookup_by_details(assay_run_details, models.AssayRunDetail, "assay_run_id", None)

        assay_runs = (
            self.db.query(models.AssayRun)
            .filter(models.AssayRun.assay_id == assay.id)
            .filter(models.AssayRun.id.in_(subq))
            .all()
        )
        return self._check_single_result(assay_runs, "assay runs")

    def build_sql(self, rows: List[Dict[str, Any]], batch_size: int = 5000):
        global_idx = 0
        for batch in self.sql_service.chunked(rows, batch_size):
            self.assay_results_to_insert = []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row)
                    batch_record = self._lookup_batch_by_details(grouped.get("batch_details"))
                    assay_run_record = self._lookup_assay_run_by_details(
                        grouped.get("assays"), grouped.get("assay_run_details")
                    )
                    self.assay_results_to_insert.extend(
                        self.property_service.build_details_records(
                            grouped.get("assay_results", {}),
                            {"batch_id": getattr(batch_record, "id"), "assay_run_id": getattr(assay_run_record, "id")},
                            False,
                        )
                    )
                    self._add_output_row(assay_run_record, grouped, "success")
                except Exception as e:
                    self.handle_row_error(row, e, global_idx, rows)
                global_idx += 1

            if self.assay_results_to_insert:
                batch_sql = self.generate_sql(self.assay_results_to_insert)
                self.sql_statements.append(batch_sql)

    def generate_sql(self, assay_results) -> str:
        details_sql = self._generate_details_sql(assay_results)
        return self.sql_service.generate_sql(details_sql, terminate_with_select=False)

    def _generate_details_sql(self, details) -> str:
        if not details:
            return ""

        cols_without_key, values_sql = self.sql_service.prepare_sql_parts(details)
        return f"""
            INSERT INTO moltrack.assay_results (batch_id, {", ".join(cols_without_key)})
            VALUES {values_sql}
        """

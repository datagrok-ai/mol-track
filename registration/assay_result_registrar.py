from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy.orm import Session
from sqlalchemy import and_, func, or_

import enums
import models
from registration.base_registrar import BaseRegistrar

import utils


class AssayResultsRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[Dict[str, str]], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.output_records: List[Dict[str, Any]] = []
        self.assay_type_records_map = self._load_reference_map(models.Assay, "name")
        self.assay_results_to_insert = []

    def _lookup_by_details(
        self,
        details: Dict[str, Any],
        details_model,
        parent_id_field: str,
        parent_id_value: Optional[int],
    ):
        """
        Generic lookup for detail-based models, e.g. batch_details or assay_run_details.

        :param details: dict of property_name -> value
        :param details_model: ORM model (BatchDetail or AssayRunDetail)
        :param parent_id_field: str, e.g. 'batch_id' or 'assay_run_id'
        :param parent_id_value: int or None, the parent entity ID to filter details by;
                                if None, don't filter by parent ID (optional)
        :return: subquery returning matching parent IDs
        """
        if not details:
            raise HTTPException(status_code=400, detail="No details provided for lookup")

        property_values = []
        for prop_name, value in details.items():
            prop = self.property_records_map.get(prop_name)
            if prop is None:
                raise HTTPException(status_code=400, detail=f"Unknown property: {prop_name}")

            value_type = getattr(prop, "value_type", None)
            if value_type not in utils.value_type_to_field or value_type not in utils.value_type_cast_map:
                raise HTTPException(
                    status_code=400,
                    detail=f"Unsupported or unknown value type for property: {prop_name}",
                )

            value_column_name = utils.value_type_to_field[value_type]
            property_values.append(
                {
                    "property_id": getattr(prop, "id"),
                    "value_column_name": value_column_name,
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
        if len(batch_matches) == 0:
            raise HTTPException(status_code=400, detail=f"No batch found matching batch details: {batch_details}")
        if len(batch_matches) > 1:
            raise HTTPException(
                status_code=400,
                detail=f"Multiple batches found matching batch details: {batch_details}",
            )
        return batch_matches[0]

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
        if len(assays) == 0:
            raise HTTPException(status_code=400, detail=f"No assays found matching: {assay_filter}")
        if len(assays) > 1:
            raise HTTPException(status_code=400, detail=f"Multiple assays found matching: {assay_filter}")

        assay = assays[0]

        subq = self._lookup_by_details(assay_run_details, models.AssayRunDetail, "assay_run_id", None)

        assay_runs = (
            self.db.query(models.AssayRun)
            .filter(models.AssayRun.assay_id == assay.id)
            .filter(models.AssayRun.id.in_(subq))
            .all()
        )

        if len(assay_runs) == 0:
            raise HTTPException(status_code=400, detail=f"No assay runs found matching details: {assay_run_details}")
        if len(assay_runs) > 1:
            raise HTTPException(
                status_code=400, detail=f"Multiple assay runs found matching details: {assay_run_details}"
            )

        return assay_runs[0]

    def build_sql(self, rows: List[Dict[str, Any]], batch_size: int = 5000):
        def chunked(lst, size):
            for i in range(0, len(lst), size):
                yield lst[i : i + size]

        global_idx = 0
        for batch in chunked(rows, batch_size):
            self.assay_results_to_insert = []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row)
                    batch_record = self._lookup_batch_by_details(grouped.get("batch_details"))
                    assay_run_record = self._lookup_assay_run_by_details(
                        grouped.get("assays"), grouped.get("assay_run_details")
                    )
                    self.assay_results_to_insert.extend(
                        self._build_details_records(
                            grouped.get("assay_results", {}),
                            {"batch_id": getattr(batch_record, "id"), "assay_run_id": getattr(assay_run_record, "id")},
                            False,
                        )
                    )

                except Exception as e:
                    self._add_output_row(row, {}, "failed", str(e))
                    if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
                        remaining_rows = rows[global_idx + 1 :]
                        for remaining_row in remaining_rows:
                            self._add_output_row(remaining_row, {}, "not_processed")
                        raise HTTPException(status_code=400, detail=self.result())

                global_idx += 1

            if self.assay_results_to_insert:
                batch_sql = self.generate_sql(self.assay_results_to_insert)
                self.sql_statements.append(batch_sql)

    def generate_sql(self, assay_results) -> str:
        details_sql = self._generate_details_sql(assay_results)
        if not details_sql:
            return ""
        return f"{details_sql};"

    def _generate_details_sql(self, details) -> str:
        if not details:
            return ""

        cols_without_key, values_sql = self._prepare_sql_parts(details)
        return f"""
            INSERT INTO moltrack.assay_results ({", ".join(cols_without_key)})
            VALUES {values_sql}
        """

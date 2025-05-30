from datetime import datetime
import random
from typing import Any, Dict, List, Optional
from fastapi import HTTPException
from pytest import Session
from sqlalchemy import text
from compound_registrar import CompoundRegistrar
import enums
import main


class BatchRegistrar(CompoundRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)

    def _build_batch_record(self, inchikey: str) -> Dict[str, Any]:
        return {
            "inchikey": inchikey,
            "notes": None,
            "created_by": main.admin_user_id,
            "updated_by": main.admin_user_id,
            "created_at": datetime.now(),
            "batch_regno": random.randint(0, 1000),
        }

    def build_batch_synonyms_records(self, batch_id: int, synonyms: Dict[str, str]) -> List[Dict[str, Any]]:
        records = []
        for synonym_name, value in synonyms.items():
            norm_synonym_name = self.normalize_key(synonym_name)
            type_id = self.synonym_type_map.get(norm_synonym_name)
            if type_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown synonym type: {synonym_name}")

            records.append(
                {
                    "batch_id": batch_id,
                    "synonym_type_id": type_id,
                    "synonym_value": value,
                    "created_by": main.admin_user_id,
                    "updated_by": main.admin_user_id,
                }
            )
            return records

    def build_batch_sql(self, rows: List[Dict[str, Any]]):
        super().build_sql(rows)
        batch_records_to_insert, batches_synonyms_to_insert = [], []
        try:
            for idx, row in enumerate(rows):
                try:
                    grouped = self._group_data(row)
                    batch_synonyms = grouped.get("batches_synonyms", {})

                    batch_record = self._build_batch_record(self.compounds_to_insert[idx]["inchikey"])
                    batch_records_to_insert.append(batch_record)

                    batches_syn_records = self._build_batch_synonym_records(
                        batch_synonyms, self.compound_record[idx]["inchikey"]
                    )
                    batches_synonyms_to_insert.extend(batches_syn_records)

                    self.successful_rows += 1
                except Exception as e:
                    self.failed_rows += 1
                    self.error_messages.append(f"Row {idx + 1} failed: {str(e)}")
                    if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
                        raise HTTPException(status_code=400, detail=self.result())

            batch_column_names = batch_records_to_insert[0].keys()

            batch_values_sql = self._values_sql(batch_records_to_insert, batch_column_names)

            self.sql += f""",
            inserted_batches AS (INSERT INTO batches (compound_id, notes, created_by, updated_by, created_at, batch_regno)
            SELECT ic.id, b.notes, b.created_by, b.updated_by, b.created_at, b.batch_regno
            FROM (VALUES {batch_values_sql}) AS b (inchikey, notes, created_by, updated_by, created_at, batch_regno)
            JOIN inserted_compounds ic ON b.inchikey = ic.inchikey
            )
            """

        except Exception as e:
            self.db.rollback()
            raise HTTPException(status_code=500, detail=f"Failed to register compounds: {str(e)}")

    def register_all(self, rows: List[Dict[str, Any]]):
        self.build_batch_sql(rows=rows)
        self.sql += """
            SELECT count(*) FROM inserted_compounds;
        """
        self.db.execute(text(self.sql))
        self.db.commit()

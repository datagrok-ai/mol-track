import csv
import io
import json
from typing import Any, Dict, List, Optional
from fastapi import HTTPException
from pytest import Session
from sqlalchemy import select

import crud
import models
import enums


class CompoundRegistrar:
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        self.db = db
        self.error_handling = error_handling
        self.user_mapping = self.parse_mapping(mapping)
        self.synonym_type_map = self.load_synonym_types()
        self.property_records_map = self.load_properties()
        self.successful_rows = 0
        self.failed_rows = 0
        self.error_messages = []

        # TODO: Implement bulk inserts using SQLAlchemy's bulk_insert_mappings
        # self.compound_data = []
        # self.synonym_data = []
        # self.property_data = []

    def normalize_key(self, name: str) -> str:
        return name.strip().lower().replace(" ", "_")

    def parse_mapping(self, mapping: Optional[str]) -> Dict[str, str]:
        try:
            return json.loads(mapping) if mapping else {}
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid JSON for mapping")

    def load_synonym_types(self) -> Dict[str, int]:
        result = self.db.execute(select(models.SynonymType)).all()
        return {self.normalize_key(st.name): st.id for (st,) in result}

    def load_properties(self) -> Dict[str, int]:
        result = self.db.execute(select(models.Property)).all()
        return {self.normalize_key(p.name): p.id for (p,) in result}

    def process_csv(self, csv_content: str) -> List[Dict[str, Any]]:
        rows = list(csv.DictReader(io.StringIO(csv_content)))
        if not rows:
            raise HTTPException(status_code=400, detail="CSV file is empty or invalid")
        self.normalized_mapping = self.user_mapping or {k: k for k in rows[0].keys()}
        return rows

    def process_row(self, row: Dict[str, Any]):
        grouped = self.group_data(row)
        compound = self.create_compound(grouped.get("compounds", {}))
        self.create_synonyms(compound.id, grouped.get("compounds_synonyms", {}))
        self.create_properties(compound.id, grouped.get("properties", {}))

    def group_data(self, row: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        grouped = {}
        for src_key, mapped_key in self.normalized_mapping.items():
            value = row.get(src_key)
            table, field = mapped_key.split(".", 1) if "." in mapped_key else ("compounds", mapped_key)
            grouped.setdefault(table, {})[field] = value
        return grouped

    def create_compound(self, data: Dict[str, Any]) -> models.Compound:
        return crud.create_compound(self.db, compound=models.CompoundCreate(**data))

    def create_synonyms(self, compound_id: int, synonyms: Dict[str, str]):
        for synonym_name, value in synonyms.items():
            # TODO: Improve lookup
            norm_synonym_name = self.normalize_key(synonym_name)
            type_id = self.synonym_type_map.get(norm_synonym_name)
            if type_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown synonym type: {synonym_name}")
            synonym_input = models.CompoundSynonymBase(
                compound_id=compound_id, synonym_type_id=type_id, synonym_value=value
            )
            crud.create_compound_synonym(self.db, compound_synonym=synonym_input)

    def create_properties(self, compound_id: int, props: Dict[str, Any]):
        for property_name, value in props.items():
            norm_property_name = self.normalize_key(property_name)
            prop_id = self.property_records_map.get(norm_property_name)
            if prop_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown property: {property_name}")

            compound_detail_input = models.CompoundDetailCreate(
                compound_id=compound_id, property_id=prop_id, value=value
            )
            crud.create_compound_detail(self.db, compound_detail=compound_detail_input)

    def register_all(self, rows: List[Dict[str, Any]]):
        for idx, row in enumerate(rows):
            try:
                self.process_row(row)
                self.successful_rows += 1
            except Exception as e:
                self.db.rollback()
                self.failed_rows += 1
                msg = f"Row {idx + 1} failed: {str(e)}"
                self.error_messages.append(msg)

                if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
                    if self.failed_rows > 0:
                        raise HTTPException(status_code=400, detail=self.result())

        self.db.commit()

    def result(self) -> Dict[str, Any]:
        return {
            "status_message": "Partial success" if self.failed_rows else "Success",
            "successful_rows": self.successful_rows,
            "failed_rows": self.failed_rows,
            "errors": self.error_messages if self.failed_rows else [],
        }

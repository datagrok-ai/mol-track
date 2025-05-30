import csv
import io
import json
import uuid
from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy import select, text
from sqlalchemy.orm import Session
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import enums
import main
import models


class CompoundRegistrar:
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        self.db = db
        self.error_handling = error_handling
        self.user_mapping = self._parse_mapping(mapping)
        self.synonym_type_map = self._load_synonym_types()
        self.property_records_map = self._load_properties()
        self.successful_rows = 0
        self.failed_rows = 0
        self.error_messages = []
        self.compounds_to_insert = []
        self.sql = None

    @staticmethod
    def _normalize_key(name: str) -> str:
        return name.strip().lower().replace(" ", "_")

    def _parse_mapping(self, mapping: Optional[str]) -> Dict[str, str]:
        try:
            return json.loads(mapping) if mapping else {}
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid JSON for mapping")

    def _load_synonym_types(self) -> Dict[str, int]:
        result = self.db.execute(select(models.SynonymType)).all()
        return {self._normalize_key(st.name): st.id for (st,) in result}

    def _load_properties(self) -> Dict[str, int]:
        result = self.db.execute(select(models.Property)).all()
        return {self._normalize_key(p.name): p.id for (p,) in result}

    def process_csv(self, csv_content: str) -> List[Dict[str, Any]]:
        rows = list(csv.DictReader(io.StringIO(csv_content)))
        if not rows:
            raise HTTPException(status_code=400, detail="CSV file is empty or invalid")
        self.normalized_mapping = self.user_mapping or {k: k for k in rows[0]}
        return rows

    def _group_data(self, row: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        grouped = {}
        for src_key, mapped_key in self.normalized_mapping.items():
            value = row.get(src_key)
            table, field = mapped_key.split(".", 1) if "." in mapped_key else ("compounds", mapped_key)
            grouped.setdefault(table, {})[field] = value
        return grouped

    def _build_compound_record(self, compound_data: Dict[str, Any], idx: int) -> Dict[str, Any]:
        mol = Chem.MolFromSmiles(compound_data.get("smiles"))
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi)

        return {
            "canonical_smiles": canonical_smiles,
            "inchi": inchi,
            "inchikey": inchikey,
            "original_molfile": compound_data.get("original_molfile", ""),
            "molregno": idx,
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "hash_mol": uuid.uuid4(),
            "hash_tautomer": uuid.uuid4(),
            "hash_canonical_smiles": uuid.uuid4(),
            "hash_no_stereo_smiles": uuid.uuid4(),
            "hash_no_stereo_tautomer": uuid.uuid4(),
            "created_at": datetime.utcnow(),
            "updated_at": datetime.utcnow(),
            "created_by": main.admin_user_id,
            "updated_by": main.admin_user_id,
            "is_archived": compound_data.get("is_archived", False),
        }

    def _build_synonym_records(self, synonyms: Dict[str, Any], inchikey: str) -> List[Dict[str, Any]]:
        records = []
        for syn_name, syn_value in synonyms.items():
            norm_name = self._normalize_key(syn_name)
            type_id = self.synonym_type_map.get(norm_name)
            if type_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown synonym type: {syn_name}")
            records.append(
                {
                    "inchikey": inchikey,
                    "synonym_type_id": type_id,
                    "synonym_value": syn_value,
                    "created_by": main.admin_user_id,
                    "updated_by": main.admin_user_id,
                }
            )
        return records

    def _build_details_records(self, properties: Dict[str, Any], inchikey: str) -> List[Dict[str, Any]]:
        records = []
        for prop_name, value in properties.items():
            norm_prop_name = self._normalize_key(prop_name)
            prop_id = self.property_records_map.get(norm_prop_name)
            if prop_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown property: {prop_name}")

            # property_record = self.db.query(models.Property).filter(models.Property.id == prop_id).first()
            # if not property_record:
            #     raise HTTPException(status_code=404, detail=f"Property with ID {prop_id} not found")

            detail = {
                "inchikey": inchikey,
                "property_id": prop_id,
                "value_string": None,
                "value_num": None,
                "value_datetime": datetime.now(),
                "created_by": main.admin_user_id,
                "updated_by": main.admin_user_id,
            }

            value_type_map = {
                "datetime": "value_datetime",
                "int": "value_num",
                "double": "value_num",
                "string": "value_string",
            }

            attr = value_type_map.get("double")
            if attr and value is not None:
                detail[attr] = float(value)
            else:
                raise HTTPException(status_code=400, detail=f"Unsupported or null value for: {prop_name}")

            records.append(detail)
        return records

    def _values_sql(self, table_data: List[Dict[str, Any]], columns: List[str]) -> str:
        def quote(val):
            if val is None:
                return "NULL"
            if isinstance(val, str):
                return f"'{val.replace("'", "''")}'"
            if isinstance(val, datetime):
                return f"'{val.isoformat()}'::timestamp with time zone"
            if isinstance(val, uuid.UUID):
                return f"'{val}'::uuid"
            return str(val)

        return ",\n".join(f"({', '.join(quote(row.get(col)) for col in columns)})" for row in table_data)

    def result(self) -> Dict[str, Any]:
        return {
            "status_message": "Partial success" if self.failed_rows else "Success",
            "successful_rows": self.successful_rows,
            "failed_rows": self.failed_rows,
            "errors": self.error_messages if self.failed_rows else [],
        }

    def build_sql(self, rows: List[Dict[str, Any]]):
        synonyms_to_insert, details_to_insert = [], []

        try:
            for idx, row in enumerate(rows):
                try:
                    grouped = self._group_data(row)
                    compound_data = grouped.get("compounds", {})
                    synonyms = grouped.get("compounds_synonyms", {})
                    properties = grouped.get("properties", {})

                    compound_record = self._build_compound_record(compound_data, idx)
                    self.compounds_to_insert.append(compound_record)

                    syn_records = self._build_synonym_records(synonyms, compound_record["inchikey"])
                    synonyms_to_insert.extend(syn_records)

                    detail_records = self._build_details_records(properties, compound_record["inchikey"])
                    details_to_insert.extend(detail_records)

                    self.successful_rows += 1
                except Exception as e:
                    self.failed_rows += 1
                    self.error_messages.append(f"Row {idx + 1} failed: {str(e)}")
                    if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
                        raise HTTPException(status_code=400, detail=self.result())

            if not self.compounds_to_insert:
                return

            # compound_column_names = [column.name for column in models.Compound.__table__.columns]
            compound_column_names = self.compounds_to_insert[0].keys()
            synonym_column_names = synonyms_to_insert[0].keys()
            detail_column_names = details_to_insert[0].keys()

            compounds_sql = self._values_sql(self.compounds_to_insert, compound_column_names)
            synonyms_sql = self._values_sql(synonyms_to_insert, synonym_column_names)
            details_sql = self._values_sql(details_to_insert, detail_column_names)

            self.sql = f"""
            WITH inserted_compounds AS (
                INSERT INTO compounds ({", ".join(compound_column_names)})
                VALUES
                {compounds_sql}
                RETURNING id, inchikey
            ),
            inserted_synonyms AS (
                INSERT INTO compound_synonyms (compound_id, synonym_type_id, synonym_value, created_by, updated_by)
                SELECT ic.id, s.synonym_type_id, s.synonym_value, s.created_by, s.updated_by
                FROM (VALUES {synonyms_sql}) AS s(inchikey, synonym_type_id, synonym_value, created_by, updated_by)
                JOIN inserted_compounds ic ON s.inchikey = ic.inchikey
            ),
            inserted_details AS (
                INSERT INTO compound_details (compound_id, property_id, value_string, value_num, value_datetime, created_by, updated_by)
                SELECT ic.id, d.property_id, d.value_string, d.value_num, d.value_datetime, d.created_by, d.updated_by
                FROM (VALUES {details_sql}) AS d(inchikey, property_id, value_string, value_num, value_datetime, created_by, updated_by)
                JOIN inserted_compounds ic ON d.inchikey = ic.inchikey
            )
            """
        except Exception as e:
            self.db.rollback()
            raise HTTPException(status_code=500, detail=f"Failed to register compounds: {str(e)}")

    def register_all(self, rows: List[Dict[str, Any]]):
        self.build_sql(rows=rows)
        self.sql += """
            SELECT count(*) FROM inserted_compounds;
        """
        self.db.execute(text(self.sql))
        self.db.commit()

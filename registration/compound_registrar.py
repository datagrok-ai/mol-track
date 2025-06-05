import uuid
from datetime import datetime
from typing import Any, Callable, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy.orm import Session
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import enums
import main
import models
from registration.base_registrar import BaseRegistrar


class CompoundRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.synonym_type_map = self._load_reference_map(models.SynonymType, "name")
        self.property_records_map = self._load_reference_map(models.Property, "name")
        self.compound_records_map = self._load_reference_map(models.Compound, "inchikey")
        self.compounds_to_insert = []
        self.output_records: List[Dict[str, Any]] = []

    def _build_compound_record(self, compound_data: Dict[str, Any], idx: int) -> Dict[str, Any]:
        mol = Chem.MolFromSmiles(compound_data.get("smiles"))
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        inchikey = Chem.InchiToInchiKey(Chem.MolToInchi(mol))

        if self.compound_records_map.get(inchikey) is not None:
            raise HTTPException(status_code=400, detail=f"Compound with InChIKey {inchikey} already exists")

        now = datetime.utcnow()
        return {
            "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True),
            "inchi": Chem.MolToInchi(mol),
            "inchikey": inchikey,
            "original_molfile": compound_data.get("original_molfile", ""),
            "molregno": idx,
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            # TODO: Replace this once hash merging logic is merged
            **{
                f"hash_{key}": uuid.uuid4()
                for key in ["mol", "tautomer", "canonical_smiles", "no_stereo_smiles", "no_stereo_tautomer"]
            },
            "created_at": now,
            "updated_at": now,
            "created_by": main.admin_user_id,
            "updated_by": main.admin_user_id,
            "is_archived": compound_data.get("is_archived", False),
        }

    def _build_synonym_records(self, synonyms: Dict[str, Any], inchikey: str) -> List[Dict[str, Any]]:
        records = []
        for name, value in synonyms.items():
            type = self.synonym_type_map.get(name)
            if type is None:
                raise HTTPException(status_code=400, detail=f"Unknown synonym type: {name}")
            records.append(
                {
                    "inchikey": inchikey,
                    "synonym_type_id": getattr(type, "id"),
                    "synonym_value": value,
                    "created_by": main.admin_user_id,
                    "updated_by": main.admin_user_id,
                }
            )
        return records

    def _build_details_records(self, properties: Dict[str, Any], inchikey: str) -> List[Dict[str, Any]]:
        records = []
        value_type_to_field = {
            "datetime": "value_datetime",
            "int": "value_num",
            "double": "value_num",
            "string": "value_string",
        }

        def cast_datetime(v):
            return v if isinstance(v, datetime) else datetime.fromisoformat(str(v))

        value_type_cast_map: Dict[str, Callable[[Any], Any]] = {
            "datetime": cast_datetime,
            "int": int,
            "double": float,
            "string": str,
        }

        for prop_name, value in properties.items():
            prop = self.property_records_map.get(prop_name)

            if prop is None:
                raise HTTPException(status_code=400, detail=f"Unknown property: {prop_name}")

            value_type = getattr(prop, "value_type", None)
            if value_type not in value_type_to_field or value_type not in value_type_cast_map:
                raise HTTPException(
                    status_code=400, detail=f"Unsupported or unknown value type for property: {prop_name}"
                )

            if value is None:
                raise HTTPException(status_code=400, detail=f"Null value for property: {prop_name}")

            field_name = value_type_to_field[value_type]
            cast_fn = value_type_cast_map[value_type]
            detail = {
                "inchikey": inchikey,
                "property_id": getattr(prop, "id"),
                "value_string": None,
                "value_num": None,
                "value_datetime": None,
                "created_by": main.admin_user_id,
                "updated_by": main.admin_user_id,
            }

            try:
                detail[field_name] = cast_fn(value)
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Error casting value for property {prop_name}: {e}")

            records.append(detail)
        return records

    def build_sql(self, rows: List[Dict[str, Any]]):
        synonyms, details = [], []

        for idx, row in enumerate(rows):
            try:
                grouped = self._group_data(row)
                compound_data = grouped.get("compounds", {})
                compound = self._build_compound_record(compound_data, idx)
                self.compounds_to_insert.append(compound)

                synonyms.extend(
                    self._build_synonym_records(grouped.get("compounds_synonyms", {}), compound["inchikey"])
                )
                details.extend(self._build_details_records(grouped.get("properties", {}), compound["inchikey"]))
                self._add_output_row(compound_data, grouped, "success")

            except Exception as e:
                self._add_output_row(row, {}, "failed", str(e))
                if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
                    raise HTTPException(status_code=400, detail=self.result())

        if not self.compounds_to_insert:
            return

        self.sql = self.generate_sql(self.compounds_to_insert, synonyms, details)

    def generate_sql(self, compounds, synonyms, details) -> str:
        c_cols = list(compounds[0].keys())
        s_cols = list(synonyms[0].keys()) if synonyms else []
        d_cols = list(details[0].keys()) if details else []

        c_sql = self._values_sql(compounds, c_cols)
        s_sql = self._values_sql(synonyms, s_cols) if synonyms else ""
        d_sql = self._values_sql(details, d_cols) if details else ""

        s_cols = s_cols[1:]
        d_cols = d_cols[1:]

        combined_sql = f"""
            WITH inserted_compounds AS (
                INSERT INTO compounds ({", ".join(c_cols)})
                VALUES {c_sql}
                RETURNING id, inchikey
            )"""

        if len(s_cols) > 0:
            combined_sql += f""",
            inserted_synonyms AS (
                INSERT INTO compound_synonyms (compound_id, {", ".join(s_cols)})
                SELECT ic.id, s.synonym_type_id, s.synonym_value, s.created_by, s.updated_by
                FROM (VALUES {s_sql}) AS s(inchikey, {", ".join(s_cols)})
                JOIN inserted_compounds ic ON s.inchikey = ic.inchikey
            )"""

        if len(d_cols) > 0:
            combined_sql += f""",
            inserted_details AS (
                INSERT INTO compound_details (compound_id, {", ".join(d_cols)})
                SELECT ic.id, d.property_id, d.value_string, d.value_num, d.value_datetime, d.created_by, d.updated_by
                FROM (VALUES {d_sql}) AS d(inchikey, {", ".join(d_cols)})
                JOIN inserted_compounds ic ON d.inchikey = ic.inchikey
            )"""

        return combined_sql

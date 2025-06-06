import random
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

    def _build_compound_record(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
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
            "molregno": random.randint(0, 100),
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

    def _build_synonym_records(self, synonyms: Dict[str, Any], entity_id: Any, id_field: str) -> List[Dict[str, Any]]:
        records = []
        for name, value in synonyms.items():
            type = self.synonym_type_map.get(name)
            if type is None:
                raise HTTPException(status_code=400, detail=f"Unknown synonym type: {name}")
            records.append(
                {
                    id_field: entity_id,
                    "synonym_type_id": getattr(type, "id"),
                    "synonym_value": value,
                    "created_by": main.admin_user_id,
                    "updated_by": main.admin_user_id,
                }
            )
        return records

    def _build_details_records(self, properties: Dict[str, Any], entity_id: Any, id_field: str) -> List[Dict[str, Any]]:
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
                id_field: entity_id,
                "property_id": getattr(prop, "id"),
                "created_by": main.admin_user_id,
                "updated_by": main.admin_user_id,
            }

            try:
                detail[field_name] = cast_fn(value)
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Error casting value for property {prop_name}: {e}")

            records.append(detail)
        return records

    def build_sql(self, rows: List[Dict[str, Any]], batch_size: int = 5000):
        def chunked(lst, size):
            for i in range(0, len(lst), size):
                yield lst[i : i + size]

        global_idx = 0
        for batch in chunked(rows, batch_size):
            self.compounds_to_insert = []
            synonyms, details = [], []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row)
                    compound_data = grouped.get("compounds", {})
                    compound = self._build_compound_record(compound_data)
                    self.compounds_to_insert.append(compound)

                    synonyms.extend(
                        self._build_synonym_records(
                            grouped.get("compounds_synonyms", {}), compound["inchikey"], "inchikey"
                        )
                    )
                    details.extend(
                        self._build_details_records(grouped.get("properties", {}), compound["inchikey"], "inchikey")
                    )
                    self._add_output_row(compound_data, grouped, "success")
                except Exception as e:
                    self._add_output_row(row, {}, "failed", str(e))
                    if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
                        remaining_rows = rows[global_idx + 1 :]
                        for remaining_row in remaining_rows:
                            self._add_output_row(remaining_row, {}, "not_processed")
                        raise HTTPException(status_code=400, detail=self.result())
                global_idx += 1

            if self.compounds_to_insert:
                # TODO: Improve error handling
                extra_sql = self.get_additional_cte(self.compounds_to_insert, rows)
                batch_sql = self.generate_sql(self.compounds_to_insert, synonyms, details, extra_sql)
                self.sql_statements.append(batch_sql)

    def generate_sql(self, compounds, synonyms, details, extra_sql) -> str:
        compound_sql = self._generate_compound_sql(compounds)
        synonym_sql = self._generate_synonym_sql(synonyms)
        details_sql = self._generate_details_sql(details)

        combined_sql = compound_sql
        if synonym_sql:
            combined_sql += ",\n" + synonym_sql
        if details_sql:
            combined_sql += ",\n" + details_sql
        if extra_sql:
            combined_sql += ",\n" + extra_sql

        combined_sql += "\nSELECT 1;"
        return combined_sql

    def _generate_compound_sql(self, compounds) -> str:
        cols = list(compounds[0].keys())
        values_sql = self._values_sql(compounds, cols)
        return f"""
            WITH inserted_compounds AS (
                INSERT INTO moltrack.compounds ({", ".join(cols)})
                VALUES {values_sql}
                RETURNING id, inchikey
            )"""

    def _generate_synonym_sql(self, synonyms) -> str:
        if not synonyms:
            return ""

        cols_without_key, values_sql = self._prepare_sql_parts(synonyms)

        return f"""
            inserted_synonyms AS (
                INSERT INTO moltrack.compound_synonyms (compound_id, {", ".join(cols_without_key)})
                SELECT ic.id, {", ".join([f"s.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS s(inchikey, {", ".join(cols_without_key)})
                JOIN inserted_compounds ic ON s.inchikey = ic.inchikey
            )"""

    def _generate_details_sql(self, details) -> str:
        if not details:
            return ""

        cols_without_key, values_sql = self._prepare_sql_parts(details)
        return f"""
            inserted_details AS (
                INSERT INTO moltrack.compound_details (compound_id, {", ".join(cols_without_key)})
                SELECT ic.id, {", ".join([f"d.{col}" for col in cols_without_key])}
                FROM (VALUES {values_sql}) AS d(inchikey, {", ".join(cols_without_key)})
                JOIN inserted_compounds ic ON d.inchikey = ic.inchikey
            )"""

    def get_additional_cte(self, compounds, rows):
        pass

    def _prepare_sql_parts(self, records: List[Dict[str, Any]]):
        cols = list(records[0].keys())
        key, *cols_without_key = cols
        values_sql = self._values_sql(records, cols)
        return cols_without_key, values_sql

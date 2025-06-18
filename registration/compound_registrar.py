from datetime import datetime
from typing import Any, Callable, Dict, List, Optional, Tuple

from fastapi import HTTPException
from sqlalchemy import func
from sqlalchemy.orm import Session
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chemistry_utils import generate_hash_layers, generate_uuid_from_string, standardize_mol
from rdkit.Chem.RegistrationHash import HashLayer, GetMolHash
import enums
import main
import models
from registration.base_registrar import BaseRegistrar


class CompoundRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.property_records_map = self._load_reference_map(models.Property, "name")
        self.compound_records_map = self._load_reference_map(models.Compound, "inchikey")
        self.compound_details_map = self._load_reference_map(models.CompoundDetail, "id")
        self.compounds_to_insert = []
        self.output_records: List[Dict[str, Any]] = []
        self._molregno_counter = (db.query(func.max(models.Compound.molregno)).scalar() or 0) + 1

    def _next_molregno(self) -> int:
        molregno = self._molregno_counter
        self._molregno_counter += 1
        return molregno

    def _build_compound_record(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
        mol = Chem.MolFromSmiles(compound_data.get("smiles"))
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        inchikey = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
        existing_compound = self.compound_records_map.get(inchikey)
        if existing_compound is not None:
            return self.model_to_dict(existing_compound)

        now = datetime.utcnow()
        standarized_mol = standardize_mol(mol)
        mol_layers = generate_hash_layers(standarized_mol)
        hash_mol = GetMolHash(mol_layers)
        canonical_smiles = mol_layers[HashLayer.CANONICAL_SMILES]
        hash_canonical_smiles = generate_uuid_from_string(mol_layers[HashLayer.CANONICAL_SMILES])
        hash_tautomer = generate_uuid_from_string(mol_layers[HashLayer.TAUTOMER_HASH])
        hash_no_stereo_smiles = generate_uuid_from_string(mol_layers[HashLayer.NO_STEREO_SMILES])
        hash_no_stereo_tautomer = generate_uuid_from_string(mol_layers[HashLayer.NO_STEREO_TAUTOMER_HASH])

        return {
            "canonical_smiles": canonical_smiles,
            "inchi": Chem.MolToInchi(mol),
            "inchikey": inchikey,
            "original_molfile": compound_data.get("original_molfile", ""),
            "molregno": self._next_molregno(),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "hash_mol": hash_mol,
            "hash_tautomer": hash_tautomer,
            "hash_canonical_smiles": hash_canonical_smiles,
            "hash_no_stereo_smiles": hash_no_stereo_smiles,
            "hash_no_stereo_tautomer": hash_no_stereo_tautomer,
            "created_at": now,
            "updated_at": now,
            "created_by": main.admin_user_id,
            "updated_by": main.admin_user_id,
            "is_archived": compound_data.get("is_archived", False),
        }

    def _build_details_records(
        self, properties: Dict[str, Any], entity_id: Any, id_field: str
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        records_to_insert, records_to_update = [], []

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

            field_name = value_type_to_field[value_type]
            cast_fn = value_type_cast_map[value_type]

            detail = {
                id_field: entity_id,
                "property_id": getattr(prop, "id"),
                "created_by": main.admin_user_id,
                "updated_by": main.admin_user_id,
                "value_datetime": datetime.now(),
                "value_num": 0,
                "value_string": None,
            }

            try:
                detail[field_name] = cast_fn(value)
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Error casting value for property {prop_name}: {e}")

            compound = self.compound_records_map.get(entity_id)
            if compound:
                compound_id = getattr(compound, "id")
                prop_id = getattr(prop, "id")

                for detail_id, compound_detail in self.compound_details_map.items():
                    detail_dict = self.model_to_dict(compound_detail)
                    if detail_dict["compound_id"] == compound_id and detail_dict["property_id"] == prop_id:
                        detail = {
                            ("compound_id" if k == id_field else k): (compound_id if k == id_field else v)
                            for k, v in detail.items()
                        }
                        records_to_update.append(detail)
            else:
                records_to_insert.append(detail)

        return records_to_insert, records_to_update

    def build_sql(self, rows: List[Dict[str, Any]], batch_size: int = 5000):
        def chunked(lst, size):
            for i in range(0, len(lst), size):
                yield lst[i : i + size]

        global_idx = 0
        for batch in chunked(rows, batch_size):
            self.compounds_to_insert = []
            details_to_insert, details_to_update = [], []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row)
                    compound_data = grouped.get("compounds", {})
                    compound = self._build_compound_record(compound_data)
                    self.compounds_to_insert.append(compound)

                    insert, update = self._build_details_records(
                        grouped.get("compounds_details", {}), compound["inchikey"], "inchikey"
                    )
                    details_to_insert.extend(insert)
                    details_to_update.extend(update)

                    self.get_additional_records(grouped, compound["inchikey"])
                    self._add_output_row(compound_data, grouped, "success")
                except Exception as e:
                    self._add_output_row(row, {}, "failed", str(e))
                    if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
                        remaining_rows = rows[global_idx + 1 :]
                        for remaining_row in remaining_rows:
                            self._add_output_row(remaining_row, {}, "not_processed")
                        raise HTTPException(status_code=400, detail=self.result())
                global_idx += 1

            extra_sql = self.get_additional_cte()
            batch_sql = self.generate_sql(self.compounds_to_insert, details_to_insert, details_to_update, extra_sql)
            self.sql_statements.append(batch_sql)

    def generate_sql(self, compounds, details_to_insert, details_to_update, extra_sql) -> str:
        parts = []
        compound_sql = self._generate_compound_sql(compounds)
        if compound_sql:
            parts.append(compound_sql)

        details_to_insert_sql = self._generate_details_sql(details_to_insert)
        if details_to_insert_sql:
            parts.append(details_to_insert_sql)

        details_to_update_sql = self._generate_details_update_sql(details_to_update)
        if details_to_update_sql:
            parts.append(details_to_update_sql)

        if extra_sql:
            parts.append(extra_sql)

        if parts:
            combined_sql = "WITH " + ",\n".join(parts)
            combined_sql += "\nSELECT 1;"
        else:
            combined_sql = "SELECT 1;"

        return combined_sql

    def _generate_compound_sql(self, compounds) -> str:
        if not compounds:
            return ""

        cols = list(compounds[0].keys())
        values_sql = self._values_sql(compounds, cols)
        insert_cte = f"""
            inserted_compounds AS (
                INSERT INTO moltrack.compounds ({", ".join(cols)})
                VALUES {values_sql}
                ON CONFLICT (inchikey) DO NOTHING
                RETURNING id, inchikey
            ),
        """

        inchikeys = [f"'{c['inchikey']}'" for c in compounds]
        inchikey_list = ", ".join(inchikeys)
        available_cte = f"""
            available_compounds AS (
                SELECT id, inchikey FROM inserted_compounds
                UNION
                SELECT id, inchikey FROM moltrack.compounds WHERE inchikey IN ({inchikey_list})
            )
        """
        return insert_cte + available_cte

    def _generate_details_update_sql(self, details) -> str:
        if not details:
            return ""

        cols = ["compound_id", "property_id", "value_datetime", "value_num", "value_string", "updated_by"]
        vals = self._values_sql(details, cols)
        return f"""updated_details AS (
            UPDATE moltrack.compound_details cd
            SET value_datetime = v.value_datetime, value_num = v.value_num, value_string = v.value_string, updated_by = v.updated_by
            FROM (VALUES {vals}) AS v(compound_id, property_id, value_datetime, value_num, value_string, updated_by)
            WHERE cd.compound_id = v.compound_id
            AND cd.property_id = v.property_id
            RETURNING cd.*
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
                JOIN available_compounds ic ON d.inchikey = ic.inchikey
            )"""

    def get_additional_cte(self):
        pass

    def get_additional_records(self, grouped, inchikey):
        pass

    def _prepare_sql_parts(self, records: List[Dict[str, Any]]):
        cols = list(records[0].keys())
        key, *cols_without_key = cols
        values_sql = self._values_sql(records, cols)
        return cols_without_key, values_sql

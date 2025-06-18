import random
from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
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
            "molregno": random.randint(0, 100000),
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

    def build_sql(self, rows: List[Dict[str, Any]], batch_size: int = 5000):
        def chunked(lst, size):
            for i in range(0, len(lst), size):
                yield lst[i : i + size]

        global_idx = 0
        for batch in chunked(rows, batch_size):
            self.compounds_to_insert = []
            details = []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row)
                    compound_data = grouped.get("compounds", {})
                    compound = self._build_compound_record(compound_data)
                    self.compounds_to_insert.append(compound)

                    details.extend(
                        self._build_details_records(
                            grouped.get("compounds_details", {}), {"inchikey": compound["inchikey"]}
                        )
                    )
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

            if self.compounds_to_insert:
                extra_sql = self.get_additional_cte()
                batch_sql = self.generate_sql(self.compounds_to_insert, details, extra_sql)
                self.sql_statements.append(batch_sql)

    def generate_sql(self, compounds, details, extra_sql) -> str:
        compound_sql = self._generate_compound_sql(compounds)
        details_sql = self._generate_details_sql(details)

        combined_sql = compound_sql
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

    def get_additional_cte(self):
        pass

    def get_additional_records(self, grouped, inchikey):
        pass

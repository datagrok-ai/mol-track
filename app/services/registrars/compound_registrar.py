from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy import func
from sqlalchemy.orm import Session
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from app.utils.chemistry_utils import generate_hash_layers, generate_uuid_from_string, standardize_mol
from rdkit.Chem.RegistrationHash import HashLayer, GetMolHash
import main
from app import models
from app.utils import enums, sql_utils
from app.services.registrars.base_registrar import BaseRegistrar


class CompoundRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.compound_records_map = self._load_reference_map(models.Compound, "inchikey")
        self.compound_details_map = self._load_reference_map(models.CompoundDetail, "id")
        self.compounds_to_insert = []
        self.output_records: List[Dict[str, Any]] = []

    def _next_molregno(self) -> int:
        db_max = self.db.query(func.max(models.Compound.molregno)).scalar() or 0
        local_max = max((c.get("molregno", 0) for c in self.compounds_to_insert), default=0)
        return max(db_max, local_max) + 1

    def _build_compound_record(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
        mol = Chem.MolFromSmiles(compound_data.get("smiles"))
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        inchikey = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
        existing_compound = self.compound_records_map.get(inchikey)
        if existing_compound is not None:
            compound_dict = self.model_to_dict(existing_compound)
            compound_dict.pop("id", None)
            return compound_dict

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

    def _compound_update_checker(self, entity_ids, detail, field_name, new_value: Any) -> models.UpdateCheckResult:
        id_field, entity_id = next(iter(entity_ids.items()))
        compound = self.compound_records_map.get(entity_id)
        if not compound:
            return models.UpdateCheckResult(action="insert")

        compound_id = getattr(compound, "id")
        prop_id = detail["property_id"]
        for compound_detail in self.compound_details_map.values():
            detail_dict = self.model_to_dict(compound_detail)
            if detail_dict["compound_id"] == compound_id and detail_dict["property_id"] == prop_id:
                if detail_dict.get(field_name) != new_value:
                    update_data = {
                        ("compound_id" if k == id_field else k): (compound_id if k == id_field else v)
                        for k, v in detail.items()
                    }
                    return models.UpdateCheckResult(action="update", update_data=update_data)
                else:
                    return models.UpdateCheckResult(action="skip")
        return models.UpdateCheckResult(action="insert")

    def build_sql(self, rows: List[Dict[str, Any]], batch_size: int = 5000):
        global_idx = 0
        for batch in sql_utils.chunked(rows, batch_size):
            self.compounds_to_insert = []
            details_to_insert, details_to_update = [], []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row)
                    compound_data = grouped.get("compounds", {})
                    compound = self._build_compound_record(compound_data)
                    self.compounds_to_insert.append(compound)

                    inserted, updated = self.property_service.build_details_records(
                        models.CompoundDetail,
                        grouped.get("compounds_details", {}),
                        {"inchikey": compound["inchikey"]},
                        enums.ScopeClass.COMPOUND,
                        True,
                        self._compound_update_checker,
                    )

                    details_to_insert.extend(inserted)
                    details_to_update.extend(updated)

                    self.get_additional_records(grouped, compound["inchikey"])
                    self._add_output_row(compound_data, grouped, "success")
                except Exception as e:
                    self.handle_row_error(row, e, global_idx, rows)
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
        values_sql = sql_utils.values_sql(compounds, cols)
        insert_cte = f"""
            inserted_compounds AS (
                INSERT INTO moltrack.compounds ({", ".join(cols)})
                VALUES {values_sql}
                ON CONFLICT (molregno) DO NOTHING
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

    def _generate_details_update_sql(self, details: List[Dict[str, Any]]) -> str:
        if not details:
            return ""

        required_cols = ["compound_id", "property_id", "updated_by"]
        value_cols = {key for detail in details for key in detail if key.startswith("value_")}
        all_cols = required_cols + sorted(value_cols)
        set_clauses = [f"{col} = v.{col}" for col in sorted(value_cols)] + ["updated_by = v.updated_by"]
        set_clause_sql = ", ".join(set_clauses)
        alias_cols_sql = ", ".join(all_cols)
        vals_sql = sql_utils.values_sql(details, all_cols)

        return f"""updated_details AS (
            UPDATE moltrack.compound_details cd
            SET {set_clause_sql}
            FROM (VALUES {vals_sql}) AS v({alias_cols_sql})
            WHERE cd.compound_id = v.compound_id
            AND cd.property_id = v.property_id
            RETURNING cd.*
        )"""

    def _generate_details_sql(self, details) -> str:
        if not details:
            return ""

        cols_without_key, values_sql = sql_utils.prepare_sql_parts(details)
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

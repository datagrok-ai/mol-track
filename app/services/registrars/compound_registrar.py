from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy.orm import Session
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from rdkit.Chem.RegistrationHash import HashLayer, GetMolHash, HashScheme
from app import main
from app import models
from app.utils import enums, sql_utils, chemistry_utils
from app.utils.logging_utils import logger
from app.services.registrars.base_registrar import BaseRegistrar
from sqlalchemy.sql import text


class CompoundRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self.compound_records_map = self._load_reference_map(models.Compound, "hash_mol")
        self.compound_details_map = self._load_reference_map(models.CompoundDetail, "id")
        self.compounds_to_insert: Dict[str, Dict[str, Any]] = {}

        self.output_records: List[Dict[str, Any]] = []
        self.matching_setting = self._load_matching_setting()

    def _next_molregno(self) -> int:
        molregno = self.db.execute(text("SELECT nextval('moltrack.molregno_seq')")).scalar()
        return molregno

    def _load_matching_setting(self) -> HashScheme:
        try:
            setting = self.db.execute(
                text("SELECT value FROM moltrack.settings WHERE name = 'Compound Matching Rule'")
            ).scalar()
            if setting is None:
                return HashScheme.ALL_LAYERS
            return HashScheme[setting]
        except Exception as e:
            logger.error(f"Error loading compound matching setting: {e}")
            return HashScheme.ALL_LAYERS

    def _build_compound_record(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
        smiles = compound_data.get("smiles")
        if not smiles:
            raise HTTPException(status_code=400, detail="SMILES value is required for compound creation.")

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES string")

        standarized_mol = chemistry_utils.standardize_mol(mol)
        mol_layers = chemistry_utils.generate_hash_layers(standarized_mol)
        hash_mol = GetMolHash(mol_layers, self.matching_setting)

        existing_compound = self.compound_records_map.get(hash_mol)
        db = True

        if existing_compound is None:
            existing_compound = self.compounds_to_insert.get(hash_mol)
            db = False
        # TODO: Implement proper uniqueness rules to ensure data integrity
        if existing_compound is not None:
            if not db:
                return existing_compound
            compound_dict = self.model_to_dict(existing_compound)
            compound_dict.pop("id", None)
            return compound_dict

        now = datetime.utcnow()

        inchikey = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
        if inchikey is None:
            raise HTTPException(status_code=400, detail="Failed to generate InChIKey: possibly invalid molecule")

        canonical_smiles = mol_layers[HashLayer.CANONICAL_SMILES]
        hash_canonical_smiles = chemistry_utils.generate_uuid_from_string(mol_layers[HashLayer.CANONICAL_SMILES])
        hash_tautomer = chemistry_utils.generate_uuid_from_string(mol_layers[HashLayer.TAUTOMER_HASH])
        hash_no_stereo_smiles = chemistry_utils.generate_uuid_from_string(mol_layers[HashLayer.NO_STEREO_SMILES])
        hash_no_stereo_tautomer = chemistry_utils.generate_uuid_from_string(
            mol_layers[HashLayer.NO_STEREO_TAUTOMER_HASH]
        )

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
        compound = next(
            (c for c in self.compound_records_map.values() if getattr(c, id_field, None) == entity_id), None
        )
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
            self.compounds_to_insert = {}
            details_to_insert, details_to_update = [], []

            for idx, row in enumerate(batch):
                try:
                    grouped = self._group_data(row)
                    compound_data = grouped.get("compound", {})
                    compound = self._build_compound_record(compound_data)
                    self.compounds_to_insert[compound["hash_mol"]] = compound

                    inserted, updated = self.property_service.build_details_records(
                        models.CompoundDetail,
                        grouped.get("compound_details", {}),
                        {"molregno": compound["molregno"]},
                        enums.ScopeClass.COMPOUND,
                        True,
                        self._compound_update_checker,
                    )

                    details_to_insert.extend(inserted)
                    details_to_update.extend(updated)

                    self.get_additional_records(grouped, compound["molregno"])
                    self._add_output_row(compound_data, grouped, "success")
                except Exception as e:
                    self.handle_row_error(row, e, global_idx, rows)
                global_idx += 1

            extra_sql = self.get_additional_cte()
            all_compounds_list = list(self.compounds_to_insert.values())
            batch_sql = self.generate_sql(all_compounds_list, details_to_insert, details_to_update, extra_sql)
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
                ON CONFLICT (hash_mol) DO NOTHING
                RETURNING id, molregno, hash_mol
            ),
        """

        hash_mols = [f"'{c['hash_mol']}'" for c in compounds]
        hash_mol_list = ", ".join(hash_mols)
        available_cte = f"""
            available_compounds AS (
                SELECT id, molregno, hash_mol FROM inserted_compounds
                UNION
                SELECT id, molregno, hash_mol FROM moltrack.compounds
                WHERE hash_mol IN ({hash_mol_list})
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
                FROM (VALUES {values_sql}) AS d(molregno, {", ".join(cols_without_key)})
                JOIN available_compounds ic ON d.molregno = ic.molregno
            )"""

    def _group_data(self, row: Dict[str, Any], entity_name: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
        grouped = super()._group_data(row, entity_name)
        value = self.property_service.institution_synonym_dict["compound_details"]
        grouped.setdefault("compound_details", {})[value] = None
        return grouped

    def get_additional_cte(self):
        pass

    def get_additional_records(self, grouped, molregno):
        pass

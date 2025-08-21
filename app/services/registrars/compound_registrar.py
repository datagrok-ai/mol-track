from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy.orm import Session
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from rdkit.Chem.RegistrationHash import HashLayer, GetMolHash, HashScheme
from app.utils.admin_utils import admin
from app import models
from app.utils import enums, sql_utils, chemistry_utils
from app.utils.logging_utils import logger
from app.services.registrars.base_registrar import BaseRegistrar
from sqlalchemy.sql import text


class CompoundRegistrar(BaseRegistrar):
    def __init__(self, db: Session, mapping: Optional[str], error_handling: str):
        super().__init__(db, mapping, error_handling)
        self._compound_records_map = None
        self._compound_details_map = None

        self.compounds_to_insert: Dict[str, Dict[str, Any]] = {}
        self.matching_setting = self._load_matching_setting()
        self.normalized_mapping = {}

        self.entity_type = enums.EntityType.COMPOUND

    @property
    def compound_records_map(self) -> Dict[str, models.Compound]:
        if self._compound_records_map is None:
            self._compound_records_map = self._load_reference_map(models.Compound, "hash_mol")
        return self._compound_records_map

    @property
    def compound_details_map(self) -> Dict[int, models.CompoundDetail]:
        if self._compound_details_map is None:
            self._compound_details_map = self._load_reference_map(models.CompoundDetail, "id")
        return self._compound_details_map

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

        mol = chemistry_utils.validate_rdkit_call(
            Chem.MolFromSmiles, smiles, err_msg_prefix=f"Invalid SMILES '{smiles}':"
        )
        standardized_mol = chemistry_utils.standardize_mol(mol, self.db)
        mol_layers = chemistry_utils.generate_hash_layers(standardized_mol)
        hash_mol = GetMolHash(mol_layers, self.matching_setting)

        existing_compound = self.compound_records_map.get(hash_mol)
        if existing_compound is None:
            existing_compound = self.compounds_to_insert.get(hash_mol)
            if existing_compound:
                return existing_compound
        else:
            compound_dict = self.model_to_dict(existing_compound)
            compound_dict.pop("id", None)
            return compound_dict

        now = datetime.now()

        inchi = chemistry_utils.validate_rdkit_call(Chem.MolToInchi, mol, err_msg_prefix="Failed to generate InChI:")
        inchikey = chemistry_utils.validate_rdkit_call(
            Chem.InchiToInchiKey, inchi, err_msg_prefix="Failed to generate InChIKey:"
        )
        canonical_smiles = mol_layers[HashLayer.CANONICAL_SMILES]
        hash_canonical_smiles = chemistry_utils.generate_uuid_from_string(mol_layers[HashLayer.CANONICAL_SMILES])
        hash_tautomer = chemistry_utils.generate_uuid_from_string(mol_layers[HashLayer.TAUTOMER_HASH])
        hash_no_stereo_smiles = chemistry_utils.generate_uuid_from_string(mol_layers[HashLayer.NO_STEREO_SMILES])
        hash_no_stereo_tautomer = chemistry_utils.generate_uuid_from_string(
            mol_layers[HashLayer.NO_STEREO_TAUTOMER_HASH]
        )

        compound = {
            "canonical_smiles": canonical_smiles,
            "inchi": inchi,
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
            "created_by": admin.admin_user_id,
            "updated_by": admin.admin_user_id,
            "is_archived": compound_data.get("is_archived", False),
        }

        del mol
        del standardized_mol

        return compound

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

    def build_sql(self, rows: List[Dict[str, Any]]) -> str:
        self.compounds_to_insert = {}
        details_to_insert, details_to_update = [], []

        for idx, row in enumerate(rows):

            def process_row(row):
                grouped = self._group_data(row)
                compound_data = grouped.get("compound", {})
                compound = self._build_compound_record(compound_data)
                molregno = compound["molregno"]

                # This step is performed here specifically to attach corporate IDs to the output row
                self.inject_corporate_property(row, grouped, molregno, enums.EntityType.COMPOUND)
                inserted, updated = self.property_service.build_details_records(
                    models.CompoundDetail,
                    grouped.get("compound_details", {}),
                    {"molregno": molregno},
                    enums.EntityType.COMPOUND,
                    True,
                    self._compound_update_checker,
                )

                self.get_additional_records(row, grouped, molregno)

                # Only add the resulting data after it has been fully processed
                # to ensure that no partial or invalid data from this row gets registered.
                self.compounds_to_insert[compound["hash_mol"]] = compound
                details_to_insert.extend(inserted)
                details_to_update.extend(updated)

            self._process_row(row, process_row)

        extra_sql = self.get_additional_cte()
        all_compounds_list = list(self.compounds_to_insert.values())
        batch_sql = self.generate_sql(all_compounds_list, details_to_insert, details_to_update, extra_sql)

        # Clear temporary data structures
        details_to_insert.clear()
        details_to_update.clear()
        return batch_sql

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

    def inject_corporate_property(
        self, row: dict[str, Any], grouped: dict[str, Any], entity_value: str, entity_type: enums.EntityType
    ):
        entity_type_lower = entity_type.value.lower()
        prop_name = f"corporate_{entity_type_lower}_id"
        props = self.property_records_map.get(prop_name, [])

        prop = next((p for p in props if p.entity_type == entity_type), None)
        if not prop:
            return

        value = prop.pattern.format(entity_value)
        row[prop_name] = value
        details_key = f"{entity_type_lower}_details"
        grouped.setdefault(details_key, {})[prop_name] = value

    def cleanup_chunk(self):
        super().cleanup_chunk()
        self.compounds_to_insert.clear()

    def cleanup(self):
        super().cleanup()
        self.cleanup_chunk()
        self._compound_records_map = None
        self._compound_details_map = None
        self.matching_setting = None

    def get_additional_cte(self):
        pass

    def get_additional_records(self, row, grouped, molregno):
        pass

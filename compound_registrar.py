import csv
import io
import json
import uuid
from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import HTTPException
from sqlalchemy import insert, select
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import Session

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import crud
import enums
import main
import models


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

    def group_data(self, row: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        grouped = {}
        for src_key, mapped_key in self.normalized_mapping.items():
            value = row.get(src_key)
            table, field = mapped_key.split(".", 1) if "." in mapped_key else ("compounds", mapped_key)
            grouped.setdefault(table, {})[field] = value
        return grouped

    def bulk_create_compounds(self, compounds_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not compounds_data:
            return []

        prepared_compounds = []
        existing_inchikeys = set()
        existing_smiles = set()

        for data in compounds_data:
            smiles = data.get("smiles")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise HTTPException(status_code=400, detail=f"Invalid SMILES string: {smiles}")

            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            inchi = Chem.MolToInchi(mol)
            inchikey = Chem.InchiToInchiKey(inchi)

            if inchikey in existing_inchikeys:
                raise HTTPException(status_code=400, detail=f"Duplicate InChIKey in batch: {inchikey}")
            if canonical_smiles in existing_smiles:
                raise HTTPException(status_code=400, detail=f"Duplicate canonical SMILES in batch: {canonical_smiles}")

            if crud.get_compound_by_inchi_key(self.db, inchikey):
                raise HTTPException(status_code=400, detail=f"Compound with InChIKey {inchikey} already exists")
            if crud.get_compound_by_canonical_smiles(self.db, canonical_smiles):
                raise HTTPException(
                    status_code=400, detail=f"Compound with canonical SMILES {canonical_smiles} already exists"
                )

            existing_inchikeys.add(inchikey)
            existing_smiles.add(canonical_smiles)

            prepared = {
                "canonical_smiles": canonical_smiles,
                "original_molfile": data.get("original_molfile"),
                "inchi": inchi,
                "inchikey": inchikey,
                "molregno": data.get("molregno", uuid.uuid4().int & (1 << 31) - 1),
                "formula": rdMolDescriptors.CalcMolFormula(mol),
                "hash_mol": uuid.uuid4(),
                "hash_tautomer": uuid.uuid4(),
                "hash_canonical_smiles": uuid.uuid4(),
                "hash_no_stereo_smiles": uuid.uuid4(),
                "hash_no_stereo_tautomer": uuid.uuid4(),
                "created_at": datetime.utcnow(),
                "updated_at": datetime.utcnow(),
                "created_by": data.get("created_by", main.admin_user_id),
                "updated_by": data.get("updated_by", main.admin_user_id),
                "is_archived": data.get("is_archived", False),
            }

            prepared_compounds.append(prepared)

        stmt = insert(models.Compound).values(prepared_compounds).returning(models.Compound.id)
        result = self.db.execute(stmt).fetchall()
        return [{"id": row.id, **prepared_compounds[idx]} for idx, row in enumerate(result)]

    def bulk_create_synonyms(self, synonyms_data: List[Dict[str, Any]]):
        if not synonyms_data:
            return
        for data in synonyms_data:
            data.setdefault("created_by", main.admin_user_id)
            data.setdefault("updated_by", main.admin_user_id)

        stmt = insert(models.CompoundSynonym).values(synonyms_data)
        self.db.execute(stmt)

    def bulk_create_properties(self, properties_data: List[Dict[str, Any]]):
        if not properties_data:
            return
        for data in properties_data:
            data.setdefault("created_by", main.admin_user_id)
            data.setdefault("updated_by", main.admin_user_id)

        stmt = insert(models.CompoundDetail).values(properties_data)
        self.db.execute(stmt)

    def prepare_synonyms_data(self, compound_id: int, synonyms: Dict[str, str]) -> List[Dict[str, Any]]:
        result = []
        for synonym_name, value in synonyms.items():
            norm_synonym_name = self.normalize_key(synonym_name)
            type_id = self.synonym_type_map.get(norm_synonym_name)
            if type_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown synonym type: {synonym_name}")
            result.append(
                {
                    "compound_id": compound_id,
                    "synonym_type_id": type_id,
                    "synonym_value": value,
                    "created_by": main.admin_user_id,
                    "updated_by": main.admin_user_id,
                }
            )
        return result

    def prepare_properties_data(self, compound_id: int, props: Dict[str, Any]) -> List[Dict[str, Any]]:
        result = []
        for property_name, value in props.items():
            norm_property_name = self.normalize_key(property_name)
            prop_id = self.property_records_map.get(norm_property_name)
            if prop_id is None:
                raise HTTPException(status_code=400, detail=f"Unknown property: {property_name}")

            property_record = self.db.query(models.Property).filter(models.Property.id == prop_id).first()
            if not property_record:
                raise HTTPException(status_code=404, detail=f"Property with ID {prop_id} not found")

            detail = {
                "compound_id": compound_id,
                "property_id": prop_id,
                "created_by": main.admin_user_id,
                "updated_by": main.admin_user_id,
            }

            value_type_map = {
                "datetime": "value_datetime",
                "int": "value_num",
                "double": "value_num",
                "string": "value_string",
            }

            attr = value_type_map.get(property_record.value_type)
            if attr:
                if value is not None:
                    detail[attr] = value
            else:
                raise HTTPException(
                    status_code=400, detail=f"Unsupported property value type: {property_record.value_type}"
                )

            result.append(detail)
        return result

    def process_row(self, row: Dict[str, Any]):
        grouped = self.group_data(row)
        compound_data = grouped.get("compounds", {})
        synonyms = grouped.get("compounds_synonyms", {})
        properties = grouped.get("properties", {})

        return compound_data, synonyms, properties

    def register_all(self, rows: List[Dict[str, Any]]):
        compounds_to_insert = []
        synonyms_to_insert = []
        properties_to_insert = []

        try:
            for idx, row in enumerate(rows):
                try:
                    compound_data, synonyms, properties = self.process_row(row)
                    compounds_to_insert.append(compound_data)
                    synonyms_to_insert.append((idx, synonyms))
                    properties_to_insert.append((idx, properties))
                    self.successful_rows += 1
                except Exception as e:
                    self.failed_rows += 1
                    msg = f"Row {idx + 1} failed: {str(e)}"
                    self.error_messages.append(msg)
                    if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
                        raise HTTPException(status_code=400, detail=self.result())

            inserted_compounds = self.bulk_create_compounds(compounds_to_insert)
            idx_to_compound_id = {idx: inserted_compounds[idx]["id"] for idx in range(len(inserted_compounds))}

            bulk_synonyms = []
            for idx, synonyms in synonyms_to_insert:
                compound_id = idx_to_compound_id.get(idx)
                if compound_id is not None and synonyms:
                    bulk_synonyms.extend(self.prepare_synonyms_data(compound_id, synonyms))

            bulk_properties = []
            for idx, properties in properties_to_insert:
                compound_id = idx_to_compound_id.get(idx)
                if compound_id is not None and properties:
                    bulk_properties.extend(self.prepare_properties_data(compound_id, properties))

            self.bulk_create_synonyms(bulk_synonyms)
            self.bulk_create_properties(bulk_properties)
            self.db.commit()

        except SQLAlchemyError as e:
            self.db.rollback()
            raise HTTPException(status_code=500, detail=f"Database error: {str(e)}")

    def result(self) -> Dict[str, Any]:
        return {
            "status_message": "Partial success" if self.failed_rows else "Success",
            "successful_rows": self.successful_rows,
            "failed_rows": self.failed_rows,
            "errors": self.error_messages if self.failed_rows else [],
        }

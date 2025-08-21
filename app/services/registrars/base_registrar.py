import csv
import io
import json
from abc import ABC, abstractmethod
from typing import Iterator, List, Dict, Any, Optional

from sqlalchemy import select, text
from fastapi import HTTPException

from app.utils import enums
from app.utils.logging_utils import logger
from app.services import property_service
from app import models

from rdkit import Chem


class BaseRegistrar(ABC):
    def __init__(self, db, mapping: Optional[str], error_handling: str):
        """
        Base class for processing and registering data to a database.
        :param db: SQLAlchemy database session.
        :param mapping: Optional JSON string defining field mappings.
        :param error_handling: Strategy for handling errors during processing.
        """
        self.db = db
        self.error_handling = error_handling
        self._property_records_map = None
        self._addition_records_map = None

        self.property_service = property_service.PropertyService(self.property_records_map)
        self.user_mapping = self._load_mapping(mapping)
        self.output_rows = []

        self.stop_registration = False
        self.entity_type = None

    @property
    def property_records_map(self):
        if self._property_records_map is None:
            self._property_records_map = self._load_reference_map(models.Property, "name", allow_list=True)
        return self._property_records_map

    @property
    def addition_records_map(self):
        if self._addition_records_map is None:
            self._addition_records_map = self._load_reference_map(models.Addition, "name")
        return self._addition_records_map

    # === Input processing methods ===

    def _load_mapping(self, mapping: Optional[str]) -> Dict[str, str]:
        if not mapping:
            return {}
        try:
            return json.loads(mapping)
        except json.JSONDecodeError:
            raise HTTPException(status_code=400, detail="Invalid JSON for mapping")

    def process_csv(self, file_stream: io.TextIOBase, chunk_size=5000) -> Iterator[List[Dict[str, Any]]]:
        reader = csv.DictReader(file_stream, skipinitialspace=True)

        try:
            first_row = next(reader)
        except StopIteration:
            raise HTTPException(status_code=400, detail="CSV file is empty or invalid")

        if self.user_mapping:
            self.normalized_mapping = self.user_mapping
        else:
            self.normalized_mapping = {}
            for col in first_row.keys():
                assigned = self._assign_column(col)
                self.normalized_mapping[col] = assigned

        chunk = [first_row]

        for row in reader:
            chunk.append(row)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []

        if chunk:
            yield chunk

    def process_sdf(self, file_stream: io.TextIOBase, chunk_size=5000) -> Iterator[List[Dict[str, Any]]]:
        chunk: List[Dict[str, Any]] = []
        current_mol_lines: List[str] = []
        current_props: Dict[str, str] = {}
        prop_name: str | None = None
        mapping_initialized = False

        for line in file_stream:
            line = line.rstrip("\n")

            if line == "$$$$":
                row = dict(current_props)
                molfile_str = "\n".join(current_mol_lines)
                row["original_molfile"] = molfile_str

                # Get smiles
                mol = Chem.MolFromMolBlock(molfile_str)
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                else:
                    smiles = None
                row["smiles"] = smiles

                if not mapping_initialized:
                    if self.user_mapping:
                        self.normalized_mapping = self.user_mapping
                    else:
                        self.normalized_mapping = {k: self._assign_column(k) for k in row.keys()}
                    mapping_initialized = True

                chunk.append(row)

                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []

                current_mol_lines = []
                current_props = {}
                prop_name = None
                continue

            if line.startswith(">  <") and line.endswith(">"):
                prop_name = line[4:-1].strip()
                current_props[prop_name] = ""
                continue

            if prop_name:
                if current_props[prop_name]:
                    current_props[prop_name] += "\n" + line
                else:
                    current_props[prop_name] = line
                continue

            current_mol_lines.append(line)

        if chunk:
            yield chunk

        if not chunk and not current_mol_lines and not current_props:
            raise HTTPException(status_code=400, detail="SDF file is empty or invalid")

    def _assign_column(self, col: str) -> str:
        prefix_base = {
            enums.EntityType.COMPOUND: "compound_details",
            enums.EntityType.BATCH: "batch_details",
            enums.EntityType.ASSAY_RUN: "assay_run_details",
            enums.EntityType.ASSAY_RESULT: "assay_result_details",
        }

        records = self.property_records_map.get(col)
        if records:
            entity_types = {r.entity_type for r in records}
            prefix_key = self.entity_type if self.entity_type in entity_types else next(iter(entity_types), None)
            prefix = prefix_base.get(prefix_key)
            return f"{prefix}.{col}" if prefix else col

        if col in self.addition_records_map:
            return f"batch_additions.{col}"

        return col

    def _group_data(self, row: Dict[str, Any], entity_name: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
        grouped = {}
        for src_key, mapped_key in self.normalized_mapping.items():
            value = row.get(src_key)
            table, field = (
                mapped_key.split(".", 1)
                if "." in mapped_key
                else (entity_name if entity_name else "compound", mapped_key)
            )
            grouped.setdefault(table, {})[field] = value
        return grouped

    # === Reference loading methods ===

    def _load_reference_map(self, model, key: str = "id", allow_list: bool = False):
        result = self.db.execute(select(model)).scalars().all()
        if allow_list:
            reference_map = {}
            for row in result:
                k = getattr(row, key)
                reference_map.setdefault(k, []).append(row)
            return reference_map
        else:
            return {getattr(row, key): row for row in result}

    def model_to_dict(self, obj):
        return {c.key: getattr(obj, c.key) for c in obj.__table__.columns}

    # === SQL construction and registration methods ===

    @abstractmethod
    def build_sql(self, rows: List[Dict[str, Any]]) -> str:
        pass

    @abstractmethod
    def generate_sql(self) -> Optional[str]:
        pass

    def register_all(self, rows: List[Dict[str, Any]]):
        batch_sql = self.build_sql(rows)

        if batch_sql:
            try:
                self.db.execute(text(batch_sql))
                self.db.commit()
            except Exception as e:
                logger.error(f"An exception occurred: {e}")
                self.db.rollback()

    def _process_row(self, row: Dict[str, Any], process_func):
        if self.stop_registration:
            self._add_output_row(row, "not_processed")

        try:
            process_func(row)
            self._add_output_row(row, "success")
        except Exception as e:
            self.handle_row_error(row, e)

    # === Output formatting methods ===

    def _add_output_row(self, row, status, error_msg=None):
        row["registration_status"] = status
        row["registration_error_message"] = error_msg or ""
        self.output_rows.append(row)

    def cleanup_chunk(self):
        self.output_rows.clear()

    def cleanup(self):
        self.cleanup_chunk()
        self.user_mapping.clear()
        self.normalized_mapping.clear()
        self._property_records_map = None
        self._addition_records_map = None
        self.property_service = None

    # === Error handling methods ===

    def handle_row_error(self, row, exception):
        self._add_output_row(row, "failed", str(exception))
        if self.error_handling == enums.ErrorHandlingOptions.reject_all.value:
            self.stop_registration = True

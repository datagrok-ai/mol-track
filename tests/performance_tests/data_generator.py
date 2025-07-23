import random
from typing import Dict, Any, List
from uuid import uuid4, uuid1
from file_writer import FileWriter
from schema_builder import SchemaBuilder
from pathlib import Path
from datetime import datetime, timedelta
from smiles_generator import SmilesGenerator
from schema_rules import SMILES, EPA_BATCH_ID, ASSAY_NAME, ASSAY_RUN_DATE


class DataGenerator:
    value_types = ["string", "int", "double"]

    def __init__(self, size: str, chunk_size: int = 10000, output_dir="output", property_config=None, seed: int = None):
        """
        Initialize DataGenerator with size, chunk_size, output_dir, and property_config.
        If a seed is provided, sets the random seed for reproducibility.
        """
        self.size = size
        self.chunk_size = chunk_size
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.property_config = property_config or {}
        self.total_compounds = self._get_compound_count()
        self.num_records: Dict[str, int] = {}
        self._non_identity_cache: Dict[str, List[Dict[str, Any]]] = {}
        if seed is not None:
            random.seed(seed)

    def generate_all(self) -> Dict[str, int]:
        """
        Generate all data: compounds, batches, assays, assay runs, and assay results.
        Writes each scope to a separate CSV or JSON file.

        Returns:
            Dict[str, int]: A dictionary of record counts per scope.
        """
        sb = SchemaBuilder(self.output_dir, self.value_types)
        self.schemas = sb.generate_schema_and_mapping()
        smiles_list = self._generate_compounds()
        batch_ids = self._generate_batches(smiles_list)
        assay_name = self._generate_assays()
        assay_runs = self._generate_assay_run(assay_name)
        self._generate_assay_results(assay_name, assay_runs, batch_ids)
        return self.num_records

    def _generate_compounds(self) -> List[str]:
        """
        Generate compound records and write to 'compound.csv'.

        Each compound includes a unique SMILES string and synthetic non-identity properties.

        Returns:
            List[str]: A list of generated SMILES strings.
        """
        scope = "compound"
        smiles_list = []
        total_generated = 0
        sg = SmilesGenerator(size=self.total_compounds)
        writer = self._get_writer(scope)

        while total_generated < self.total_compounds:
            chunk = []
            for _ in range(min(self.chunk_size, self.total_compounds - total_generated)):
                smi = sg.generate_smile()
                smiles_list.append(smi)
                details = self._generate_details(scope)
                row = {SMILES: smi, **details}
                chunk.append(row)
                total_generated += 1
            writer.write_csv_chunk(f"{scope}.csv", chunk)

        self.num_records[scope] = total_generated + total_generated * (
            len(self._get_non_identity_properties(scope)) + 1
        )  # +1 for corporate_batch_id
        return smiles_list

    def _generate_batches(self, smiles_list: List[str]) -> List[str]:
        """
        Generate batch records for each compound and write to 'batch.csv'.

        Args:
            smiles_list (List[str]): List of SMILES strings from generated compounds.

        Returns:
            List[str]: A list of generated EPA Batch IDs.
        """
        scope = "batch"
        writer = self._get_writer(scope)
        batch_ids = []
        total_generated = 0

        for chunk_start in range(0, len(smiles_list), self.chunk_size):
            chunk_smiles = smiles_list[chunk_start : chunk_start + self.chunk_size]
            chunk = []
            for smiles in chunk_smiles:
                num_batches = self._get_record_count(scope)
                total_generated += num_batches
                for _ in range(num_batches):
                    epa_batch_id = f"batch_id_{uuid1()}"
                    batch_ids.append(epa_batch_id)
                    details = self._generate_details(scope)
                    row = {SMILES: smiles, EPA_BATCH_ID: epa_batch_id, **details}
                    chunk.append(row)
            writer.write_csv_chunk(f"{scope}.csv", chunk)

        self.num_records[scope] = total_generated + total_generated * (
            len(self._get_non_identity_properties(scope)) + 1
        )  # +1 for batch_compound_id
        return batch_ids

    def _generate_assays(self) -> str:
        """
        Generate a single assay and write it to 'assays.json'.

        Returns:
            str: Name of the generated assay.
        """
        scope = "assay"
        properties = self.schemas[scope]["properties"]
        ar_properties = self.schemas["assay_result"]["properties"]
        assay = [{prop["name"]: self._generate_value(prop["value_type"]) for prop in properties}]
        assay[0][ASSAY_NAME] = f"assay_name_{uuid1()}"
        assay[0]["assay_result_properties"] = [{"name": p["name"], "required": False} for p in ar_properties]

        writer = FileWriter(self.output_dir)
        writer.write_json("assays.json", assay)
        self.num_records[scope] = 1
        return assay[0][ASSAY_NAME]

    def _generate_assay_run(self, assay_name: str) -> List[str]:
        """
        Generate 10 unique assay runs and write them to 'assay_run.csv'.

        Args:
            assay_name (str): The name of the assay to associate with the assay runs.

        Returns:
            List[str]: A list of unique assay run dates.
        """
        scope = "assay_run"
        assay_runs = []
        unique_dates = set()
        assay_run_dates = []
        writer = self._get_writer(scope)
        total_generated = 0

        for _ in range(10):
            total_generated += 1
            assay_run_date = self._generate_value("datetime")
            while assay_run_date in unique_dates:
                assay_run_date = self._generate_value("datetime")
            unique_dates.add(assay_run_date)
            assay_run_dates.append(assay_run_date)

            details = self._generate_details(scope)
            row = {ASSAY_NAME: assay_name, ASSAY_RUN_DATE: assay_run_date, **details}
            assay_runs.append(row)
        writer.write_csv_chunk(f"{scope}.csv", assay_runs)

        self.num_records[scope] = total_generated + total_generated * len(self._get_non_identity_properties(scope))
        return assay_run_dates

    def _generate_assay_results(self, assay_name: str, assay_runs: List[datetime], batch_ids: List[str]) -> None:
        """
        Generate assay results and write them to 'assay_result.csv'.

        Args:
            assay_name (str): The name of the assay to associate with the assay results.
            assay_runs (List[str]): List of unique assay run dates.
            batch_ids (List[str]): List of unique EPA Batch IDs.
        """
        scope = "assay_result"
        writer = self._get_writer(scope)
        total_generated = 0

        for chunk_start in range(0, len(batch_ids), self.chunk_size):
            chunk_batches = batch_ids[chunk_start : chunk_start + self.chunk_size]
            chunk = []
            for batch_id in chunk_batches:
                for ar in assay_runs:
                    num_records = self._get_record_count(scope)
                    total_generated += num_records
                    for _ in range(num_records):
                        details = self._generate_details(scope)
                        row = {ASSAY_NAME: assay_name, ASSAY_RUN_DATE: ar, EPA_BATCH_ID: batch_id, **details}
                        chunk.append(row)
            writer.write_csv_chunk(f"{scope}.csv", chunk)

        self.num_records[scope] = total_generated * len(self._get_non_identity_properties(scope))

    # ==================== Utility Functions ==================== #
    def _get_writer(self, scope: str) -> FileWriter:
        """Returns a writer for the given scope."""
        return FileWriter(self.output_dir, self._get_headers(scope))

    def _get_headers(self, scope: str) -> List[str]:
        """
        Constructs the CSV header for a given schema scope. It excludes
        identity columns and includes only synthetic property names,
        as defined by the schema rules for that scope.

        Args:
            scope (str): The name of the schema scope (e.g., 'compound').

        Returns:
            List[str]: A list of column headers for the CSV row.
        """
        schema = self.schemas[scope]
        identity = schema["identity"]
        properties = schema["properties"]
        return [p["name"] for p in identity + properties]

    def _get_non_identity_properties(self, scope: str) -> List[Dict[str, Any]]:
        """Return and cache non-identity properties for a given scope."""
        if scope in self._non_identity_cache:
            return self._non_identity_cache[scope]
        schema = self.schemas[scope]
        identity_cols = set(schema["identity_cols"])
        non_identity = [prop for prop in schema["properties"] if prop["name"] not in identity_cols]
        self._non_identity_cache[scope] = non_identity
        return non_identity

    def _generate_details(self, scope: str) -> Dict[str, Any]:
        """
        Generates a dictionary of synthetic values for non-identity properties.

        Args:
            scope (str): The schema scope (e.g., 'compound').

        Returns:
            Dict[str, Any]: A dictionary mapping property names to values.
        """
        return {
            prop["name"]: self._generate_value(prop["value_type"]) for prop in self._get_non_identity_properties(scope)
        }

    def _generate_value(self, value_type) -> Any:
        """Generate a mock value for the given type."""
        if value_type in ["float", "double"]:
            return round(random.uniform(0.1, 1000.0), 2)
        elif value_type in ["int", "integer"]:
            return random.randint(0, 1000)
        elif value_type == "string":
            return f"mock_{uuid4().hex[:6]}"
        elif value_type == "datetime":
            return self.random_date_from_past_50_years()
        return None

    def random_date_from_past_50_years(self) -> datetime:
        """Generate a random date in the past 50 years."""

        def random_date(start_date, end_date):
            delta = end_date - start_date
            return start_date + timedelta(days=random.randint(0, delta.days))

        today = datetime.now()
        past_date = today - timedelta(days=50 * 365)  # ~100 years ago
        return random_date(past_date, today)

    def _get_compound_count(self) -> int:
        """Get total number of compounds to generate based on size label."""
        return {"xsmall": 5, "small": 1000, "large": 1_000_000, "xlarge": 100_000_000}[self.size]

    def _get_record_count(self, scope: str) -> int:
        """Get the number of records to generate for a scope, using weighted probabilities."""
        distributions = {
            "batch": {"choices": [1, 2, 3, 4, 5], "weights": [80, 5, 5, 5, 5]},
            "assay_result": {"choices": [2, 3, 4, 5], "weights": [79, 7, 7, 7]},
        }

        if scope not in distributions:
            raise ValueError(f"No distribution defined for scope: {scope}")

        dist = distributions[scope]
        return random.choices(dist["choices"], weights=dist["weights"], k=1)[0]

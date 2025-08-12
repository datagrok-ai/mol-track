import random
from typing import Dict, Any, List, NamedTuple
from uuid import uuid4, uuid1
from tests.performance_tests.file_writer import FileWriter
from tests.performance_tests.schema_builder import SchemaBuilder
from pathlib import Path
from datetime import datetime, timedelta
from tests.performance_tests.smiles_generator import SmilesGenerator
from tests.performance_tests.schema_rules import SMILES, EPA_BATCH_ID, ASSAY_NAME, ASSAY_RUN_DATE
from tests.performance_tests.enums import DatasetSize
from app.utils.enums import ValueType, EntityType


class AssayRunKey(NamedTuple):
    assay_name: str
    assay_run_date: str


class DataGenerator:
    def __init__(self, config: Dict[str, Any], chunk_size: int = 10000, output_dir="output", seed: int = None):
        """
        Initialize DataGenerator with size, chunk_size, output_dir, and property_config.
        If a seed is provided, sets the random seed for reproducibility.
        """

        self.chunk_size = chunk_size
        self.output_dir = Path(output_dir)
        self._create_directory(self.output_dir)
        schema_directory = self.output_dir / "schemas"
        self._create_directory(schema_directory)
        size = DatasetSize(config.get("dataset_size", "small"))
        self.total_compounds = self._get_compound_count(size)
        self.num_records: Dict[str, int] = {}
        self._non_identity_cache: Dict[str, List[Dict[str, Any]]] = {}
        self._create_value_types_for_entity_types()
        self.distributions = config.get("distributions", {})
        self.config = config
        if seed is not None:
            random.seed(seed)

    def generate_all(self) -> Dict[str, int]:
        """
        Generate all data: compounds, batches, assays, assay runs, and assay results.
        Writes each entity_type to a separate CSV or JSON file.

        Returns:
            Dict[str, int]: A dictionary of record counts per entity_type.
        """

        sb = SchemaBuilder(self.output_dir, self.config, self.value_types)
        self.schemas = sb.generate_schema_and_mapping()
        smiles_list = self._generate_compounds()
        batch_ids = self._generate_batches(smiles_list)
        assay_names = self._generate_assays()
        assay_runs_keys = self._generate_assay_runs(assay_names)
        self._generate_assay_results(assay_runs_keys, batch_ids)
        return self.num_records

    def _generate_compounds(self) -> List[str]:
        """
        Generate compound records and write to 'compound.csv'.

        Each compound includes a unique SMILES string and synthetic non-identity properties.

        Returns:
            List[str]: A list of generated SMILES strings.
        """

        entity_type = EntityType.COMPOUND
        smiles_list = []
        sg = SmilesGenerator(data_path=self.output_dir)
        writer = self._get_writer(entity_type)

        while len(smiles_list) < self.total_compounds:
            chunk = []
            for _ in range(min(self.chunk_size, self.total_compounds - len(smiles_list))):
                smi = sg.generate_smile()
                smiles_list.append(smi)
                details = self._generate_details(entity_type)
                row = {SMILES: smi, **details}
                chunk.append(row)
            writer.write_csv_chunk(f"{entity_type.value}.csv", chunk)

        self.num_records[entity_type] = self.total_compounds + self.total_compounds * (
            len(self._get_non_identity_properties(entity_type)) + 1
        )  # +1 for corporate_compound_id
        return smiles_list

    def _generate_batches(self, smiles_list: List[str]) -> List[str]:
        """
        Generate batch records for each compound and write to 'batch.csv'.

        Args:
            smiles_list (List[str]): List of SMILES strings from generated compounds.

        Returns:
            List[str]: A list of generated EPA Batch IDs.
        """

        entity_type = EntityType.BATCH
        writer = self._get_writer(entity_type)
        batch_ids = []
        total_generated = 0

        for chunk_start in range(0, len(smiles_list), self.chunk_size):
            chunk_smiles = smiles_list[chunk_start : chunk_start + self.chunk_size]
            chunk = []
            for smiles in chunk_smiles:
                num_batches = self._get_record_count(entity_type)
                total_generated += num_batches
                for _ in range(num_batches):
                    epa_batch_id = f"batch_id_{uuid4()}"
                    batch_ids.append(epa_batch_id)
                    details = self._generate_details(entity_type)
                    row = {SMILES: smiles, EPA_BATCH_ID: epa_batch_id, **details}
                    chunk.append(row)
            writer.write_csv_chunk(f"{entity_type.value}.csv", chunk)

        self.num_records[entity_type] = total_generated + total_generated * (
            len(self._get_non_identity_properties(entity_type)) + 2
        )  # +2 for corporate_batch_id and EPA Batch ID
        return batch_ids

    def _generate_assays(self) -> List[str]:
        """
        Generate a single assay and write it to 'ASSAY.json'.

        Returns:
            str: Name of the generated assay.
        """

        entity_type = EntityType.ASSAY
        num_of_assays = self.config.get("number_of_assays", 1)
        assays = []
        assay_names = []
        for i in range(num_of_assays):
            properties = self.schemas[entity_type]["properties"]
            ar_properties = self.schemas[EntityType.ASSAY_RESULT]["properties"]
            assays.append({prop["name"]: self._generate_value(prop["value_type"]) for prop in properties})
            assays[i][ASSAY_NAME] = f"assay_name_{uuid1()}"
            assay_names.append(assays[i][ASSAY_NAME])
            assays[i]["assay_result_properties"] = [{"name": p["name"], "required": False} for p in ar_properties]

        writer = FileWriter(self.output_dir)
        writer.write_json(f"{entity_type.value}.json", assays)
        self.num_records[entity_type] = 1
        return assay_names

    def _generate_assay_runs(self, assay_names: List[str]) -> List[AssayRunKey]:
        """
        Generate assay runs and write them to 'assay_run.csv'.

        Args:
            assay_names (List[str]): The name of the assay to associate with the assay runs.

        Returns:
            List[str]: A list of unique assay run dates.
        """

        entity_type = EntityType.ASSAY_RUN
        assay_runs = []
        unique_dates = set()
        assay_run_dates = []
        writer = self._get_writer(entity_type)
        assay_run_keys = []
        total_generated = 0

        for assay_name in assay_names:
            num_assay_runs = self._get_record_count(entity_type)
            total_generated += num_assay_runs
            for _ in range(num_assay_runs):
                assay_run_date = self.generate_unique_date(unique_dates)
                assay_run_dates.append(assay_run_date)

                details = self._generate_details(entity_type)
                row = {ASSAY_NAME: assay_name, ASSAY_RUN_DATE: assay_run_date, **details}
                assay_run_keys.append(AssayRunKey(assay_name, assay_run_date))
                assay_runs.append(row)
        writer.write_csv_chunk(f"{entity_type.value}.csv", assay_runs)

        self.num_records[entity_type] = total_generated + total_generated * len(
            self._get_non_identity_properties(entity_type)
        )
        return assay_run_keys

    def _generate_assay_results(self, assay_run_keys: List[AssayRunKey], batch_ids: List[str]) -> None:
        """
        Generate assay results and write them to 'assay_result.csv'.

        Args:
            assay_name (str): The name of the assay to associate with the assay results.
            assay_runs (List[str]): List of unique assay run dates.
            batch_ids (List[str]): List of unique EPA Batch IDs.
        """

        entity_type = EntityType.ASSAY_RESULT
        writer = self._get_writer(entity_type)
        total_generated = 0

        for chunk_start in range(0, len(batch_ids), self.chunk_size):
            chunk_batches = batch_ids[chunk_start : chunk_start + self.chunk_size]
            chunk = []
            for batch_id in chunk_batches:
                for ark in assay_run_keys:
                    num_records = self._get_record_count(entity_type)
                    total_generated += num_records
                    for _ in range(num_records):
                        details = self._generate_details(entity_type)
                        row = {
                            ASSAY_NAME: ark.assay_name,
                            ASSAY_RUN_DATE: ark.assay_run_date,
                            EPA_BATCH_ID: batch_id,
                            **details,
                        }
                        chunk.append(row)
            writer.write_csv_chunk(f"{entity_type.value}.csv", chunk)

        self.num_records[entity_type] = total_generated + total_generated * len(
            self._get_non_identity_properties(entity_type)
        )

    # ==================== Utility Functions ==================== #
    def _get_writer(self, entity_type: EntityType) -> FileWriter:
        """Returns a writer for the given entity_type."""

        return FileWriter(self.output_dir, self._get_headers(entity_type.value))

    def _get_headers(self, entity_type: EntityType) -> List[str]:
        """
        Constructs the CSV header for a given schema entity_type. It excludes
        identity columns and includes only synthetic property names,
        as defined by the schema rules for that entity_type.

        Args:
            entity_type (str): The name of the schema entity_type (e.g., 'compound').

        Returns:
            List[str]: A list of column headers for the CSV row.
        """

        schema = self.schemas[entity_type]
        identity = schema["identity"]
        properties = schema["properties"]
        return [p["name"] for p in identity + properties]

    def _get_non_identity_properties(self, entity_type: EntityType) -> List[Dict[str, Any]]:
        """Return and cache non-identity properties for a given entity_type."""

        if entity_type in self._non_identity_cache:
            return self._non_identity_cache[entity_type]
        schema = self.schemas[entity_type]
        identity_cols = set(schema["identity_cols"])
        non_identity = [prop for prop in schema["properties"] if prop["name"] not in identity_cols]
        self._non_identity_cache[entity_type] = non_identity
        return non_identity

    def _generate_details(self, entity_type: EntityType) -> Dict[str, Any]:
        """
        Generates a dictionary of synthetic values for non-identity properties.

        Args:
            entity_type (str): The schema entity_type (e.g., 'compound').

        Returns:
            Dict[str, Any]: A dictionary mapping property names to values.
        """

        return {
            prop["name"]: self._generate_value(prop["value_type"])
            for prop in self._get_non_identity_properties(entity_type)
        }

    def _generate_value(self, value_type) -> Any:
        """Generate a mock value for the given type."""

        if value_type == "double":
            return round(random.uniform(0.1, 1000.0), 2)
        elif value_type == "int":
            return random.randint(0, 1000)
        elif value_type == "string":
            return f"mock_{uuid4().hex[:6]}"
        elif value_type == "datetime":
            return str(self.random_date_from_past_50_years())
        elif value_type == "bool":
            return random.choice([True, False])
        elif value_type == "uuid":
            return str(uuid4())
        return None

    def random_date_from_past_50_years(self) -> datetime:
        """Generate a random date in the past 50 years."""

        def random_date(start_date, end_date):
            delta = end_date - start_date
            return start_date + timedelta(days=random.randint(0, delta.days))

        today = datetime.now()
        past_date = today - timedelta(days=50 * 365)  # ~100 years ago
        return random_date(past_date, today)

    def _get_compound_count(self, size: DatasetSize) -> int:
        """Get total number of compounds to generate based on size label."""

        return {
            DatasetSize.XSMALL: 5,
            DatasetSize.SMALL: 1000,
            DatasetSize.LARGE: 1_000_000,
            DatasetSize.XLARGE: 100_000_000,
        }[size]

    def _get_record_count(self, entity_type: EntityType) -> int:
        """Get the number of records to generate for a entity_type, using weighted probabilities."""

        if entity_type not in self.distributions:
            raise ValueError(f"No distribution defined for entity_type: {entity_type.value}")

        dist = self.distributions[entity_type.value]
        return random.choices(dist["choices"], weights=dist["weights"], k=1)[0]

    def _create_directory(self, path: Path) -> None:
        """Create a directory if it does not exist."""

        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)

    def _create_value_types_for_entity_types(self) -> None:
        """Create value types and entity_types for the schema builder."""

        value_types = {}
        shared_list = [ValueType.string, ValueType.int, ValueType.double, ValueType.datetime, ValueType.uuid]
        for entity_type in (EntityType.COMPOUND, EntityType.BATCH, EntityType.ASSAY_RUN, EntityType.ASSAY):
            value_types[entity_type] = shared_list
        value_types[EntityType.ASSAY_RESULT] = [ValueType.string, ValueType.int, ValueType.double, ValueType.bool]
        self.value_types = value_types

    def generate_unique_date(self, unique_dates: set):
        assay_run_date = self._generate_value("datetime")
        while assay_run_date in unique_dates:
            assay_run_date = self._generate_value("datetime")
            unique_dates.add(assay_run_date)
        return assay_run_date

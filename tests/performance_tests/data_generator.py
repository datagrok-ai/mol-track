import random
import json
from typing import Dict, Any, Tuple


class DataGenerator:
    value_types = ["string", "integer", "float"]
    scopes = {  # Mapping schema file names to corresponding scopes
        "batch": {
            "num_properties": 10,
            "property_names": [],
            "file_name": "batches_schema.json",
            "details": "batch_details",
        },
        "compound": {
            "num_properties": 5,
            "property_names": [],
            "file_name": "compounds_schema.json",
            "details": "compound_details",
        },
        # "assays_schema.json": { "name": "compound", "num_properties": 5},
        # "assay_runs_schema.json": { "name": "compound", "num_properties": 5},
        "assay_result": {
            "num_properties": 5,
            "property_names": [],
            "file_name": "assay_results_schema.json",
            "details": "assay_results_details",
        },
    }

    def __init__(self):
        pass

    def get_number_of_records(self):
        # TODO: This will be tailored based on the type of record
        choices = [1, 2, 3, 4, 5]  # Possible batch counts
        probabilities = [0.80, 0.05, 0.05, 0.05, 0.05]  # Your given probabilities

        return random.choices(choices, probabilities, k=1)[0]

    def generate_schema_and_mapping_files(self, output_dir="."):
        """Generate schema files and save to the output directory."""
        for scope, data in self.scopes.items():
            self._create_property_names(scope)
            dicts = self._generate_schema_and_mapping(scope)
            # Save schema as a JSON file
            file_paths = [f"{output_dir}/{scope}_schema.json", f"{output_dir}/{scope}_mapping.json"]
            for i, file_path in enumerate(file_paths):
                with open(file_path, "w") as f:
                    json.dump(dicts[i], f, indent=4)
                print(f"Generated schema file: {file_path}")

    def _generate_schema_and_mapping(self, scope) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """Generate a full schema for a given file and scope."""
        properties = [self._generate_property(scope) for _ in range(self.scopes[scope]["num_properties"])]
        mapping = {p["name"]: f"{self.scopes[scope]['details']}.{p['name']}" for p in properties}

        return {"properties": properties}, mapping

    def _create_property_names(self, scope):
        lst = self.scopes[scope]["property_names"]

        for i in range(self.scopes[scope]["num_properties"]):
            lst.append(f"{scope}_detail_name_{i + 1}")

    def _generate_property(self, scope):
        """Generate a single property based on the given scope."""
        value_type = random.choice(self.value_types)  # Randomly select a value type
        unique_name = self.scopes[scope]["property_names"].pop()  # Unique name for the property
        unit = ""  # Default to empty string for unit

        return {
            "name": unique_name,
            "scope": scope,
            "property_class": "measured",
            "value_type": value_type,
            "unit": unit,
        }

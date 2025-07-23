import random
from typing import Dict, Any
from file_writer import FileWriter
from pathlib import Path
from schema_rules import SCHEMA_GENERATING_RULES


class SchemaBuilder:
    def __init__(self, output_dir, value_types):
        self.output_dir = Path(output_dir)
        self.schema_generating_rules = SCHEMA_GENERATING_RULES
        self.value_types = value_types
        self.schemas: Dict[str, Dict[str, Any]] = {}

    def generate_schema_and_mapping(self) -> None:
        """Generate schema and mapping files for every scope."""
        writer = FileWriter(self.output_dir)
        for scope, rule in self.schema_generating_rules.items():
            schema = {"properties": [], "identity": [], "identity_cols": []}
            mapping = {}

            self._generate_identity_properties(scope, rule, schema, mapping)
            self._generate_random_properties(scope, rule, schema, mapping)

            writer.write_json(f"schemas/{scope}_schema.json", {"properties": schema["properties"]})
            writer.write_json(f"{scope}_mapping.json", mapping)

            self.schemas[scope] = schema

        return self.schemas

    def _generate_identity_properties(
        self, scope: str, rule: Dict[str, Any], schema: Dict[str, Any], mapping: Dict[str, str]
    ) -> None:
        """Generate identity properties and update schema and mapping dictionaries."""
        identity = rule["identity_columns"]
        for i, prop_name in enumerate(identity["property_names"]):
            prop_type = identity["property_type"][i]
            prop_mapping = identity["property_mapping"][i]

            prop = self._generate_property(scope, prop_name, prop_type)

            if prop_mapping == rule["alias"]:
                schema["properties"].append(prop)
            else:
                schema["identity"].append(prop)

            schema["identity_cols"].append(prop["name"])
            mapping[prop_name] = f"{prop_mapping}.{prop_name}" if prop_mapping else prop_name

    def _generate_random_properties(
        self, scope: str, rule: Dict[str, Any], schema: Dict[str, Any], mapping: Dict[str, str]
    ) -> None:
        """Generate random properties and update schema and mapping dictionaries."""
        for i in range(rule["num_properties"]):
            prop_name = f"{scope}_detail_name_{i + 1}"
            prop = self._generate_property(scope, prop_name)
            schema["properties"].append(prop)
            mapping[prop_name] = f"{rule['alias']}.{prop_name}"

    def _generate_property(self, scope: str, name: str, value_type: str = None) -> Dict[str, Any]:
        """Generate a single property dictionary."""
        return {
            "name": name,
            "scope": scope,
            "property_class": "measured",
            "value_type": value_type or random.choice(self.value_types),
            "unit": "",
        }

import random
from typing import Dict, Any
from tests.performance_tests.file_writer import FileWriter
from pathlib import Path
from tests.performance_tests.schema_rules import SCHEMA_GENERATING_RULES
from app.utils.enums import PropertyClass, ScopeClass, ValueType


class SchemaBuilder:
    def __init__(self, output_dir: Path, config: Dict[str, Any], value_types: Dict[ScopeClass, list[ValueType]]):
        self.output_dir = output_dir
        self.schema_generating_rules = SCHEMA_GENERATING_RULES
        self.value_types = value_types
        self.schemas: Dict[str, Dict[str, Any]] = {}
        self.details_columns = config.get("details_columns", {})

    def generate_schema_and_mapping(self) -> Dict[str, Dict[str, Any]]:
        """Generate schema and mapping files for every scope."""

        writer = FileWriter(self.output_dir)
        for scope, rule in self.schema_generating_rules.items():
            schema = {"properties": [], "identity": [], "identity_cols": []}
            mapping = {}

            self._generate_identity_properties(scope, rule, schema, mapping)
            self._generate_random_properties(scope, rule, schema, mapping)

            writer.write_json(f"schemas/{scope.value}_schema.json", {"properties": schema["properties"]})
            writer.write_json(f"{scope.value}_mapping.json", mapping)

            self.schemas[scope] = schema

        return self.schemas

    def _generate_identity_properties(
        self, scope: ScopeClass, rule: Dict[str, Any], schema: Dict[str, Any], mapping: Dict[str, str]
    ) -> None:
        """Generate identity properties and update schema and mapping dictionaries."""

        identity = rule.get("identity_columns", {})
        property_names = identity.get("property_names", [])
        property_types = identity.get("property_type", [])
        property_mappings = identity.get("property_mapping", [])
        alias = rule.get("alias", "")

        for i, prop_name in enumerate(property_names):
            prop_type = property_types[i]
            prop_mapping = property_mappings[i]

            prop = self._generate_property(scope, prop_name, prop_type)

            if prop_mapping == alias:
                schema["properties"].append(prop)
            else:
                schema["identity"].append(prop)

            schema["identity_cols"].append(prop["name"])
            mapping[prop_name] = f"{prop_mapping}.{prop_name}" if prop_mapping else prop_name

    def _generate_random_properties(
        self, scope: ScopeClass, rule: Dict[str, Any], schema: Dict[str, Any], mapping: Dict[str, str]
    ) -> None:
        """Generate random properties and update schema and mapping dictionaries."""

        num_properties = self.details_columns.get(scope.value, 1)
        alias = rule["alias"]
        for i in range(num_properties):
            prop_name = f"{scope.value}_detail_name_{i + 1}"
            prop = self._generate_property(scope, prop_name)
            schema["properties"].append(prop)
            mapping[prop_name] = f"{alias}.{prop_name}"

    def _generate_property(self, scope: ScopeClass, name: str, value_type: ValueType = None) -> Dict[str, Any]:
        """Generate a single property dictionary."""

        return {
            "name": name,
            "scope": scope,
            "property_class": PropertyClass.MEASURED.value,
            "value_type": value_type or random.choice(self.value_types[scope]),
            "unit": "",
        }

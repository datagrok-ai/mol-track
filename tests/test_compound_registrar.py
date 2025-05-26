import json
import pytest

import enums
from tests.conftest import load_csv_file

TEST_SCHEMA_DATA = {
    "properties": [
        {"name": "MolLogP", "scope": "COMPOUND", "property_class": "CALCULATED", "value_type": "double", "unit": ""}
    ],
    "synonym_types": [
        {"name": "batch_corporate_id", "synonym_level": "BATCH", "pattern": ""},
        {"name": "corporate_id", "synonym_level": "COMPOUND", "pattern": ""},
        {"name": "CAS", "synonym_level": "COMPOUND", "pattern": r"\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b"},
        {"name": "common name", "synonym_level": "COMPOUND", "pattern": ""},
        {"name": "USAN", "synonym_level": "COMPOUND", "pattern": ""},
    ],
}

DEFAULT_MAPPING = {
    "smiles": "smiles",
    "common_name": "compounds_synonyms.common_name",
    "cas": "compounds_synonyms.cas",
    "usan": "compounds_synonyms.usan",
    "MolLogP": "properties.MolLogP",
}


@pytest.fixture
def preload_schema(client):
    return client.post("/schema/", json=TEST_SCHEMA_DATA)


def post_compound_data(client, file_path, error_handling):
    csv_data = load_csv_file(file_path)
    files = {"csv_file": (file_path, csv_data, "text/csv")}
    data = {
        "error_handling": error_handling.value,
        "mapping": json.dumps(DEFAULT_MAPPING),
    }
    return client.post("/v1/compounds/", files=files, data=data)


def test_preload_schema(client):
    response = client.post("/schema/", json=TEST_SCHEMA_DATA)
    assert response.status_code == 200

    data = response.json()
    assert data["status"] == "success"

    expected_synonyms = ["batch_corporate_id", "corporate_id", "common name", "CAS", "USAN"]
    expected_properties = ["MolLogP"]
    actual_synonyms = [s["name"] for s in data["created"]["synonym_types"]]
    actual_properties = [p["name"] for p in data["created"]["property_types"]]

    assert "created" in data
    assert sorted(actual_synonyms) == sorted(expected_synonyms)
    assert actual_properties == expected_properties


def test_register_compounds_reject_all(client, preload_schema):
    response = post_compound_data(client, "demo-data/black/2_batches.csv", enums.ErrorHandlingOptions.reject_all)
    assert response.status_code == 400

    result = response.json()["detail"]
    assert "errors" in result
    assert result["failed_rows"] == 1
    assert result["successful_rows"] == 0


def test_register_compounds_reject_row(client, preload_schema):
    response = post_compound_data(client, "demo-data/black/2_batches.csv", enums.ErrorHandlingOptions.reject_row)
    assert response.status_code == 200

    result = response.json()
    assert result["failed_rows"] == 1
    assert result["successful_rows"] == 2

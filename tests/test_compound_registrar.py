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
    return client.post("v1/schema/", json=TEST_SCHEMA_DATA)


def post_compound_data(client, file_path, error_handling):
    csv_data = load_csv_file(file_path)
    files = {"csv_file": (file_path, csv_data, "text/csv")}
    data = {
        "error_handling": error_handling.value,
        "mapping": json.dumps(DEFAULT_MAPPING),
    }
    return client.post("/v1/compounds/", files=files, data=data)


def test_preload_schema(client):
    response = client.post("v1/schema/", json=TEST_SCHEMA_DATA)
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
    assert result["status"] == "Success"

    data = result["data"]
    assert len(data) == 3

    item0 = data[0]
    assert item0["batch_corporate_id"] == "EPA-001-002"
    assert item0["registration_status"] == "failed"
    assert item0["registration_error_message"] == "400: Invalid SMILES string"

    item1 = data[1]
    assert item1["batch_corporate_id"] == "EPA-001-001"
    assert item1["registration_status"] == "not_processed"
    assert item1["registration_error_message"] is None

    item2 = data[2]
    assert item2["batch_corporate_id"] == "EPA-120-001"
    assert item2["registration_status"] == "not_processed"
    assert item2["registration_error_message"] is None


def test_register_compounds_reject_row(client, preload_schema):
    response = post_compound_data(client, "demo-data/black/2_batches.csv", enums.ErrorHandlingOptions.reject_row)
    assert response.status_code == 200

    result = response.json()
    data = result["data"]
    assert isinstance(data, list)
    assert len(data) == 3

    failed = [item for item in data if item["registration_status"] == "failed"]
    successful = [item for item in data if item["registration_status"] == "successful"]

    assert len(failed) == 3
    assert len(successful) == 0

    assert data[0]["batch_corporate_id"] == "EPA-001-002"
    assert data[0]["registration_status"] == "failed"
    assert data[0]["registration_error_message"] == "400: Invalid SMILES string"

    assert data[1]["batch_corporate_id"] == "EPA-001-001"
    assert data[1]["registration_status"] == "failed"
    assert data[1]["registration_error_message"] == "400: Unknown synonym type: common_name"

    assert data[2]["batch_corporate_id"] == "EPA-120-001"
    assert data[2]["registration_status"] == "failed"
    assert data[2]["registration_error_message"] == "400: Unknown synonym type: common_name"

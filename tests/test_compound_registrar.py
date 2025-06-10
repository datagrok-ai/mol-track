import json
from tests.conftest import TEST_SCHEMA_DATA

import enums
from tests.conftest import load_csv_file

DEFAULT_MAPPING = {
    "smiles": "smiles",
    "common_name": "compounds_synonyms.common_name",
    "cas": "compounds_synonyms.cas",
    "usan": "compounds_synonyms.usan",
    "MolLogP": "properties.MolLogP",
}


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

    expected_synonyms = ["batch_corporate_id", "corporate_id", "common_name", "cas", "usan"]
    expected_properties = ["MolLogP", "acquired_date"]
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
    assert item0["registration_status"] == "failed"
    assert item0["registration_error_message"] == "400: Invalid SMILES string"

    item1 = data[1]
    assert item1["registration_status"] == "not_processed"
    assert item1["registration_error_message"] is None

    item2 = data[2]
    assert item2["registration_status"] == "not_processed"
    assert item2["registration_error_message"] is None


def test_register_compounds_reject_row(client, preload_schema):
    response = post_compound_data(client, "demo-data/black/2_batches.csv", enums.ErrorHandlingOptions.reject_row)
    assert response.status_code == 200

    result = response.json()
    data = result["data"]
    assert isinstance(data, list)
    assert len(data) == 3

    item0 = data[0]
    assert item0["registration_status"] == "failed"
    assert item0["registration_error_message"] == "400: Invalid SMILES string"

    item1 = data[1]
    assert item1["registration_status"] == "success"
    assert item1["registration_error_message"] is None

    item2 = data[2]
    assert item2["registration_status"] == "success"
    assert item2["registration_error_message"] is None

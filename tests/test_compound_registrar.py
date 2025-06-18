import pytest
import enums
from tests.conftest import BLACK_DIR, read_json, _preload_compounds


def extract_name_scope(items):
    return sorted([{"name": item["name"], "scope": item["scope"].upper()} for item in items], key=lambda x: x["name"])


def assert_name_scope_equal(actual, expected):
    actual_name_scope = extract_name_scope(actual)
    expected_name_scope = extract_name_scope(expected)
    assert actual_name_scope == expected_name_scope


@pytest.mark.parametrize(
    "endpoint,schema_file,response_key,expected_keys",
    [
        ("/v1/schema/compounds", "compounds_schema.json", None, ["properties", "synonym_types"]),
        ("/v1/schema/compounds/synonyms", "compounds_schema.json", None, ["synonym_types"]),
        ("/v1/schema/batches", "batches_schema.json", "properties", ["properties", "synonym_types"]),
        ("/v1/schema/batches/synonyms", "batches_schema.json", "synonym_types", ["synonym_types"]),
    ],
)
def test_schema(client, endpoint, schema_file, response_key, expected_keys, preload_schema):
    response = client.get(endpoint)
    assert response.status_code == 200

    response_data = response.json()
    actual = response_data if response_key is None else response_data[response_key]
    expected_data = read_json(BLACK_DIR / schema_file)
    expected = []
    for key in expected_keys:
        expected.extend(expected_data.get(key, []))

    assert_name_scope_equal(actual, expected)

    if "/v1/schema/batches" in endpoint:
        assert "additions" in response_data


def test_register_compounds_reject_all(client, preload_schema):
    response = _preload_compounds(
        client, BLACK_DIR / "compounds.csv", BLACK_DIR / "compounds_mapping.json", enums.ErrorHandlingOptions.reject_all
    )
    assert response.status_code == 400

    result = response.json()["detail"]
    assert result["status"] == "Success"

    data = result["data"]
    assert len(data) == 54

    item8 = data[8]
    assert item8["registration_status"] == "failed"
    assert item8["registration_error_message"] == "400: Invalid SMILES string"

    for item in data[9:]:
        assert item["registration_status"] == "not_processed"
        assert item["registration_error_message"] is None


def test_register_compounds_reject_row(client, preload_schema):
    response = _preload_compounds(client, BLACK_DIR / "compounds.csv", BLACK_DIR / "compounds_mapping.json")
    assert response.status_code == 200

    result = response.json()
    data = result["data"]
    assert isinstance(data, list)
    assert len(data) == 54

    item8 = data[8]
    assert item8["registration_status"] == "failed"
    assert item8["registration_error_message"] == "400: Invalid SMILES string"

    item9 = data[9]
    assert item9["registration_status"] == "success"
    assert item9["registration_error_message"] is None

import pytest
import app.utils.enums as enums
from tests.conftest import BLACK_DIR, _preload_batches


def test_register_batches_without_mapping(client, preload_schema):
    register_response = _preload_batches(client, BLACK_DIR / "batches.csv")
    assert register_response.status_code == 200

    register_data = register_response.json()
    assert len(register_data) == 54

    get_response = client.get("/v1/batches/")
    assert get_response.status_code == 200

    batches = get_response.json()
    expected_properties = {
        "Purity",
        "Synthesized Date",
        "corporate_batch_id",
        "Responsible Party",
        "ELN Reference",
        "Project",
        "EPA Batch ID",
        "Source Batch Code",
        "Source",
    }

    properties = batches[0].get("properties", [])
    assert len(properties) == len(expected_properties), (
        f"Expected {len(expected_properties)} properties, got {len(properties)}"
    )

    prop_names = {p["name"] for p in properties}
    assert prop_names == expected_properties, f"Property names mismatch: {prop_names} != {expected_properties}"


@pytest.mark.skip(reason="No test datasets contain invalid records to validate 'reject all' behaviour.")
def test_register_batches_reject_all(client, preload_schema, preload_additions):
    response = _preload_batches(
        client, BLACK_DIR / "batches.csv", BLACK_DIR / "batches_mapping.json", enums.ErrorHandlingOptions.reject_all
    )
    assert response.status_code == 400

    result = response.json()["detail"]
    assert result["status"] == "Success"

    data = result["data"]
    assert len(data) == 54

    assert data[8]["registration_status"] == "failed"
    assert data[8]["registration_error_message"] == "400: Invalid SMILES string"

    for item in data[9:]:
        assert item["registration_status"] == "not_processed"
        assert item["registration_error_message"] is None


@pytest.mark.skip(reason="No test datasets contain invalid records to validate 'reject row' behaviour.")
def test_register_batches_reject_row(client, preload_schema, preload_additions):
    response = _preload_batches(client, BLACK_DIR / "batches.csv", BLACK_DIR / "batches_mapping.json")
    assert response.status_code == 200

    result = response.json()
    data = result["data"]
    assert isinstance(data, list)
    assert len(data) == 54

    assert data[8]["registration_status"] == "failed"
    assert data[8]["registration_error_message"] == "400: Invalid SMILES string"

    assert data[9]["registration_status"] == "success"
    assert data[9]["registration_error_message"] is None

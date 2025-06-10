import json
import enums
from tests.conftest import load_csv_file

DEFAULT_MAPPING = {
    "smiles": "smiles",
    "batch_corporate_id": "batches_synonyms.batch_corporate_id",
    "acquistion_date": "batches_properties.acquired_date",
    "HCl": "batches_additions.HCl",
    "common_name": "compounds_synonyms.common_name",
    "cas": "compounds_synonyms.cas",
    "usan": "compounds_synonyms.usan",
    "MolLogP": "properties.MolLogP",
}


def post_batch_data(client, file_path, error_handling):
    csv_data = load_csv_file(file_path)
    files = {"csv_file": (file_path, csv_data, "text/csv")}
    data = {
        "error_handling": error_handling.value,
        "mapping": json.dumps(DEFAULT_MAPPING),
    }
    return client.post("/v1/batches/", files=files, data=data)


def test_register_batches_reject_all(client, preload_schema):
    response = post_batch_data(
        client, "demo-data/black/2_batches_with_additions.csv", enums.ErrorHandlingOptions.reject_all
    )
    assert response.status_code == 400

    result = response.json()["detail"]
    assert result["status"] == "Success"

    data = result["data"]
    assert len(data) == 2

    item0 = data[0]
    assert item0["registration_status"] == "failed"
    assert item0["registration_error_message"] == "400: Unknown addition: HCl"

    item1 = data[1]
    assert item1["registration_status"] == "not_processed"
    assert item1["registration_error_message"] is None


def test_register_batches_reject_row(client, preload_schema, uploaded_additions):
    response = post_batch_data(
        client, "demo-data/black/2_batches_with_additions.csv", enums.ErrorHandlingOptions.reject_all
    )
    assert response.status_code == 200

    result = response.json()
    data = result["data"]
    assert isinstance(data, list)
    assert len(data) == 2

    item0 = data[0]
    assert item0["registration_status"] == "success"
    assert item0["registration_error_message"] is None

    item1 = data[1]
    assert item1["registration_status"] == "success"
    assert item1["registration_error_message"] is None

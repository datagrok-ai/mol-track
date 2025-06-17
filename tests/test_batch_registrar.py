import enums
from tests.conftest import BLACK_DIR, _preload_batches


def test_register_batches_reject_all(client, preload_schema, preload_additions, preload_batches):
    response = _preload_batches(
        client, BLACK_DIR / "batches.csv", BLACK_DIR / "batches_mapping.json", enums.ErrorHandlingOptions.reject_all
    )
    assert response.status_code == 400

    result = response.json()["detail"]
    assert result["status"] == "Success"

    data = result["data"]
    assert len(data) == 54

    item0 = data[0]
    assert item0["registration_status"] == "failed"
    assert (
        item0["registration_error_message"] == "400: Compound with InChIKey QQVIHTHCMHWDBS-UHFFFAOYSA-N already exists"
    )

    for item in data[1:]:
        assert item["registration_status"] == "not_processed"
        assert item["registration_error_message"] is None


def test_register_batches_reject_row(client, preload_schema, preload_additions, preload_batches):
    response = _preload_batches(client, BLACK_DIR / "batches.csv", BLACK_DIR / "batches_mapping.json")
    assert response.status_code == 200

    result = response.json()
    data = result["data"]
    assert isinstance(data, list)
    assert len(data) == 54

    item0 = data[0]
    assert item0["registration_status"] == "failed"
    assert (
        item0["registration_error_message"] == "400: Compound with InChIKey QQVIHTHCMHWDBS-UHFFFAOYSA-N already exists"
    )

    for item in data[1:]:
        assert item["registration_status"] == "failed"
        assert item["registration_error_message"] is not None

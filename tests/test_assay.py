from tests.conftest import _preload_assays, _preload_assay_runs, BLACK_DIR


def test_create_assay(client, preload_schema):
    response = _preload_assays(client, BLACK_DIR / "assays.json")
    assert response.status_code == 200

    assays = response.json()["created"]
    assert isinstance(assays, list)
    assert len(assays) == 1
    assert assays[0]["name"] == "Hepatocyte Stability"


def test_get_all_assays(client, preload_schema, preload_assays):
    response = client.get("/v1/assays/")
    assert response.status_code == 200

    assays = response.json()
    assert isinstance(assays, list)
    assert len(assays) == 1

    assay = assays[0]
    assert assay["name"] == "Hepatocyte Stability"
    assert "properties" in assay and isinstance(assay["properties"], list) and len(assay["properties"]) > 0
    assert "assay_details" in assay and isinstance(assay["assay_details"], list) and len(assay["assay_details"]) > 0
    assert (
        "property_requirements" in assay
        and isinstance(assay["property_requirements"], list)
        and len(assay["property_requirements"]) > 0
    )


def test_create_assay_run(client, preload_schema, preload_assays):
    response = _preload_assay_runs(client, BLACK_DIR / "assay_runs.csv", BLACK_DIR / "assay_runs_mapping.json")
    assert response.status_code == 200
    result = response.json()

    assert isinstance(result, list)
    assert len(result) == 9
    assert all(run.get("registration_status") == "success" for run in result)


def test_get_assay_run_by_id(client, preload_schema, preload_assays, preload_assay_runs):
    response = client.get("/v1/assay_runs/2")
    assert response.status_code == 200

    data = response.json()
    assert data["id"] == 2
    assert data["name"] == "Hepatocyte Stability2024-02-02"
    assert data["assay_id"] == 1

    assay = data.get("assay")
    assert assay and assay["id"] == 1
    assert assay["name"] == "Hepatocyte Stability"

    details = data.get("assay_run_details")
    assert isinstance(details, list)
    assert len(details) == 3

    detail_22 = next((d for d in details if d["property_id"] == 22), None)
    assert detail_22 is not None
    assert detail_22["value_string"] == "Human"

    props = data.get("properties")
    assert isinstance(props, list)
    assert any(p["name"] == "Cell Species" for p in props)


def test_get_assay_runs(client, preload_schema, preload_assays, preload_assay_runs):
    response = client.get("/v1/assay_runs/")
    assert response.status_code == 200

    data = response.json()
    assert len(data) == 9

    first = data[0]
    assert first["name"] == "Hepatocyte Stability2024-02-01"
    assert first["assay_id"] == 1
    assert first["id"] == 1

    assay = first.get("assay")
    assert assay is not None
    assert assay.get("name") == "Hepatocyte Stability"

    details = first.get("assay_run_details", [])
    assert isinstance(details, list)
    assert len(details) > 0

    properties = first.get("properties", [])
    assert isinstance(properties, list)
    assert len(properties) > 0

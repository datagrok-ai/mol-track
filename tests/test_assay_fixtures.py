def test_get_all_assays(client, preload_schema, preload_assays):
    response = client.get("/v1/assays/")
    assert response.status_code == 200

    assays = response.json()
    assert isinstance(assays, list)
    assert len(assays) >= 1

    assay = assays[0]
    assert assay["name"] == "Hepatocyte Stability"
    assert "properties" in assay and isinstance(assay["properties"], list) and len(assay["properties"]) > 0
    assert "assay_details" in assay and isinstance(assay["assay_details"], list) and len(assay["assay_details"]) > 0
    assert (
        "property_requirements" in assay
        and isinstance(assay["property_requirements"], list)
        and len(assay["property_requirements"]) > 0
    )


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

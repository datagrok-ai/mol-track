import pytest
from app import models
from tests.conftest import _preload_assay_results, _preload_assays, _preload_assay_runs, BLACK_DIR
from tests.utils.test_base_registrar import BaseRegistrarTest

# === Tests for Assay endpoints ===


def test_create_assay(client, preload_schema, api_headers):
    response = _preload_assays(client, BLACK_DIR / "assays.json", api_headers)
    assert response.status_code == 200

    assays = response.json()["created"]
    assert isinstance(assays, list)
    assert len(assays) == 1
    assert assays[0]["name"] == "Hepatocyte Stability"


def test_get_all_assays(client, preload_schema, preload_assays, api_headers):
    response = client.get("/v1/assays/", headers=api_headers)
    assert response.status_code == 200

    assays = response.json()
    assert isinstance(assays, list)
    assert len(assays) == 1

    assay = assays[0]
    assert assay["name"] == "Hepatocyte Stability"
    assert "properties" in assay and isinstance(assay["properties"], list) and len(assay["properties"]) > 0
    assert (
        "property_requirements" in assay
        and isinstance(assay["property_requirements"], list)
        and len(assay["property_requirements"]) > 0
    )


def test_get_assay_by_id(client, preload_schema, preload_assays, api_headers):
    all_response = client.get("/v1/assays/", headers=api_headers)
    assert all_response.status_code == 200

    assays = all_response.json()
    assert isinstance(assays, list) and len(assays) > 0
    expected_assay = assays[0]
    assay_id = expected_assay["id"]

    response = client.get(f"/v1/assays/{assay_id}", headers=api_headers)
    assert response.status_code == 200, f"Failed to fetch assay {assay_id}"
    assay = response.json()

    assert assay["id"] == expected_assay["id"]
    assert assay["name"] == expected_assay["name"]
    assert assay["properties"] == expected_assay["properties"]
    assert assay["property_requirements"] == expected_assay["property_requirements"]

    missing = client.get("/v1/assays/999999", headers=api_headers)
    assert missing.status_code == 404, "Expected 404 for non-existent assay"


# === Tests for Assay Run endpoints ===


@pytest.mark.usefixtures("preload_assays")
class TestAssayRunRegistrar(BaseRegistrarTest):
    entity_name = "assay_runs"
    expected_properties = {"Cell Species", "Cell Lot", "Assay Run Date"}
    expected_first = {
        "cell species": "Human",
        "cell lot": " H1",
        "assay date": "2024-02-01",
    }
    preload_func = staticmethod(_preload_assay_runs)
    get_response_model = models.AssayRunResponse
    default_entity_count = 9

    def test_get_assay_run_by_id(self, client, preload_schema, preload_assay_runs, api_headers):
        all_response = client.get("/v1/assay_runs/", headers=api_headers)
        assert all_response.status_code == 200

        assay_runs = all_response.json()
        assert isinstance(assay_runs, list) and len(assay_runs) > 0
        expected_run = assay_runs[0]
        run_id = expected_run["id"]

        response = client.get(f"/v1/assay_runs/{run_id}", headers=api_headers)
        assert response.status_code == 200, f"Failed to fetch assay run {run_id}"

        assay_run = response.json()
        assert assay_run["id"] == expected_run["id"]
        assert assay_run["name"] == expected_run["name"]
        assert assay_run["properties"] == expected_run["properties"]
        assert assay_run["assay_run_details"] == expected_run["assay_run_details"]

        missing = client.get("/v1/assay_runs/999999", headers=api_headers)
        assert missing.status_code == 404, "Expected 404 for non-existent assay run"

    @pytest.mark.skip(reason="Registering assay run requires a mapping to determine assay names.")
    def test_register_without_mapping(self, **kwargs):
        pass


# === Tests for Assay Result endpoints ===


@pytest.mark.usefixtures("preload_compounds", "preload_batches", "preload_assays", "preload_assay_runs")
class TestAssayResultRegistrar(BaseRegistrarTest):
    entity_name = "assay_results"
    expected_properties = {}
    expected_first = {
        "Concentration (µM)": "0.1",
        "Mean HTC recovery (%; n=2)": "59.51",
        "SD HTC recovery (%; n=2)": "84.16",
    }
    preload_func = staticmethod(_preload_assay_results)
    get_response_model = models.AssayResultResponse
    default_entity_count = 238

    @pytest.mark.skip(reason="Registering assay result requires a mapping to determine assay and batch.")
    def test_register_without_mapping(self, **kwargs):
        pass

import pytest
from tests.conftest import _preload_batches
from tests.utils.test_base_registrar import BaseRegistrarTest
from tests.utils.test_mixins import CRUDTestsMixin, SynonymTestsMixin
from app import models


@pytest.fixture
def first_batch_with_synonyms(client, preload_schema, preload_batches):
    response = client.get("/v1/batches/")
    assert response.status_code == 200
    batches = response.json()
    first = batches[0]
    synonyms = [p for p in first["properties"] if p["semantic_type_id"] == 1]
    return first, synonyms


class TestBatchesRegistrar(BaseRegistrarTest, SynonymTestsMixin, CRUDTestsMixin):
    entity_name = "batches"
    allow_put = False
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
    expected_first = {
        "EPA Batch ID": "EPA-001-001",
        "Source": "EPA",
        "Source Batch Code": "EPA-001",
        "Responsible Party": "chemist 1",
        "ELN Reference": "cchemist_123",
        "Purity": "93",
        "Project": "Project 1",
    }
    preload_func = staticmethod(_preload_batches)
    get_response_model = models.BatchResponse
    first_entity_fixture_name = "first_batch_with_synonyms"

    @pytest.mark.skip(reason="PUT endpoint not implemented for batches")
    def test_put_entity(self, *args, **kwargs):
        pass

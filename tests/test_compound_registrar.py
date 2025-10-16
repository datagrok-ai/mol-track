import pytest
from app.utils import enums
from app import models
from tests.conftest import BLACK_DIR, read_json, _preload_compounds
from tests.utils.test_base_registrar import BaseRegistrarTest


def extract_name_entity_type(items):
    return sorted(
        [{"name": item["name"], "entity_type": item["entity_type"].upper()} for item in items], key=lambda x: x["name"]
    )


def assert_name_entity_type_equal(actual, expected):
    actual_name_entity_type = extract_name_entity_type(actual)
    expected_name_entity_type = extract_name_entity_type(expected)
    assert actual_name_entity_type == expected_name_entity_type


@pytest.fixture
def first_compound_with_synonyms(client, preload_schema, preload_compounds):
    response = client.get("/v1/compounds/")
    assert response.status_code == 200
    compounds = response.json()
    assert compounds, "No compounds returned"

    first = compounds[0]
    synonyms = [p for p in first["properties"] if p["semantic_type_id"] == 1]
    assert synonyms, f"No synonym properties found for compound {first['canonical_smiles']}"

    return first, synonyms


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

    if "/compounds" in endpoint:
        expected.append({"name": "corporate_compound_id", "entity_type": enums.EntityType.COMPOUND})
    if "/batches" in endpoint:
        expected.append({"name": "corporate_batch_id", "entity_type": enums.EntityType.BATCH})

    assert_name_entity_type_equal(actual, expected)

    if "/v1/schema/batches" in endpoint:
        assert "additions" in response_data


class TestCompoundsRegistrar(BaseRegistrarTest):
    entity_name = "compounds"
    expected_properties = {
        "EPA Compound ID",
        "corporate_compound_id",
        "MolLogP",
        "Source Compound Code",
        "CAS",
        "Source",
        "Common Name",
    }
    expected_first = {
        "Common Name": "1,3-Benzenedicarboxylic acid",
        "CAS": "121-91-5",
        "EPA Compound ID": "EPA-001",
        "MolLogP": 1.082999945,
        "Source": "EPA",
        "Source Compound Code": "EPA-001",
    }
    preload_func = staticmethod(_preload_compounds)
    get_response_model = models.CompoundResponse
    first_entity_fixture_name = "first_compound_with_synonyms"

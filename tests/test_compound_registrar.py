import pytest
from app.utils import enums
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

    if "/compounds" in endpoint:
        expected.append({"name": "corporate_compound_id", "scope": enums.ScopeClass.COMPOUND})
    if "/batches" in endpoint:
        expected.append({"name": "corporate_batch_id", "scope": enums.ScopeClass.BATCH})

    assert_name_scope_equal(actual, expected)

    if "/v1/schema/batches" in endpoint:
        assert "additions" in response_data


def test_register_compounds_without_mapping(client, preload_schema):
    register_response = _preload_compounds(client, BLACK_DIR / "compounds.csv")
    assert register_response.status_code == 200

    register_data = register_response.json().get("data", [])
    assert len(register_data) == 54

    get_response = client.get("/v1/compounds/")
    assert get_response.status_code == 200
    compounds = get_response.json()

    def assert_properties(compound, expected_props, index):
        properties = compound.get("properties", [])
        assert len(properties) == len(expected_props), (
            f"[Compound {index}] Expected {len(expected_props)} properties, got {len(properties)}"
        )
        names = {p["name"] for p in properties}
        assert names == expected_props, f"[Compound {index}] Property names mismatch: {names} != {expected_props}"

    assert_properties(compounds[0], {"epa_compound_id", "corporate_compound_id", "MolLogP"}, index=0)
    assert_properties(compounds[8], {"epa_compound_id", "corporate_compound_id"}, index=8)


@pytest.mark.skip(reason="No test datasets contain invalid records to validate 'reject all' behaviour.")
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


@pytest.mark.skip(reason="No test datasets contain invalid records to validate 'reject row' behaviour.")
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


def test_get_compounds_list(client, preload_schema, preload_compounds):
    response = client.get("/v1/compounds/")
    assert response.status_code == 200

    result = response.json()
    assert isinstance(result, list)
    assert len(result) == 54

    first = result[0]
    assert first["id"] == 1
    assert first["canonical_smiles"] == "O=C(O)c1cccc(C(=O)O)c1"
    assert first["inchikey"] == "QQVIHTHCMHWDBS-UHFFFAOYSA-N"

    props = {p["name"]: p for p in first["properties"]}
    assert props["epa_compound_id"]["value_string"] == "EPA-001"
    assert props["cas"]["value_string"] == "121-91-5"
    assert props["common_name"]["value_string"].strip() == "1,3-Benzenedicarboxylic acid"
    assert abs(props["MolLogP"]["value_num"] - 1.083) < 1e-3


def test_get_compound_by_id(client, preload_schema, preload_compounds):
    response = client.get("/v1/compounds/2")
    assert response.status_code == 200

    result = response.json()
    assert result["id"] == 2
    assert result["canonical_smiles"] == "Oc1ccc(C(c2ccc(O)cc2)(C(F)(F)F)C(F)(F)F)cc1"
    assert result["inchikey"] == "ZFVMWEVVKGLCIJ-UHFFFAOYSA-N"
    assert result["inchi"] == (
        "InChI=1S/C15H10F6O2/c16-14(17,18)13(15(19,20)21,9-1-5-11(22)6-2-9)10-3-7-12(23)8-4-10/h1-8,22-23H"
    )
    assert result["is_archived"] is False
    assert result["batches"] == []

    props = {p["name"]: p for p in result["properties"]}

    assert props["epa_compound_id"]["value_string"] == "EPA-002"
    assert props["cas"]["value_string"] == "1478-61-1"
    assert props["common_name"]["value_string"].strip() == "Bisphenol AF"
    assert abs(props["MolLogP"]["value_num"] - 4.5085) < 1e-3

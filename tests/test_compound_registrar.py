import pytest
from app.utils import enums
from tests.conftest import BLACK_DIR, preload_entity, read_json, _preload_compounds


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


def test_register_compounds_without_mapping(client, preload_schema):
    register_response = _preload_compounds(client, BLACK_DIR / "compounds.csv")
    assert register_response.status_code == 200

    register_data = register_response.json()
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

    assert_properties(
        compounds[0],
        {"EPA Compound ID", "corporate_compound_id", "MolLogP", "Source Compound Code", "CAS", "Source", "Common Name"},
        index=0,
    )
    assert_properties(
        compounds[8],
        {"EPA Compound ID", "corporate_compound_id", "Source Compound Code", "CAS", "Source", "Common Name"},
        index=8,
    )


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

    assert "batches" in first
    assert isinstance(first["batches"], list)

    assert "properties" in first
    assert isinstance(first["properties"], list)
    props = {p["name"]: p for p in first["properties"]}
    assert props["EPA Compound ID"]["value_string"] == "EPA-001"
    assert props["CAS"]["value_string"] == "121-91-5"
    assert props["Common Name"]["value_string"].strip() == "1,3-Benzenedicarboxylic acid"
    assert abs(props["MolLogP"]["value_num"] - 1.083) < 1e-3


def test_get_compound_by_any_synonym(client, preload_schema, preload_compounds, first_compound_with_synonyms):
    first, synonyms = first_compound_with_synonyms

    for prop in synonyms:
        synonym_name = prop["name"]
        synonym_value = prop["value_string"]

        resp_val = client.get(f"/v1/compounds?property_value={synonym_value}")
        assert resp_val.status_code == 200
        assert resp_val.json()["id"] == first["id"]

        resp_name = client.get(f"/v1/compounds?property_value={synonym_value}&property_name={synonym_name}")
        assert resp_name.status_code == 200
        assert resp_name.json()["id"] == first["id"]


def test_get_compound_synonyms(client, preload_schema, preload_compounds, first_compound_with_synonyms):
    _, synonyms = first_compound_with_synonyms
    for prop in synonyms:
        resp = client.get(f"/v1/compounds/synonyms?property_value={prop['value_string']}")
        assert resp.status_code == 200
        props = resp.json()
        assert all(p["semantic_type_id"] == 1 for p in props), "Non-synonym returned"


def test_get_compound_properties(client, preload_schema, preload_compounds, first_compound_with_synonyms):
    first, synonyms = first_compound_with_synonyms
    for prop in synonyms:
        resp = client.get(f"/v1/compounds/properties?property_value={prop['value_string']}")
        assert resp.status_code == 200
        props = resp.json()
        returned_names = {p["name"] for p in props}
        original_names = {p["name"] for p in first["properties"]}
        assert returned_names == original_names


@pytest.mark.parametrize(
    "update_payload",
    [
        ({"is_archived": True}),
        ({"canonical_smiles": "CCC"}),
    ],
)
def test_put_compound(client, preload_schema, preload_compounds, first_compound_with_synonyms, update_payload):
    first, synonyms = first_compound_with_synonyms
    corporate_compound_prop = next((p for p in first["properties"] if p["name"] == "corporate_compound_id"), None)
    assert corporate_compound_prop is not None, "No Corporate Compound ID found"
    corporate_compound_id = corporate_compound_prop["value_string"]

    response = client.put(f"/v1/compounds/{corporate_compound_id}", json=update_payload)
    assert response.status_code == 200
    data = response.json()

    for key, value in update_payload.items():
        assert data[key] == value, f"Expected {key} to be {value}, but got {data[key]}"


def test_delete_compound(client, preload_schema, preload_compounds, first_compound_with_synonyms):
    first, synonyms = first_compound_with_synonyms
    corporate_compound_prop = next((p for p in first["properties"] if p["name"] == "corporate_compound_id"), None)
    assert corporate_compound_prop is not None, "No Corporate Compound ID found"
    corporate_compound_id = corporate_compound_prop["value_string"]

    response_delete = client.delete(f"/v1/compounds/{corporate_compound_id}")
    assert response_delete.status_code == 200

    response_get = client.get(f"/v1/compounds/properties?property_value={corporate_compound_id}")
    assert response_get.status_code == 404

    response_get_data = response_get.json()
    assert "detail" in response_get_data
    assert response_get_data["detail"] == "Compound not found"


@pytest.mark.parametrize(
    "input_file, mime_type, expected_length, expected_first",
    [
        (
            BLACK_DIR / "compounds.csv",
            "text/csv",
            54,
            {
                "Common Name": "1,3-Benzenedicarboxylic acid",
                "CAS": "121-91-5",
                "EPA Compound ID": "EPA-001",
                "MolLogP": "1.082999945",
                "Source": "EPA",
                "Source Compound Code": "EPA-001",
                "registration_status": "success",
                "registration_error_message": "",
            },
        ),
        (
            BLACK_DIR / "compounds.sdf",
            "chemical/x-mdl-sdfile",
            54,
            {
                "Common Name": "1,3-Benzenedicarboxylic acid",
                "CAS": "121-91-5",
                "EPA Compound ID": "EPA-001",
                "MolLogP": "1.082999945",
                "Source": "EPA",
                "Source Compound Code": "EPA-001",
                "registration_status": "success",
                "registration_error_message": "",
            },
        ),
    ],
)
def test_compounds_input(client, input_file, mime_type, expected_length, expected_first):
    response = preload_entity(
        client,
        "/v1/compounds/",
        input_file,
        mapping_path=None,
        error_handling=enums.ErrorHandlingOptions.reject_row,
        mime_type=mime_type,
    )

    assert response.status_code == 200
    data = response.json()
    assert len(data) == expected_length

    first = data[0]
    for k, v in expected_first.items():
        actual = first.get(k)
        assert actual.strip() == v.strip(), f"Mismatch in field {k}: expected {v!r}, got {actual!r}"


@pytest.mark.parametrize(
    "input_file, output_format, expected_content_snippet",
    [
        (BLACK_DIR / "compounds.csv", enums.OutputFormat.sdf, "$$$$"),
        (BLACK_DIR / "compounds.sdf", enums.OutputFormat.csv, "Common Name,CAS"),
    ],
)
def test_compounds_output(client, input_file, output_format, expected_content_snippet):
    response = preload_entity(
        client,
        "/v1/compounds/",
        input_file,
        mapping_path=None,
        error_handling=enums.ErrorHandlingOptions.reject_row,
        output_format=output_format,
    )

    assert response.status_code == 200
    content = response.text
    assert expected_content_snippet in content, f"Expected snippet not found: {expected_content_snippet}"

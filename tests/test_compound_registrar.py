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


class BaseRegistrarTest:
    """Base class for registrar-type entity tests (compounds, batches)."""

    entity_name = None
    expected_properties = None
    expected_first = None
    preload_func = None
    first_entity_fixture_name = None
    allow_put = True

    # --- helpers ---
    def _preload(self, client, file=None, mapping=None, error_handling=enums.ErrorHandlingOptions.reject_row):
        file = file or BLACK_DIR / f"{self.entity_name}.csv"
        return self.preload_func(client, file, mapping, error_handling)

    def _get_entities(self, client):
        resp = client.get(f"/v1/{self.entity_name}/")
        assert resp.status_code == 200
        return resp.json()

    def _get_first_entity(self, request):
        return request.getfixturevalue(self.first_entity_fixture_name)

    def test_register_without_mapping(self, client, preload_schema):
        response = self._preload(client)
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 54

        entities = self._get_entities(client)

        def assert_properties(entity, expected_props, index):
            properties = entity.get("properties", [])
            assert len(properties) == len(expected_props), (
                f"[Entity {index}] Expected {len(expected_props)} properties, got {len(properties)}"
            )
            names = {p["name"] for p in properties}
            assert names == expected_props, f"[Entity {index}] Property names mismatch: {names} != {expected_props}"

        assert_properties(entities[0], self.expected_properties, 0)
        assert_properties(entities[8], self.expected_properties - {"MolLogP"}, 8)

    @pytest.mark.skip(reason="No test datasets contain invalid records to validate 'reject all' behaviour.")
    def test_register_reject_all(self, client, preload_schema):
        response = self._preload(
            client,
            mapping=BLACK_DIR / f"{self.entity_name}_mapping.json",
            error_handling=enums.ErrorHandlingOptions.reject_all,
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
    def test_register_reject_row(self, client, preload_schema):
        response = self._preload(client, mapping=BLACK_DIR / f"{self.entity_name}_mapping.json")
        assert response.status_code == 200
        data = response.json()["data"]

        assert len(data) == 54
        assert data[8]["registration_status"] == "failed"
        assert data[8]["registration_error_message"] == "400: Invalid SMILES string"
        assert data[9]["registration_status"] == "success"
        assert data[9]["registration_error_message"] is None

    def test_get_list(self, client, preload_schema):
        response = self._preload(
            client,
            mapping=BLACK_DIR / f"{self.entity_name}_mapping.json",
            error_handling=enums.ErrorHandlingOptions.reject_all,
        )
        assert response.status_code == 200
        entities = self._get_entities(client)

        first = entities[0]
        assert first["id"] == 1
        assert first["canonical_smiles"] == "O=C(O)c1cccc(C(=O)O)c1"
        assert first["inchikey"] == "QQVIHTHCMHWDBS-UHFFFAOYSA-N"

        assert "properties" in first
        assert isinstance(first["properties"], list)
        props = {p["name"]: p for p in first["properties"]}
        for k, v in self.expected_first.items():
            if isinstance(v, float):
                assert abs(props[k]["value_num"] - v) < 1e-3
            else:
                assert props[k]["value_string"].strip() == str(v).strip()

    def test_get_by_any_synonym(self, client, preload_schema, request):
        first, synonyms = self._get_first_entity(request)
        for prop in synonyms:
            synonym_name = prop["name"]
            synonym_value = prop["value_string"]

            resp_val = client.get(f"/v1/{self.entity_name}?property_value={synonym_value}")
            assert resp_val.status_code == 200
            assert resp_val.json()["id"] == first["id"]

            resp_name = client.get(
                f"/v1/{self.entity_name}?property_value={synonym_value}&property_name={synonym_name}"
            )
            assert resp_name.status_code == 200
            assert resp_name.json()["id"] == first["id"]

    def test_get_properties(self, client, preload_schema, request):
        first, synonyms = self._get_first_entity(request)
        for prop in synonyms:
            resp = client.get(f"/v1/{self.entity_name}/properties?property_value={prop['value_string']}")
            assert resp.status_code == 200
            props = resp.json()
            returned_names = {p["name"] for p in props}
            original_names = {p["name"] for p in first["properties"]}
            assert returned_names == original_names

    def test_get_synonyms(self, client, preload_schema, request):
        _, synonyms = self._get_first_entity(request)
        for prop in synonyms:
            resp = client.get(f"/v1/{self.entity_name}/synonyms?property_value={prop['value_string']}")
            assert resp.status_code == 200
            props = resp.json()
            assert all(p["semantic_type_id"] == 1 for p in props), "Non-synonym returned"

    @pytest.mark.parametrize("update_payload", [{"is_archived": True}, {"canonical_smiles": "CCC"}])
    def test_put_entity(self, client, preload_schema, request, update_payload):
        if not self.allow_put:
            pytest.skip(f"PUT endpoint not implemented for {self.entity_name}")
        first, _ = self._get_first_entity(request)
        corporate_id = next(p["value_string"] for p in first["properties"] if "corporate" in p["name"])
        response = client.put(f"/v1/{self.entity_name}/{corporate_id}", json=update_payload)
        assert response.status_code == 200
        data = response.json()
        for k, v in update_payload.items():
            assert data[k] == v, f"Expected {k} to be {v}, but got {data[k]}"

    def test_delete_entity(self, client, preload_schema, request):
        first, _ = self._get_first_entity(request)
        corporate_id = next(p["value_string"] for p in first["properties"] if "corporate" in p["name"])
        response_delete = client.delete(f"/v1/{self.entity_name}/{corporate_id}")
        assert response_delete.status_code == 200

        response_get = client.get(f"/v1/{self.entity_name}/properties?property_value={corporate_id}")
        assert response_get.status_code == 404
        assert "detail" in response_get.json()

    def test_input_files(self, client):
        test_cases = [
            ("csv", "text/csv"),
            ("sdf", "chemical/x-mdl-sdfile"),
        ]

        for ext, mime_type in test_cases:
            input_file = BLACK_DIR / f"{self.entity_name}.{ext}"
            response = preload_entity(
                client,
                f"/v1/{self.entity_name}/",
                input_file,
                mapping_path=None,
                error_handling=enums.ErrorHandlingOptions.reject_row,
                mime_type=mime_type,
            )

            assert response.status_code == 200
            data = response.json()
            assert len(data) == 54

            first = data[0]
            for k, v in self.expected_first.items():
                actual = first.get(k)
                assert str(actual).strip() == str(v).strip(), f"Mismatch in field {k}: expected {v!r}, got {actual!r}"

    def test_output_formats(self, client):
        test_cases = [
            ("csv", enums.OutputFormat.sdf, "$$$$"),
            ("sdf", enums.OutputFormat.csv, "Common Name,CAS"),
        ]

        for ext, output_format, expected_snippet in test_cases:
            input_file = BLACK_DIR / f"{self.entity_name}.{ext}"
            response = preload_entity(
                client,
                f"/v1/{self.entity_name}/",
                input_file,
                mapping_path=None,
                error_handling=enums.ErrorHandlingOptions.reject_row,
                output_format=output_format,
            )

            assert response.status_code == 200
            content = response.text
            assert expected_snippet in content, f"Expected snippet not found: {expected_snippet}"


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
    first_entity_fixture_name = "first_compound_with_synonyms"

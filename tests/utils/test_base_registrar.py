import pytest
from app.utils import enums
from tests.conftest import BLACK_DIR, preload_entity


class BaseRegistrarTest:
    """Base class for registrar-type entity tests (compounds, batches)."""

    entity_name = None
    expected_properties = None
    expected_first = None
    preload_func = None
    get_response_model = None
    first_entity_fixture_name = None
    allow_put = True
    default_entity_count = 54

    # --- helpers ---
    def _preload(
        self,
        client,
        file=None,
        mapping=None,
        error_handling=enums.ErrorHandlingOptions.reject_row,
    ):
        file = file or self._default_file()
        return self.preload_func(client, file, mapping, error_handling)

    def _default_file(self):
        return BLACK_DIR / f"{self.entity_name}.csv"

    def _mapping_path(self):
        return BLACK_DIR / f"{self.entity_name}_mapping.json"

    def _get_entities(self, client):
        resp = client.get(f"/v1/{self.entity_name}/")
        assert resp.status_code == 200
        return resp.json()

    def _get_first_entity(self, request):
        return request.getfixturevalue(self.first_entity_fixture_name)

    def _get_corporate_id(self, entity):
        return next(p["value_string"] for p in entity["properties"] if "corporate" in p["name"])

    # --- tests ---
    def test_register_without_mapping(self, client, preload_schema):
        response = self._preload(client)
        assert response.status_code == 200
        data = response.json()
        assert len(data) == self.default_entity_count

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
            mapping=self._mapping_path(),
            error_handling=enums.ErrorHandlingOptions.reject_all,
        )
        assert response.status_code == 400
        result = response.json()["detail"]
        assert result["status"] == "Success"

        data = result["data"]
        assert len(data) == self.default_entity_count
        assert data[8]["registration_status"] == "failed"
        assert data[8]["registration_error_message"] == "400: Invalid SMILES string"
        for item in data[9:]:
            assert item["registration_status"] == "not_processed"
            assert item["registration_error_message"] is None

    @pytest.mark.skip(reason="No test datasets contain invalid records to validate 'reject row' behaviour.")
    def test_register_reject_row(self, client, preload_schema):
        response = self._preload(client, mapping=self._mapping_path())
        assert response.status_code == 200
        data = response.json()["data"]

        assert len(data) == self.default_entity_count
        assert data[8]["registration_status"] == "failed"
        assert data[8]["registration_error_message"] == "400: Invalid SMILES string"
        assert data[9]["registration_status"] == "success"
        assert data[9]["registration_error_message"] is None

    def test_get_list(self, client, preload_schema):
        response = self._preload(
            client,
            mapping=self._mapping_path(),
            error_handling=enums.ErrorHandlingOptions.reject_all,
        )
        assert response.status_code == 200
        entities = self._get_entities(client)
        assert len(entities) == self.default_entity_count

        for e in entities:
            self.get_response_model(**e)

    def test_get_by_any_synonym(self, client, preload_schema, request):
        first_entity, synonyms = self._get_first_entity(request)
        for prop in synonyms:
            synonym_name = prop["name"]
            synonym_value = prop["value_string"]

            resp_val = client.get(f"/v1/{self.entity_name}?property_value={synonym_value}")
            assert resp_val.status_code == 200
            assert resp_val.json()["id"] == first_entity["id"]

            resp_name = client.get(
                f"/v1/{self.entity_name}?property_value={synonym_value}&property_name={synonym_name}"
            )
            assert resp_name.status_code == 200
            assert resp_name.json()["id"] == first_entity["id"]

    def test_get_properties(self, client, preload_schema, request):
        first_entity, synonyms = self._get_first_entity(request)
        for prop in synonyms:
            resp = client.get(f"/v1/{self.entity_name}/properties?property_value={prop['value_string']}")
            assert resp.status_code == 200
            props = resp.json()
            returned_names = {p["name"] for p in props}
            original_names = {p["name"] for p in first_entity["properties"]}
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
        first_entity, _ = self._get_first_entity(request)
        corporate_id = self._get_corporate_id(first_entity)
        response = client.put(f"/v1/{self.entity_name}/{corporate_id}", json=update_payload)
        assert response.status_code == 200
        data = response.json()
        for k, v in update_payload.items():
            assert data[k] == v, f"Expected {k} to be {v}, but got {data[k]}"

    def test_delete_entity(self, client, preload_schema, request):
        first_entity, _ = self._get_first_entity(request)
        corporate_id = self._get_corporate_id(first_entity)
        response_delete = client.delete(f"/v1/{self.entity_name}/{corporate_id}")
        assert response_delete.status_code == 200

        response_get = client.get(f"/v1/{self.entity_name}/properties?property_value={corporate_id}")
        assert response_get.status_code == 404
        assert "detail" in response_get.json()

    @pytest.mark.parametrize(
        "ext, mime_type",
        [
            ("csv", "text/csv"),
            ("sdf", "chemical/x-mdl-sdfile"),
        ],
    )
    def test_input_files(self, client, ext, mime_type):
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
        assert len(data) == self.default_entity_count

        first = data[0]
        for k, v in self.expected_first.items():
            actual = first.get(k)
            assert str(actual).strip() == str(v).strip(), f"Mismatch in field {k}: expected {v!r}, got {actual!r}"

    @pytest.mark.parametrize(
        "ext, output_format, expected_snippet",
        [
            ("csv", enums.OutputFormat.sdf, "$$$$"),
            ("sdf", enums.OutputFormat.csv, "Common Name,CAS"),
        ],
    )
    def test_output_formats(self, client, ext, output_format, expected_snippet):
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

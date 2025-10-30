import pytest
from app.utils import enums
from tests.conftest import BLACK_DIR, preload_entity


class BaseRegistrarTest:
    """Shared logic for registrar-type entity tests."""

    entity_name = None
    expected_properties = None
    expected_first = None
    preload_func = None
    get_response_model = None
    first_entity_fixture_name = None
    default_entity_count = 54

    # --- helpers ---
    def _preload(self, client, file=None, mapping=None, error_handling=enums.ErrorHandlingOptions.reject_row):
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

    # --- common tests ---
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

    def test_get_list(self, client, preload_schema):
        response = self._preload(
            client, mapping=self._mapping_path(), error_handling=enums.ErrorHandlingOptions.reject_all
        )
        assert response.status_code == 200
        entities = self._get_entities(client)

        assert len(entities) == self.default_entity_count
        for e in entities:
            self.get_response_model(**e)

    @pytest.mark.parametrize(
        "ext, mime_type",
        [
            ("csv", "text/csv"),
            ("sdf", "chemical/x-mdl-sdfile"),
        ],
    )
    def test_input_files(self, client, ext, mime_type):
        input_file = BLACK_DIR / f"{self.entity_name}.{ext}"
        if not input_file.exists():
            return

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
        if not input_file.exists():
            return

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

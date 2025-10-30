import pytest


class SynonymTestsMixin:
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


class CRUDTestsMixin:
    @pytest.mark.parametrize("update_payload", [{"is_archived": True}, {"canonical_smiles": "CCC"}])
    def test_put_entity(self, client, preload_schema, request, update_payload):
        first_entity, _ = self._get_first_entity(request)
        corporate_id = self._get_corporate_id(first_entity)
        response = client.put(f"/v1/{self.entity_name}/{corporate_id}", json=update_payload)
        assert response.status_code == 200
        data = response.json()
        for k, v in update_payload.items():
            assert data[k] == v, f"Expected {k}={v}, got {data[k]}"

    def test_delete_entity(self, client, preload_schema, request):
        first_entity, _ = self._get_first_entity(request)
        corporate_id = self._get_corporate_id(first_entity)
        response_delete = client.delete(f"/v1/{self.entity_name}/{corporate_id}")
        assert response_delete.status_code == 200

        response_get = client.get(f"/v1/{self.entity_name}/properties?property_value={corporate_id}")
        assert response_get.status_code == 404
        assert "detail" in response_get.json()

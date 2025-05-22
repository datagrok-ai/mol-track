from tests.conftest import synonym_type_data


def test_create_synonym_type(client, synonym_type):
    assert synonym_type["name"] == synonym_type_data["name"]
    assert synonym_type["synonym_level"] == synonym_type_data["synonym_level"]
    assert synonym_type["pattern"] == synonym_type_data["pattern"]
    assert "id" in synonym_type


def test_list_synonym_types(client, synonym_type):
    """Test listing all synonym types"""
    response = client.get("/synonym-types/")
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) > 0
    for synonym_type in data:
        assert "id" in synonym_type
        assert "name" in synonym_type
        assert "synonym_level" in synonym_type
        assert "pattern" in synonym_type

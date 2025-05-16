def test_create_synonym_type(client):
    """Test creating a synonym type"""
    synonym_type_data = {
        "name": "CAS Number",
        "synonym_level": "COMPOUND",
        "pattern": r"\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b"
    }
    response = client.post("/synonym-types/", json=synonym_type_data)
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == synonym_type_data["name"]
    assert data["synonym_level"] == synonym_type_data["synonym_level"]
    assert data["pattern"] == synonym_type_data["pattern"]
    assert "id" in data

def test_list_synonym_types(client):
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
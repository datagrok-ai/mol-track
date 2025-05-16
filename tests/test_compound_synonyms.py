def test_create_compound_synonym(client):
    """Test creating a compound synonym"""
    # Create a synonym type
    synonym_type_data = {
        "name": "CAS Number",
        "synonym_level": "COMPOUND",
        "pattern": r"\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b"
    }
    synonym_type = client.post("/synonym-types/", json=synonym_type_data).json()

    # Create a compound synonym
    compound_synonym_data = {
        "compound_id": 1,  # Replace with a valid compound ID
        "synonym_type_id": synonym_type["id"],
        "synonym_value": "12345-67-8"
    }
    response = client.post("/compound-synonyms/", json=compound_synonym_data)
    assert response.status_code == 200
    data = response.json()
    assert data["compound_id"] == compound_synonym_data["compound_id"]
    assert data["synonym_type_id"] == compound_synonym_data["synonym_type_id"]
    assert data["synonym_value"] == compound_synonym_data["synonym_value"]
    assert "id" in data

def test_create_invalid_compound_synonym(client):
    """Test creating a compound synonym with an invalid value"""
    # Create a synonym type
    synonym_type_data = {
        "name": "CAS Number",
        "synonym_level": "COMPOUND",
        "pattern": r"\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b"
    }
    synonym_type = client.post("/synonym-types/", json=synonym_type_data).json()

    # Attempt to create a compound synonym with an invalid value
    invalid_compound_synonym_data = {
        "compound_id": 1,  # Replace with a valid compound ID
        "synonym_type_id": synonym_type["id"],
        "synonym_value": "invalid-value"
    }
    response = client.post("/compound-synonyms/", json=invalid_compound_synonym_data)
    assert response.status_code == 422
    assert "does not match the required pattern" in response.json()["detail"]
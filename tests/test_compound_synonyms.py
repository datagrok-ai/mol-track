def test_create_compound_synonym(client, synonym_type, compound):
    # Create a compound synonym
    compound_synonym_data = {
        "compound_id": compound["id"],
        "synonym_type_id": synonym_type["id"],
        "synonym_value": "ABC-1234",
    }
    response = client.post("/compound-synonyms/", json=compound_synonym_data)
    assert response.status_code == 200
    data = response.json()
    assert data["compound_id"] == compound_synonym_data["compound_id"]
    assert data["synonym_type_id"] == compound_synonym_data["synonym_type_id"]
    assert data["synonym_value"] == compound_synonym_data["synonym_value"]
    assert "id" in data


def test_create_invalid_compound_synonym(client, synonym_type, compound):
    """Test creating a compound synonym with an invalid value"""
    invalid_compound_synonym_data = {
        "compound_id": compound["id"],  # Replace with a valid compound ID
        "synonym_type_id": synonym_type["id"],
        "synonym_value": "12345-67-8",
    }
    response = client.post("/compound-synonyms/", json=invalid_compound_synonym_data)
    assert response.status_code == 422
    assert "does not match the required pattern" in response.json()["detail"]

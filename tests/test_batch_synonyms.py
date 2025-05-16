def test_create_batch_synonym(client):
    """Test creating a batch synonym"""
    # Create a synonym type
    synonym_type_data = {
        "name": "Batch Code",
        "synonym_level": "BATCH",
        "pattern": r"[A-Z]{3}-\d{4}"
    }
    synonym_type = client.post("/synonym-types/", json=synonym_type_data).json()

    # Create a batch synonym
    batch_synonym_data = {
        "batch_id": 1,  # Replace with a valid batch ID
        "synonym_type_id": synonym_type["id"],
        "synonym_value": "ABC-1234"
    }
    response = client.post("/batch-synonyms/", json=batch_synonym_data)
    assert response.status_code == 200
    data = response.json()
    assert data["batch_id"] == batch_synonym_data["batch_id"]
    assert data["synonym_type_id"] == batch_synonym_data["synonym_type_id"]
    assert data["synonym_value"] == batch_synonym_data["synonym_value"]
    assert "id" in data

def test_create_invalid_batch_synonym(client):
    """Test creating a batch synonym with an invalid value"""
    # Create a synonym type
    synonym_type_data = {
        "name": "Batch Code",
        "synonym_level": "BATCH",
        "pattern": r"[A-Z]{3}-\d{4}"
    }
    synonym_type = client.post("/synonym-types/", json=synonym_type_data).json()

    # Attempt to create a batch synonym with an invalid value
    invalid_batch_synonym_data = {
        "batch_id": 1,  # Replace with a valid batch ID
        "synonym_type_id": synonym_type["id"],
        "synonym_value": "invalid-code"
    }
    response = client.post("/batch-synonyms/", json=invalid_batch_synonym_data)
    assert response.status_code == 422
    assert "does not match the required pattern" in response.json()["detail"]
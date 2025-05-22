test_schema_data = {
    "properties": [
        {"name": "MolLogP", "scope": "COMPOUND", "property_class": "CALCULATED", "value_type": "double", "unit": ""}
    ],
    "synonym_types": [
        {"name": "batch_corporate_id", "synonym_level": "BATCH", "pattern": ""},
        {"name": "corporate_id", "synonym_level": "COMPOUND", "pattern": ""},
        {"name": "CAS number", "synonym_level": "COMPOUND", "pattern": r"\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b"},
        {"name": "common name", "synonym_level": "COMPOUND", "pattern": ""},
        {"name": "USAN", "synonym_level": "COMPOUND", "pattern": ""},
    ],
}


def test_preload_schema(client):
    response = client.post("/schema/", json=test_schema_data)
    assert response.status_code == 200

    data = response.json()
    assert data["status"] == "success"

    expected_synonyms = ["batch_corporate_id", "corporate_id", "common name", "CAS number", "USAN"]
    assert "created" in data
    assert "synonym_types" in data["created"]
    assert sorted(data["created"]["synonym_types"]) == sorted(expected_synonyms)

    expected_properties = ["MolLogP"]
    assert "property_types" in data["created"]
    assert data["created"]["property_types"] == expected_properties

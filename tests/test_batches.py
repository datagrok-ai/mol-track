import pytest
from fastapi.testclient import TestClient

# Import common data from conftest.py (fixtures are automatically available)
from tests.conftest import aspirin_smiles, caffeine_smiles

# Test data for batch creation
test_batch_data = {
    "compounds": [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"   # Ibuprofen
    ]
}

# Test data for single compound
test_compound_data = {
    "smiles": aspirin_smiles,  # Aspirin SMILES
    "is_archived": False
}

def test_create_compounds_batch(client):
    """Test creating multiple compounds via the batch API endpoint"""
    response = client.post("/compounds/batch/", json=test_batch_data)
    
    # Check response status code and content
    assert response.status_code == 200
    data = response.json()
    
    # Verify we got the right number of compounds
    assert len(data) == len(test_batch_data["compounds"])
    
    # Verify each compound has the expected fields
    for compound in data:
        assert "id" in compound
        assert "canonical_smiles" in compound
        assert "inchi" in compound
        assert "inchikey" in compound
        assert "is_archived" in compound
        assert isinstance(compound["id"], int)

def test_create_compounds_batch_with_invalid_smiles(client):
    """Test that batch creation fails with invalid SMILES"""
    invalid_batch_data = {
        "compounds": [
            "CC(=O)OC1=CC=CC=C1C(=O)O",  # Valid
            "INVALID_SMILES",             # Invalid
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # Valid
        ]
    }
    
    response = client.post("/compounds/batch/", json=invalid_batch_data)
    
    # Check if the API properly validates the request
    assert response.status_code == 400
    assert "Invalid SMILES" in response.json()["detail"]

def test_create_compounds_batch_with_duplicate(client):
    """Test that batch creation fails when a compound already exists"""
    # First create a single compound
    client.post("/compounds/", json=test_compound_data)
    
    # Then try to create a batch that includes the same compound
    duplicate_batch_data = {
        "compounds": [
            "CC(=O)OC1=CC=CC=C1C(=O)O",  # Already exists (Aspirin)
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # New compound (Caffeine)
        ]
    }
    
    response = client.post("/compounds/batch/", json=duplicate_batch_data)
    
    # Check if the API handles duplicates appropriately
    assert response.status_code == 400
    assert "already exist" in response.json()["detail"].lower() 
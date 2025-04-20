import pytest
from fastapi.testclient import TestClient

# Import common data from conftest.py (fixtures are automatically available)
from tests.conftest import aspirin_smiles, aspirin_smiles_noncanonical, caffeine_smiles

# Test data
test_compound_data = {
    "smiles": aspirin_smiles,  # Aspirin SMILES
    "is_archived": False
}

def test_create_compound(client):
    """Test creating a new compound via the API endpoint"""
    response = client.post("/compounds/", json=test_compound_data)
    
    # Check response status code and content
    assert response.status_code == 200
    data = response.json()
    
    # Verify the returned data contains the expected values
    assert "canonical_smiles" in data
    assert "inchi" in data
    assert "inchikey" in data
    assert data["is_archived"] == test_compound_data["is_archived"]
    
    # Verify the compound has an ID (was saved to the database)
    assert "id" in data
    assert isinstance(data["id"], int)

def test_create_compound_duplicate_smiles(client):
    """Test that creating a compound with a duplicate SMILES string fails appropriately"""
    # First create a compound
    client.post("/compounds/", json=test_compound_data)
    
    # Try to create another compound with the same SMILES
    response = client.post("/compounds/", json=test_compound_data)
    
    # Check if the API handles duplicates appropriately
    assert response.status_code == 400
    assert "already exists" in response.json()["detail"].lower()

def test_create_compound_invalid_data(client):
    """Test creating a compound with invalid data"""
    invalid_data = {
        "smiles": "invalid data"
    }
    response = client.post("/compounds/", json=invalid_data)
    
    # Check if the API properly validates the request
    assert response.status_code == 400  # Unprocessable Entity
    assert "Invalid SMILES string" in response.json()["detail"]

def test_get_compound_after_creation(client):
    """Test that a created compound can be retrieved"""
    # Create a compound
    create_response = client.post("/compounds/", json=test_compound_data)
    compound_id = create_response.json()["id"]
    
    # Then try to get it
    get_response = client.get(f"/compounds/{compound_id}")
    
    # Check response
    assert get_response.status_code == 200
    data = get_response.json()
    
    # Verify the returned data matches what we created
    assert data["id"] == compound_id
    assert "canonical_smiles" in data
    assert "inchi" in data
    assert "inchikey" in data

def test_list_compounds(client):
    """Test listing all compounds"""
   
    client.post("/compounds/", json={"smiles": aspirin_smiles})
    client.post("/compounds/", json={"smiles": caffeine_smiles})
    
    # Get all compounds
    response = client.get("/compounds/")
    
    # Check response
    assert response.status_code == 200
    data = response.json()
    
    # Should have at least 2 compounds
    assert len(data) >= 2
    
    # Check if compounds are in the list (we can't check exact SMILES since they're canonicalized)
    assert len(data) >= 2

def test_update_compound(client):
    """Test updating a compound"""
    # First create a compound
    create_response = client.post("/compounds/", json=test_compound_data)
    compound_id = create_response.json()["id"]
    
    # Update data
    update_data = {
        "is_archived": True
    }
    
    # Update the compound
    update_response = client.put(f"/compounds/{compound_id}", json=update_data)
    
    # Check response
    assert update_response.status_code == 200
    updated_data = update_response.json()
    
    # Verify the updated fields
    assert updated_data["is_archived"] == True
    
    # Verify other fields remain unchanged
    assert "canonical_smiles" in updated_data
    assert "inchi" in updated_data
    assert "inchikey" in updated_data

    # Verify that dates are appropriate
    assert "created_at" in updated_data
    assert "updated_at" in updated_data
    assert updated_data["updated_at"] > updated_data["created_at"]

def test_delete_compound(client):
    """Test deleting a compound"""
    # First create a compound
    create_response = client.post("/compounds/", json=test_compound_data)
    compound_id = create_response.json()["id"]
    
    # Delete the compound
    delete_response = client.delete(f"/compounds/{compound_id}")
    
    # Check response
    assert delete_response.status_code == 200
    
    # Try to get the deleted compound
    get_response = client.get(f"/compounds/{compound_id}")
    
    # Should return 404 Not Found
    assert get_response.status_code == 404 

def test_create_compound_with_equivalent_smiles(client):
    """Test that creating a compound with equivalent but different SMILES representation fails"""
    # First create a compound with canonical SMILES
    client.post("/compounds/", json={"smiles": aspirin_smiles})
    
    # Try to create another compound with non-canonical but chemically equivalent SMILES
    response = client.post("/compounds/", json={"smiles": aspirin_smiles_noncanonical})
    
    # This should fail since they represent the same molecule
    assert response.status_code == 400
    assert "already exists" in response.json()["detail"].lower() 
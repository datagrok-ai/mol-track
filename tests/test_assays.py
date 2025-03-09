import pytest
from fastapi.testclient import TestClient
import sys
import os

# Add the parent directory to sys.path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Fixtures are automatically available from conftest.py

# Test data for properties
test_property_data = {
    "name": "LogP",
    "value_type": "NUMBER",
    "property_class": "CALCULATED",
    "unit": "log units"
}

test_property_data2 = {
    "name": "Solubility",
    "value_type": "NUMBER",
    "property_class": "MEASURED",
    "unit": "mg/mL"
}

# Test data for assays
test_assay_data = {
    "name": "Lipophilicity Assay",
    "description": "Measures the lipophilicity of compounds",
    "property_ids": []  # Will be populated with actual property IDs in each test
}

test_assay_update_data = {
    "name": "Updated Lipophilicity Assay",
    "description": "Updated description for the lipophilicity assay",
    "property_ids": []  # Will be populated with actual property IDs in each test
}

def test_create_property(client):
    """Test creating a property that will be used in assay tests"""
    response = client.post("/properties/", json=test_property_data)
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == test_property_data["name"]
    assert data["value_type"] == test_property_data["value_type"]
    assert data["property_class"] == test_property_data["property_class"]
    assert data["unit"] == test_property_data["unit"]
    assert "id" in data

def test_create_assay(client):
    """Test creating a new assay via the API endpoint"""
    # Create a property first
    prop_response = client.post("/properties/", json=test_property_data)
    assert prop_response.status_code == 200
    prop_data = prop_response.json()
    
    # Create assay data with the new property ID
    assay_data = {
        "name": test_assay_data["name"],
        "description": test_assay_data["description"],
        "property_ids": [prop_data["id"]]
    }
    
    # Now create the assay
    response = client.post("/assays/", json=assay_data)
    
    # Print response for debugging
    print(f"Response status: {response.status_code}")
    print(f"Response body: {response.text}")
    
    # Check response status code and content
    assert response.status_code == 200
    data = response.json()
    
    # Verify the returned data contains the expected values
    assert data["name"] == assay_data["name"]
    assert data["description"] == assay_data["description"]
    
    # Verify the assay has an ID (was saved to the database)
    assert "id" in data
    assert isinstance(data["id"], int)
    
    # Verify the properties were associated correctly
    assert len(data["properties"]) == len(assay_data["property_ids"])
    if assay_data["property_ids"]:
        assert data["properties"][0]["id"] == assay_data["property_ids"][0]

def test_create_assay_with_invalid_property(client):
    """Test creating an assay with an invalid property ID"""
    invalid_assay_data = {
        "name": "Invalid Assay",
        "description": "This assay has an invalid property ID",
        "property_ids": [9999]  # Assuming this ID doesn't exist
    }
    
    response = client.post("/assays/", json=invalid_assay_data)
    
    # Check if the API properly validates the request
    assert response.status_code == 404
    assert "Property with ID" in response.json()["detail"]
    assert "not found" in response.json()["detail"]

def test_get_assay_after_creation(client):
    """Test that a created assay can be retrieved"""
    # Create a property first
    prop_response = client.post("/properties/", json=test_property_data)
    assert prop_response.status_code == 200
    prop_data = prop_response.json()
    
    # Create assay data with the new property ID
    assay_data = {
        "name": test_assay_data["name"],
        "description": test_assay_data["description"],
        "property_ids": [prop_data["id"]]
    }
    
    # Create an assay
    create_response = client.post("/assays/", json=assay_data)
    assert create_response.status_code == 200
    assay_id = create_response.json()["id"]
    
    # Retrieve the assay
    response = client.get(f"/assays/{assay_id}")
    
    # Check response
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == assay_id
    assert data["name"] == assay_data["name"]
    assert data["description"] == assay_data["description"]

def test_list_assays(client):
    """Test listing all assays"""
    # Create an assay if not already created
    client.post("/assays/", json={"name": "List Test Assay", "description": "For testing list endpoint"})
    
    # Get all assays
    response = client.get("/assays/")
    
    # Check response
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) > 0
    
    # Check that each item has the expected structure
    for assay in data:
        assert "id" in assay
        assert "name" in assay
        assert "description" in assay
        assert "properties" in assay

def test_update_assay(client):
    """Test updating an assay"""
    # Create first property
    prop1_response = client.post("/properties/", json=test_property_data)
    assert prop1_response.status_code == 200
    prop1_data = prop1_response.json()
    
    # Create second property for update
    prop2_response = client.post("/properties/", json=test_property_data2)
    assert prop2_response.status_code == 200
    prop2_data = prop2_response.json()
    
    # Create assay data with the first property
    assay_data = {
        "name": test_assay_data["name"],
        "description": test_assay_data["description"],
        "property_ids": [prop1_data["id"]]
    }
    
    # Create update data with the second property
    update_data = {
        "name": test_assay_update_data["name"],
        "description": test_assay_update_data["description"],
        "property_ids": [prop2_data["id"]]
    }
    
    # First create an assay
    create_response = client.post("/assays/", json=assay_data)
    assert create_response.status_code == 200
    assay_id = create_response.json()["id"]
    
    # Update the assay
    response = client.put(f"/assays/{assay_id}", json=update_data)
    
    # Check response
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == assay_id
    assert data["name"] == update_data["name"]
    assert data["description"] == update_data["description"]
    
    # Check that properties were updated
    if update_data["property_ids"]:
        assert len(data["properties"]) == len(update_data["property_ids"])
        property_ids = [prop["id"] for prop in data["properties"]]
        for prop_id in update_data["property_ids"]:
            assert prop_id in property_ids

def test_delete_assay(client):
    """Test deleting an assay"""
    # First create an assay
    create_response = client.post("/assays/", json={"name": "Delete Test Assay", "description": "For testing delete endpoint"})
    assay_id = create_response.json()["id"]
    
    # Delete the assay
    response = client.delete(f"/assays/{assay_id}")
    
    # Check response
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == assay_id
    
    # Verify the assay was deleted
    get_response = client.get(f"/assays/{assay_id}")
    assert get_response.status_code == 404 
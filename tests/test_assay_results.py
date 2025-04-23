import pytest
from fastapi.testclient import TestClient
import sys
import os

# Add the parent directory to sys.path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Fixtures are automatically available from conftest.py

# Define individual property test data
ic50prop = {
    "name": "IC50",
    "value_type": "NUMBER",
    "property_class": "MEASURED",
    "unit": "nM"
}

solProp = {
    "name": "Solubility",
    "value_type": "NUMBER",
    "property_class": "MEASURED",
    "unit": "mg/mL"
}

activityProp = {
    "name": "Activity Score",
    "value_type": "NUMBER",
    "property_class": "CALCULATED",
    "unit": ""
}

# Test properties list
test_properties = [ic50prop, solProp, activityProp]

# Helper functions for creating test data
def create_test_property(client, property_data):
    """Helper function to create a property and return its data"""
    response = client.post("/properties/", json=property_data)
    assert response.status_code == 200
    return response.json()

def create_test_properties(client, property_data_list):
    """Helper function to create multiple properties and return their IDs"""
    property_ids = []
    for property_data in property_data_list:
        response = client.post("/properties/", json=property_data)
        assert response.status_code == 200
        property_ids.append(response.json()["id"])
    return property_ids

def create_test_assay_type(client, name, description, property_ids):
    """Helper function to create an assay type and return its data"""
    assay_type_data = {
        "name": name,
        "description": description,
        "property_ids": property_ids
    }
    response = client.post("/assay-types/", json=assay_type_data)
    assert response.status_code == 200
    return response.json()

def create_test_assay(client, name, description, assay_type_id, property_ids=None):
    """Helper function to create an assay and return its data"""
    if property_ids is None:
        property_ids = []
    
    assay_data = {
        "name": name,
        "description": description,
        "assay_type_id": assay_type_id,
        "property_ids": property_ids
    }
    response = client.post("/assays/", json=assay_data)
    assert response.status_code == 200
    return response.json()

def create_test_compound(client, smiles="CC(=O)OC1=CC=CC=C1C(=O)O"):
    """Helper function to create a test compound"""
    compound_data = {
        "smiles": smiles,
        "original_molfile": ""
    }
    response = client.post("/compounds/", json=compound_data)
    assert response.status_code == 200
    return response.json()

def create_test_batch(client, compound_id, batch_number="BATCH-001"):
    """Helper function to create a test batch"""
    batch_data = {
        "compound_id": compound_id,
        "batch_number": batch_number,
        "amount": 100.0,
        "amount_unit": "mg",
        "purity": 98.5
    }
    response = client.post("/batches/", json=batch_data)
    assert response.status_code == 200
    return response.json()

def create_test_assay_result(client, assay_id, batch_id, property_id, result_value):
    """Helper function to create a single assay result"""
    assay_result_data = {
        "assay_id": assay_id,
        "batch_id": batch_id,
        "property_id": property_id,
        "result_value": result_value
    }
    response = client.post("/assay-results/", json=assay_result_data)
    assert response.status_code == 200
    return response.json()

def create_test_batch_assay_results(client, assay_id, batch_id, measurements):
    """Helper function to create multiple assay results for a batch"""
    batch_results_data = {
        "assay_id": assay_id,
        "batch_id": batch_id,
        "measurements": measurements
    }
    response = client.post("/batch-assay-results/", json=batch_results_data)
    assert response.status_code == 200
    return response.json()

def test_create_assay_result(client):
    """Test creating a single assay result"""
    # Create a property
    property = create_test_property(client, ic50prop)
    
    # Create an assay type with the property
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", [property["id"]])
    
    # Create an assay
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], [property["id"]])
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create a single assay result
    result_value = 5.2
    assay_result = create_test_assay_result(client, assay["id"], batch["id"], property["id"], result_value)
    
    # Verify data
    assert assay_result["assay_id"] == assay["id"]
    assert assay_result["batch_id"] == batch["id"]
    assert assay_result["property_id"] == property["id"]
    assert assay_result["result_value"] == result_value
    assert "id" in assay_result

def test_create_batch_assay_results(client):
    """Test creating multiple assay results for a batch at once"""
    # Create properties
    properties = [create_test_property(client, prop) for prop in test_properties]
    property_ids = [prop["id"] for prop in properties]
    property_names = [prop["name"] for prop in properties]
    
    # Create an assay type with properties
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", property_ids)
    
    # Create an assay with all properties
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], property_ids)
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create measurements for all properties
    measurements = {
        property_names[0]: 5.2,
        property_names[1]: 0.89,
        property_names[2]: 75
    }
    
    # Create batch assay results
    results = create_test_batch_assay_results(client, assay["id"], batch["id"], measurements)
    
    # Verify data
    assert results["assay_id"] == assay["id"]
    assert results["batch_id"] == batch["id"]
    assert results["assay_name"] == assay["name"]
    assert results["measurements"] == measurements

def test_create_assay_result_with_invalid_property(client):
    """Test creating an assay result with invalid property"""
    # Create a property
    property = create_test_property(client, ic50prop)
    
    # Create an assay type with the property
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", [property["id"]])
    
    # Create an assay but don't associate the property with it
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], [])
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Try to create an assay result with the property not associated with the assay
    assay_result_data = {
        "assay_id": assay["id"],
        "batch_id": batch["id"],
        "property_id": property["id"],
        "result_value": 5.2
    }
    
    response = client.post("/assay-results/", json=assay_result_data)
    
    # Should fail
    assert response.status_code == 400
    assert "not associated with assay" in response.json()["detail"]

def test_batch_assay_results_with_invalid_property(client):
    """Test creating batch assay results with invalid property name"""
    # Create properties
    properties = [create_test_property(client, prop) for prop in test_properties[:1]]
    property_ids = [prop["id"] for prop in properties]
    
    # Create an assay type and assay with only one property
    assay_type = create_test_assay_type(client, "Limited Assay", "Assay with limited properties", property_ids)
    assay = create_test_assay(client, "Limited Assay", "Limited assay", assay_type["id"], property_ids)
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create measurements with invalid property
    measurements = {
        "InvalidProperty": 5.2
    }
    
    # Try to create batch assay results
    batch_results_data = {
        "assay_id": assay["id"],
        "batch_id": batch["id"],
        "measurements": measurements
    }
    
    response = client.post("/batch-assay-results/", json=batch_results_data)
    
    # Should fail with 400 Bad Request
    assert response.status_code == 400
    assert "not associated with assay" in response.json()["detail"]

def test_get_assay_result(client):
    """Test retrieving a single assay result"""
    # Create a property
    property = create_test_property(client, ic50prop)
    
    # Create an assay type with the property
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", [property["id"]])
    
    # Create an assay
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], [property["id"]])
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create a single assay result
    result_value = 5.2
    assay_result = create_test_assay_result(client, assay["id"], batch["id"], property["id"], result_value)
    
    # Get the assay result
    response = client.get(f"/assay-results/{assay_result['id']}")
    assert response.status_code == 200
    
    get_result = response.json()
    assert get_result["id"] == assay_result["id"]
    assert get_result["assay_id"] == assay["id"]
    assert get_result["batch_id"] == batch["id"]
    assert get_result["property_id"] == property["id"]
    assert get_result["result_value"] == result_value

def test_get_batch_assay_results(client):
    """Test retrieving all assay results for a batch"""
    # Create properties
    properties = [create_test_property(client, prop) for prop in test_properties]
    property_ids = [prop["id"] for prop in properties]
    property_names = [prop["name"] for prop in properties]
    
    # Create different assay types
    assay_type1 = create_test_assay_type(client, "Assay Type 1", "First assay type", property_ids[:2])
    assay_type2 = create_test_assay_type(client, "Assay Type 2", "Second assay type", property_ids[1:])
    
    # Create assays
    assay1 = create_test_assay(client, "Assay 1", "First assay", assay_type1["id"], property_ids[:2])
    assay2 = create_test_assay(client, "Assay 2", "Second assay", assay_type2["id"], property_ids[1:])
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create measurements for each assay
    measurements1 = {
        property_names[0]: 5.2,
        property_names[1]: 0.89
    }
    
    measurements2 = {
        property_names[1]: 1.45,
        property_names[2]: 85
    }
    
    # Create batch assay results
    create_test_batch_assay_results(client, assay1["id"], batch["id"], measurements1)
    create_test_batch_assay_results(client, assay2["id"], batch["id"], measurements2)
    
    # Get all assay results for the batch
    response = client.get(f"/batches/{batch['id']}/assay-results")
    assert response.status_code == 200
    
    results = response.json()
    assert len(results) == 2
    
    # Verify both assays are included
    assay_ids = [r["assay_id"] for r in results]
    assert assay1["id"] in assay_ids
    assert assay2["id"] in assay_ids
    
    # Verify measurements for each assay
    for result in results:
        if result["assay_id"] == assay1["id"]:
            assert result["measurements"][property_names[0]] == measurements1[property_names[0]]
            assert result["measurements"][property_names[1]] == measurements1[property_names[1]]
        elif result["assay_id"] == assay2["id"]:
            assert result["measurements"][property_names[1]] == measurements2[property_names[1]]
            assert result["measurements"][property_names[2]] == measurements2[property_names[2]]

def test_update_assay_result(client):
    """Test updating a single assay result"""
    # Create a property
    property = create_test_property(client, ic50prop)
    
    # Create an assay type with the property
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", [property["id"]])
    
    # Create an assay
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], [property["id"]])
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create a single assay result
    initial_value = 5.2
    assay_result = create_test_assay_result(client, assay["id"], batch["id"], property["id"], initial_value)
    
    # Update the assay result
    updated_value = 6.8
    update_data = {
        "result_value": updated_value
    }
    
    response = client.put(f"/assay-results/{assay_result['id']}", json=update_data)
    assert response.status_code == 200
    
    updated_result = response.json()
    assert updated_result["id"] == assay_result["id"]
    assert updated_result["result_value"] == updated_value
    
    # Verify the update persisted
    get_response = client.get(f"/assay-results/{assay_result['id']}")
    assert get_response.status_code == 200
    get_result = get_response.json()
    assert get_result["result_value"] == updated_value

def test_update_batch_assay_results(client):
    """Test updating multiple assay results for a batch"""
    # Create properties
    properties = [create_test_property(client, prop) for prop in test_properties]
    property_ids = [prop["id"] for prop in properties]
    property_names = [prop["name"] for prop in properties]
    
    # Create an assay type with properties
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", property_ids)
    
    # Create an assay
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], property_ids)
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create initial measurements
    initial_measurements = {
        property_names[0]: 5.2,
        property_names[1]: 0.89,
        property_names[2]: 75
    }
    
    # Create batch assay results
    create_test_batch_assay_results(client, assay["id"], batch["id"], initial_measurements)
    
    # Update measurements
    updated_measurements = {
        property_names[0]: 6.1,
        property_names[1]: 1.05,
        property_names[2]: 82
    }
    
    # Update the batch assay results
    response = client.put(
        f"/assays/{assay['id']}/batches/{batch['id']}/results", 
        json=updated_measurements
    )
    assert response.status_code == 200
    
    updated_results = response.json()
    assert updated_results["assay_id"] == assay["id"]
    assert updated_results["batch_id"] == batch["id"]
    assert updated_results["measurements"] == updated_measurements
    
    # Verify the update persisted
    get_response = client.get(f"/batches/{batch['id']}/assay-results")
    assert get_response.status_code == 200
    
    results = get_response.json()
    for result in results:
        if result["assay_id"] == assay["id"]:
            assert result["measurements"][property_names[0]] == updated_measurements[property_names[0]]
            assert result["measurements"][property_names[1]] == updated_measurements[property_names[1]]
            assert result["measurements"][property_names[2]] == updated_measurements[property_names[2]]

def test_delete_assay_result(client):
    """Test deleting a single assay result"""
    # Create a property
    property = create_test_property(client, ic50prop)
    
    # Create an assay type with the property
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", [property["id"]])
    
    # Create an assay
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], [property["id"]])
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create a single assay result
    assay_result = create_test_assay_result(client, assay["id"], batch["id"], property["id"], 5.2)
    
    # Delete the assay result
    response = client.delete(f"/assay-results/{assay_result['id']}")
    assert response.status_code == 200
    
    # Verify it was deleted
    get_response = client.get(f"/assay-results/{assay_result['id']}")
    assert get_response.status_code == 404

def test_delete_batch_assay_results(client):
    """Test deleting all assay results for a batch and assay"""
    # Create properties
    properties = [create_test_property(client, prop) for prop in test_properties]
    property_ids = [prop["id"] for prop in properties]
    property_names = [prop["name"] for prop in properties]
    
    # Create an assay type with properties
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", property_ids)
    
    # Create an assay
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], property_ids)
    
    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])
    
    # Create measurements
    measurements = {
        property_names[0]: 5.2,
        property_names[1]: 0.89,
        property_names[2]: 75
    }
    
    # Create batch assay results
    create_test_batch_assay_results(client, assay["id"], batch["id"], measurements)
    
    # Delete all assay results for the batch and assay
    response = client.delete(f"/assays/{assay['id']}/batches/{batch['id']}/results")
    assert response.status_code == 200
    
    # Verify they were deleted
    get_response = client.get(f"/batches/{batch['id']}/assay-results")
    assert get_response.status_code == 200
    assert len(get_response.json()) == 0 
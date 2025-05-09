import sys
import os

# Add the parent directory to sys.path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Fixtures are automatically available from conftest.py

# Define individual property test data
ic50prop = {
    "name": "IC50",
    "value_type": "double",
    "property_class": "MEASURED",
    "unit": "nM",
}

solProp = {
    "name": "Solubility",
    "value_type": "double",
    "property_class": "MEASURED",
    "unit": "mg/mL",
}

activityProp = {
    "name": "Activity Score",
    "value_type": "double",
    "property_class": "CALCULATED",
    "unit": "",
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
        "property_ids": property_ids,
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
        "property_ids": property_ids,
    }
    response = client.post("/assays/", json=assay_data)
    assert response.status_code == 200
    return response.json()


def create_test_compound(client, smiles="CC(=O)OC1=CC=CC=C1C(=O)O"):
    """Helper function to create a test compound"""
    compound_data = {
        "smiles": smiles,
        "original_molfile": "",
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
        "purity": 98.5,
    }
    response = client.post("/batches/", json=batch_data)
    assert response.status_code == 200
    return response.json()


def create_test_assay_result(client, assay_id, batch_id, property_id, result_value):
    """Helper function to create a single assay result"""
    # Get the property to determine its type
    response = client.get(f"/properties/{property_id}")
    assert response.status_code == 200
    property_data = response.json()
    property_type = property_data["value_type"]

    # Create the appropriate payload based on property type
    assay_result_data = {
        "assay_id": assay_id,
        "batch_id": batch_id,
        "property_id": property_id,
        "value_qualifier": 0,  # 0 for "=" (default)
    }

    # Set the appropriate value field based on property type
    if property_type == "DOUBLE" or property_type == "double" or property_type == "INT" or property_type == "int":
        assay_result_data["value_num"] = result_value
    elif property_type == "STRING" or property_type == "string":
        assay_result_data["value_string"] = str(result_value)
    elif property_type == "BOOL" or property_type == "bool":
        assay_result_data["value_bool"] = bool(result_value)

    response = client.post("/assay-results/", json=assay_result_data)
    assert response.status_code == 200

    result = response.json()

    # For backwards compatibility with the test expectations, add a result_value field
    if "value_num" in result and result["value_num"] is not None:
        result["result_value"] = result["value_num"]
    elif "value_string" in result and result["value_string"] is not None:
        result["result_value"] = result["value_string"]
    elif "value_bool" in result and result["value_bool"] is not None:
        result["result_value"] = result["value_bool"]

    return result


def create_test_batch_assay_results(client, assay_id, batch_id, measurements):
    """Helper function to create multiple assay results for a batch"""
    batch_results_data = {
        "assay_id": assay_id,
        "batch_id": batch_id,
        "measurements": measurements,
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
        property_names[2]: 75,
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
    # Create a property for the assay type
    property1 = create_test_property(client, ic50prop)

    # Create a different property that won't be associated with the assay
    newProp = {
        "name": "Unrelated Property",
        "value_type": "double",
        "property_class": "MEASURED",
        "unit": "units",
    }
    property2 = create_test_property(client, newProp)

    # Create an assay type with property1
    assay_type = create_test_assay_type(client, "Enzyme Assay", "Enzyme activity assay", [property1["id"]])

    # Create an assay with property1
    assay = create_test_assay(client, "Kinase Assay", "Kinase inhibition assay", assay_type["id"], [property1["id"]])

    # Create compound and batch
    compound = create_test_compound(client)
    batch = create_test_batch(client, compound["id"])

    # Try to create an assay result with property2 that is not associated with the assay
    assay_result_data = {
        "assay_id": assay["id"],
        "batch_id": batch["id"],
        "property_id": property2["id"],  # Use the unrelated property
        "value_num": 5.2,
    }

    response = client.post("/assay-results/", json=assay_result_data)

    # Should fail
    assert response.status_code == 400
    # The message should be about property not being associated with the assay type
    error_message = response.json()["detail"]
    assert "not associated with assay" in error_message


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
    measurements = {"InvalidProperty": 5.2}

    # Try to create batch assay results
    batch_results_data = {
        "assay_id": assay["id"],
        "batch_id": batch["id"],
        "measurements": measurements,
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

    # Verify data
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == assay_result["id"]
    assert data["assay_id"] == assay["id"]
    assert data["batch_id"] == batch["id"]
    assert data["property_id"] == property["id"]
    assert data["result_value"] == result_value


def test_get_batch_assay_results(client):
    """Test retrieving all assay results for a batch"""
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
        property_names[2]: 75,
    }

    # Create batch assay results
    create_test_batch_assay_results(client, assay["id"], batch["id"], measurements)

    # Get the batch assay results
    response = client.get(f"/batches/{batch['id']}/assay-results")

    # Verify data
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) == 1  # Should have one assay with results

    # Check the first assay result
    assert data[0]["assay_id"] == assay["id"]
    assert data[0]["batch_id"] == batch["id"]
    assert data[0]["assay_name"] == assay["name"]

    # Check measurements
    for prop_name, value in measurements.items():
        assert prop_name in data[0]["measurements"]
        assert data[0]["measurements"][prop_name] == value

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
    "scope": "COMPOUND"
}

solProp = {
    "name": "Solubility",
    "value_type": "double",
    "property_class": "MEASURED",
    "unit": "mg/mL",
    "scope": "COMPOUND"
}

activityProp = {
    "name": "Activity Score",
    "value_type": "double",
    "property_class": "CALCULATED",
    "unit": "",
    "scope": "COMPOUND"
}

# Test properties list
test_properties = [ic50prop, solProp, activityProp]

test_semantic_type_data = {
    "name": "Absorption",
    "description": "Describes properties related to compound absorption"
}

# Test data for assay types
test_assay_type_data = {
    "name": "Enzyme Inhibition Assay Type",
    "description": "Template for enzyme inhibition assays",
    "property_ids": [],  # Will be populated with actual property IDs in each test
}

# Test data for assays
test_assay_data = {
    "name": "Kinase Inhibition Assay",
    "description": "Measures inhibition of kinase enzymes",
    "assay_type_id": 1,  # Will be populated with an actual assay type ID in each test
    "property_ids": [],  # Will be populated with actual property IDs in each test
}

test_assay_update_data = {
    "name": "Updated Kinase Inhibition Assay",
    "description": "Updated description for the kinase inhibition assay",
    "property_ids": [],  # Will be populated with actual property IDs in each test
}


# Helper functions for creating test data
def create_test_property(client, property_data):
    """Helper function to create a property and return its data"""
    response = client.post("/properties/", json=property_data)
    assert response.status_code == 200
    return response.json()

def create_test_semantic_type(client, semantic_type_data):
    response = client.post("/semantic-types/", json=semantic_type_data)
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

def test_create_semantic_type(client):
    """Test creating a semantic type"""
    data = create_test_semantic_type(client, test_semantic_type_data)
    assert data["name"] == test_semantic_type_data["name"]
    assert data["description"] == test_semantic_type_data["description"]
    assert "id" in data

    for prop in test_properties:
        prop["semantic_type_id"] = data["id"]

def test_create_property(client):
    """Test creating a property"""
    data = create_test_property(client, ic50prop)
    assert data["name"] == ic50prop["name"]
    assert data["value_type"] == ic50prop["value_type"]
    assert data["property_class"] == ic50prop["property_class"]
    assert data["unit"] == ic50prop["unit"]
    assert "id" in data


def test_create_assay(client):
    """Test creating a new assay via the API endpoint"""
    # Create a property
    prop_data = create_test_property(client, ic50prop)

    # Create an assay type with the property
    assay_type = create_test_assay_type(client, "Test Assay Type", "For testing assay creation", [prop_data["id"]])

    # Create an assay
    assay = create_test_assay(
        client, test_assay_data["name"], test_assay_data["description"], assay_type["id"], [prop_data["id"]]
    )

    # Verify the returned data
    assert assay["name"] == test_assay_data["name"]
    assert assay["description"] == test_assay_data["description"]
    assert assay["assay_type_id"] == assay_type["id"]
    assert "id" in assay
    assert isinstance(assay["id"], int)

    # Verify the properties were associated correctly
    assert len(assay["properties"]) == 1
    assert assay["properties"][0]["id"] == prop_data["id"]


def test_create_assay_with_invalid_property(client):
    """Test creating an assay with an invalid property ID"""
    # Create an assay type
    assay_type = create_test_assay_type(client, "Test Assay Type", "For testing with invalid property", [])

    # Try to create an assay with an invalid property ID
    invalid_assay_data = {
        "name": "Invalid Assay",
        "description": "This assay has an invalid property ID",
        "assay_type_id": assay_type["id"],
        "property_ids": [9999],  # Assuming this ID doesn't exist
    }

    response = client.post("/assays/", json=invalid_assay_data)

    # Check if the API properly validates the request
    assert response.status_code == 404
    assert "Property with ID" in response.json()["detail"]
    assert "not found" in response.json()["detail"]


def test_get_assay_after_creation(client):
    """Test that a created assay can be retrieved"""
    # Create a property
    prop_data = create_test_property(client, ic50prop)

    # Create an assay type
    assay_type = create_test_assay_type(client, "Test Assay Type", "For testing assay retrieval", [prop_data["id"]])

    # Create an assay
    assay = create_test_assay(
        client, test_assay_data["name"], test_assay_data["description"], assay_type["id"], [prop_data["id"]]
    )

    # Retrieve the assay
    response = client.get(f"/assays/{assay['id']}")

    # Check response
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == assay["id"]
    assert data["name"] == assay["name"]
    assert data["description"] == assay["description"]
    assert data["assay_type_id"] == assay_type["id"]


def test_list_assays(client):
    """Test listing all assays"""

    # Create an assay type
    assay_type = create_test_assay_type(client, "List Test Assay Type", "For testing list endpoint", [])

    # Create an assay
    create_test_assay(client, "List Test Assay", "For testing list endpoint", assay_type["id"])

    # Get all assays
    response = client.get("/assays/")

    # Check response
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) > 0

    # Check structure of returned items
    for assay in data:
        assert "id" in assay
        assert "name" in assay
        assert "description" in assay
        assert "assay_type_id" in assay
        assert "properties" in assay


def test_create_assay_with_multiple_properties(client):
    """Test creating an assay with multiple properties"""
    # Create properties
    property_ids = create_test_properties(client, test_properties)

    # Create an assay type
    assay_type = create_test_assay_type(
        client, "Multiple Property Assay Type", "For testing multiple properties", property_ids
    )

    # Create an assay with all properties
    assay = create_test_assay(
        client, "Multiple Property Assay", "Assay with multiple properties", assay_type["id"], property_ids
    )

    # Verify the properties were associated correctly
    assert len(assay["properties"]) == len(property_ids)

    # Check that all property IDs are in the response
    response_property_ids = [prop["id"] for prop in assay["properties"]]
    for prop_id in property_ids:
        assert prop_id in response_property_ids

    # Get the assay by ID to verify it persisted correctly
    get_response = client.get(f"/assays/{assay['id']}")
    assert get_response.status_code == 200
    get_data = get_response.json()

    # Verify the properties still match after retrieval
    assert len(get_data["properties"]) == len(property_ids)
    get_property_ids = [prop["id"] for prop in get_data["properties"]]
    for prop_id in property_ids:
        assert prop_id in get_property_ids


def test_create_assay_type_with_multiple_properties(client):
    """Test creating an assay type with multiple properties"""
    # Create properties
    property_ids = create_test_properties(client, test_properties)

    # Create an assay type
    assay_type = create_test_assay_type(
        client, "Multiple Property Assay Type", "For testing multiple properties", property_ids
    )

    # Verify the properties were associated correctly
    assert len(assay_type["properties"]) == len(property_ids)
    assay_type_property_ids = [prop["id"] for prop in assay_type["properties"]]
    for prop_id in property_ids:
        assert prop_id in assay_type_property_ids


def test_create_assay_with_assay_type(client):
    """Test creating an assay based on an assay type"""
    # Create properties
    property_ids = create_test_properties(client, test_properties)

    # Create an assay type
    assay_type = create_test_assay_type(
        client, "Template Assay Type", "For testing assay creation from template", property_ids
    )

    # Create an assay based on the assay type
    assay = create_test_assay(client, "Derived Assay", "Assay derived from template", assay_type["id"])

    # Verify the assay was created correctly
    assert assay["name"] == "Derived Assay"
    assert assay["description"] == "Assay derived from template"
    assert assay["assay_type_id"] == assay_type["id"]

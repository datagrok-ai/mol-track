import csv
import io
import pytest
import json
from app.utils import enums


def send_registration_request(client, rows):
    csv_buffer = io.StringIO()
    writer = csv.DictWriter(csv_buffer, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
    csv_buffer.seek(0)

    files = {"file": ("data.csv", csv_buffer.getvalue(), "text/csv")}
    data = {"error_handling": enums.ErrorHandlingOptions.reject_row.value}

    response = client.post("v1/compounds/", files=files, data=data)
    return response


validated_property_validator = {
    "properties": [
        {
            "name": "val_int_prop_1",
            "description": "val_int_prop_1 description",
            "value_type": "int",
            "property_class": "DECLARED",
            "unit": "",
            "entity_type": "COMPOUND",
            "pattern": "",
            "validators": '[">5", "<15"]',
        }
    ]
}


rows_validators = [
    {"smiles": "Brc1cc(Br)c2ccn(Br)c2c1", "val_int_prop_1": 10},
    {"smiles": "Cc1cc(Br)cc2c1ccn2Br", "val_int_prop_1": 16},
]


@pytest.mark.usefixtures("preload_simple_data")
def test_single_property_validators(client):
    schema_response = client.post("/v1/schema/", json=validated_property_validator)

    schema_response_content = schema_response.content.decode("utf-8")
    schema_response_content = json.loads(schema_response_content)

    assert schema_response_content["status"] == "success"
    assert schema_response_content["property_types"][0]["validators"] is not None

    response = send_registration_request(client, rows_validators)

    assert response.status_code == 200
    assert response.headers["content-disposition"].startswith("attachment;")

    content = response.content.decode("utf-8")
    content = json.loads(content)

    assert len(content) == 2
    assert content[0]["registration_status"] == "success"
    assert content[1]["registration_status"] == "failed"


validated_property_min_max = {
    "properties": [
        {
            "name": "val_int_prop_2",
            "description": "val_int_prop_2 description",
            "value_type": "int",
            "property_class": "DECLARED",
            "unit": "",
            "entity_type": "COMPOUND",
            "pattern": "",
            "min": 5,
            "max": 15,
        }
    ]
}


rows_min_max = [
    {"smiles": "Brc1cc(Br)c2ccn(Br)c2c1", "val_int_prop_2": 10},
    {"smiles": "Cc1cc(Br)cc2c1ccn2Br", "val_int_prop_2": 16},
]


@pytest.mark.usefixtures("preload_simple_data")
def test_single_property_min_max(client):
    schema_response = client.post("/v1/schema/", json=validated_property_min_max)

    schema_response_content = schema_response.content.decode("utf-8")
    schema_response_content = json.loads(schema_response_content)

    assert schema_response_content["status"] == "success"
    assert schema_response_content["property_types"][0]["min"] is not None
    assert schema_response_content["property_types"][0]["max"] is not None

    response = send_registration_request(client, rows_min_max)

    assert response.status_code == 200
    assert response.headers["content-disposition"].startswith("attachment;")

    content = response.content.decode("utf-8")
    content = json.loads(content)

    assert len(content) == 2
    assert content[0]["registration_status"] == "success"
    assert content[1]["registration_status"] == "failed"


validated_property_choices = {
    "properties": [
        {
            "name": "val_string_prop_1",
            "description": "val_string_prop_1 description",
            "value_type": "string",
            "property_class": "DECLARED",
            "unit": "",
            "entity_type": "COMPOUND",
            "pattern": "",
            "choices": '["choice1", "choice2", "choice3"]',
        }
    ]
}


rows_choices = [
    {"smiles": "Brc1cc(Br)c2ccn(Br)c2c1", "val_string_prop_1": "choice2"},
    {"smiles": "Cc1cc(Br)cc2c1ccn2Br", "val_string_prop_1": "choice4"},
]


@pytest.mark.usefixtures("preload_simple_data")
def test_single_property_choices(client):
    schema_response = client.post("/v1/schema/", json=validated_property_choices)

    schema_response_content = schema_response.content.decode("utf-8")
    schema_response_content = json.loads(schema_response_content)

    assert schema_response_content["status"] == "success"
    assert schema_response_content["property_types"][0]["choices"] is not None

    response = send_registration_request(client, rows_choices)

    assert response.status_code == 200
    assert response.headers["content-disposition"].startswith("attachment;")

    content = response.content.decode("utf-8")
    content = json.loads(content)

    assert len(content) == 2
    assert content[0]["registration_status"] == "success"
    assert content[1]["registration_status"] == "failed"


properties = {
    "properties": [
        {
            "name": "val_string_prop_1",
            "description": "val_string_prop_1 description",
            "value_type": "string",
            "property_class": "DECLARED",
            "unit": "",
            "entity_type": "COMPOUND",
            "pattern": "",
        },
        {
            "name": "val_int_prop_1",
            "description": "val_int_prop_1 description",
            "value_type": "int",
            "property_class": "DECLARED",
            "unit": "",
            "entity_type": "COMPOUND",
            "pattern": "",
        },
    ]
}


rows_record_validation = [
    {"smiles": "Brc1cc(Br)c2ccn(Br)c2c1", "val_string_prop_1": "input_1", "val_int_prop_1": 11},
    {"smiles": "Cc1cc(Br)cc2c1ccn2Br", "val_string_prop_1": "inputnotok", "val_int_prop_1": 16},
    {"smiles": "CC(=O)CC(=O)c1cc(Br)cc2c1ccn2Br", "val_string_prop_1": "input_ok", "val_int_prop_1": 2},
]


valid_rule = "matches(${val_string_prop_1}, r'^input_') && ${val_int_prop_1} > 10"


@pytest.mark.usefixtures("preload_simple_data")
def test_validator_registration_and_record_validation_valid(client):
    # Register properties
    schema_response = client.post("/v1/schema/", json=properties)

    schema_response_content = schema_response.content.decode("utf-8")
    schema_response_content = json.loads(schema_response_content)

    assert schema_response_content["status"] == "success"

    # Register validator
    validator_payload = {
        "name": "test_validator_1",
        "entity_type": enums.EntityType.COMPOUND.value,
        "expression": valid_rule,
        "description": "test validator 1 description",
    }
    register_rule_response = client.post("/v1/validators/", data=validator_payload)

    register_rule_response_content = register_rule_response.content.decode("utf-8")
    register_rule_response_content = json.loads(register_rule_response_content)

    assert register_rule_response_content["status"] == "success"
    assert register_rule_response_content["added_validator"] == valid_rule

    response = send_registration_request(client, rows_record_validation)

    assert response.status_code == 200
    assert response.headers["content-disposition"].startswith("attachment;")

    content = response.content.decode("utf-8")
    content = json.loads(content)

    assert len(content) == 3
    assert content[0]["registration_status"] == "success"
    assert content[1]["registration_status"] == "failed"
    assert content[2]["registration_status"] == "failed"


invalid_rule = "matches(${val_string_prop_1}, r'^input_') & ${val_int_prop_1} > 10"


@pytest.mark.usefixtures("preload_simple_data")
def test_validator_registration_invalid_property(client):
    # Register validator
    validator_payload = {
        "name": "test_validator_1",
        "entity_type": enums.EntityType.COMPOUND.value,
        "expression": invalid_rule,
        "description": "test validator 1 description",
    }
    register_rule_response = client.post("/v1/validators/", data=validator_payload)

    register_rule_response_content = register_rule_response.content.decode("utf-8")
    register_rule_response_content = json.loads(register_rule_response_content)

    assert register_rule_response_content["detail"]["status"] == "failed"
    assert register_rule_response_content["detail"]["message"].startswith("Error adding validators:")


@pytest.mark.usefixtures("preload_simple_data")
def test_validator_registration_invalid_rule(client):
    # Register properties
    schema_response = client.post("/v1/schema/", json=properties)

    schema_response_content = schema_response.content.decode("utf-8")
    schema_response_content = json.loads(schema_response_content)

    assert schema_response_content["status"] == "success"

    # Register validator
    validator_payload = {
        "name": "test_validator_1",
        "entity_type": enums.EntityType.COMPOUND.value,
        "expression": invalid_rule,
        "description": "test validator 1 description",
    }
    register_rule_response = client.post("/v1/validators/", data=validator_payload)

    register_rule_response_content = register_rule_response.content.decode("utf-8")
    register_rule_response_content = json.loads(register_rule_response_content)

    assert register_rule_response_content["detail"]["status"] == "failed"
    assert register_rule_response_content["detail"]["message"].startswith("Error adding validators:")

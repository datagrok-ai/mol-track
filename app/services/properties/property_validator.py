import json
from typing import Any
from app import models
from app.services.properties.numeric_constraint import NumericConstraint


class PropertyValidator:
    """
    Validates properties based on their type and associated constraints.
    """

    @classmethod
    def validate_value(cls, value: Any, property: models.Property) -> None:
        """
        Validate a string against a regex pattern and/or a list of choices.
        """

        if property.value_type in ("int", "double"):
            cls.check_numeric_constraints(value, property)

        if property.value_type == "string":
            cls.check_string_constraints(value, property)

        if property.validators:
            cls.check_validators(value, property)

    @classmethod
    def check_validators(cls, value: Any, property: models.Property) -> None:
        """
        Check if the value passes all validators defined for the property.
        """

        try:
            coerced_value = float(value) if property.value_type == "double" else int(value)
        except (ValueError, TypeError):
            raise ValueError(f"Value '{value}' must be a {property.value_type}")

        validators = json.loads(property.validators.replace("'", '"'))
        for validator in validators:
            constraint = NumericConstraint.parse(validator)
            if constraint:
                if not constraint.is_satisfied_for(coerced_value):
                    raise ValueError(f"Value '{coerced_value}' does not satisfy the validator: {validator}")

    @classmethod
    def check_numeric_constraints(cls, value: Any, property: models.Property) -> None:
        """
        Check if the numeric value is within the defined constraints.
        """

        try:
            coerced_value = float(value) if property.value_type == "double" else int(value)
        except (ValueError, TypeError):
            raise ValueError(f"Value '{value}' must be a {property.value_type}")

        if property.min is not None and coerced_value < property.min:
            raise ValueError(f"Value {coerced_value} is less than the minimum allowed {property.min}")

        if property.max is not None and coerced_value > property.max:
            raise ValueError(f"Value {coerced_value} is greater than the maximum allowed {property.max}")

    @classmethod
    def check_string_constraints(cls, value: str, property: models.Property) -> None:
        """
        Check if the string value matches the defined pattern.
        """

        if property.choices:
            choices = json.loads(property.choices.replace("'", '"'))
            if value not in choices:
                raise ValueError(f"Value '{value}' is not in the allowed choices: {choices}")

from fastapi import HTTPException
from typing import List, Optional
from pytest import Session
from sqlalchemy import insert

from app import models
from app.services.properties.complex_validator import ComplexValidator
from app.utils import enums
from app.utils.admin_utils import admin


def get_validators_for_entity(db: Session, entity: enums.EntityType) -> List[str]:
    results = db.query(models.Validator).filter(models.Validator.entity_type == entity).all()
    return results


def delete_validator_by_name(db: Session, name: str) -> None:
    validator = db.query(models.Validator).filter(models.Validator.name == name).first()
    if not validator:
        raise HTTPException(status_code=404, detail="Validator not found")
    db.delete(validator)
    db.commit()
    return validator


def create_validator(
    db: Session, entity: enums.EntityType, name: str, expression: str, description: Optional[str] = None
) -> models.Validator:
    if not expression:
        raise HTTPException(status_code=400, detail="Validator is not provided")

    properties = (
        db.query(models.Property.name, models.Property.value_type).filter(models.Property.entity_type == entity).all()
    )
    properties = {name: val_type.value for (name, val_type) in properties}

    try:
        ComplexValidator.validate_rule(expression, properties)
    except Exception as e:
        raise HTTPException(
            status_code=400, detail={"status": "failed", "message": f"Error adding validators: {str(e)}"}
        )

    try:
        insert_stmt = insert(models.Validator).values(
            [
                {
                    "name": name,
                    "entity_type": entity.value,
                    "expression": expression,
                    "description": description,
                    "created_by": admin.admin_user_id,
                    "updated_by": admin.admin_user_id,
                }
            ]
        )
        db.execute(insert_stmt)
        db.commit()
        return {"status": "success", "added_validator": expression}
    except Exception as e:
        db.rollback()
        raise HTTPException(
            status_code=500, detail={"status": "failed", "message": f"Error adding validators: {str(e)}"}
        )

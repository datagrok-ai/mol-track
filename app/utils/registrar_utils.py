from typing import Dict, Any, Type
from app import models
from app.utils import enums
from sqlalchemy.orm import Session


def get_validation_prefix(entity_type: enums.EntityType) -> str:
    validation_name_map = {
        enums.EntityType.COMPOUND: "compound_details",
        enums.EntityType.BATCH: "batch_details",
        enums.EntityType.ASSAY: "assay_details",
        enums.EntityType.ASSAY_RUN: "assay_run_details",
        enums.EntityType.ASSAY_RESULT: "assay_result_details",
    }
    return validation_name_map.get(entity_type, "")


def get_details_for_entity(
    entity_id: int,
    cache: Dict[int, Dict[str, Any]],
    entity_type: enums.EntityType,
    db: Session,
    detail_model: Type,
    foreign_key_field: str,
) -> Dict[str, Any]:
    if entity_id in cache:
        return cache[entity_id]

    properties = (
        db.query(detail_model, models.Property)
        .join(models.Property, detail_model.property_id == models.Property.id)
        .filter(getattr(detail_model, foreign_key_field) == entity_id)
        .all()
    )

    prop_dict: Dict[str, Any] = {}
    for detail, prop in properties:
        entity_scope = f"{get_validation_prefix(entity_type)}"
        prop_name = f"{prop.name}"

        scope_dict = prop_dict.setdefault(entity_scope, {})
        match prop.value_type:
            case enums.ValueType.string:
                scope_dict[prop_name] = detail.value_string
            case enums.ValueType.int | enums.ValueType.double:
                scope_dict[prop_name] = detail.value_num
            case enums.ValueType.datetime:
                scope_dict[prop_name] = detail.value_datetime
            case enums.ValueType.uuid:
                scope_dict[prop_name] = detail.value_uuid

    cache[entity_id] = prop_dict
    return prop_dict

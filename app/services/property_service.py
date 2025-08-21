from fastapi import HTTPException
from sqlalchemy import inspect
from app.utils import type_casting_utils, enums
from app.utils.admin_utils import admin
from typing import Callable, Dict, Any, List, Optional, Tuple, Type


class PropertyService:
    def __init__(self, property_records_map: Dict[str, Any]):
        self.property_records_map = property_records_map

    def get_property_info(self, prop_name: str, entity_type: enums.EntityType) -> Dict[str, Any]:
        props = self.property_records_map.get(prop_name)
        if not props:
            raise HTTPException(status_code=400, detail=f"Unknown property: {prop_name}")

        matching_prop = next((p for p in props if getattr(p, "entity_type", None) == entity_type), None)

        if not matching_prop:
            existing_types = [getattr(p, "entity_type", None) for p in props]
            raise HTTPException(
                status_code=400,
                detail=f"Property '{prop_name}' has entity_type(s) {existing_types}, expected '{entity_type.value}'",
            )

        value_type = getattr(matching_prop, "value_type", None)
        if (
            value_type not in type_casting_utils.value_type_to_field
            or value_type not in type_casting_utils.value_type_cast_map
        ):
            raise HTTPException(status_code=400, detail=f"Unsupported or unknown value type for property: {prop_name}")

        return {
            "property": matching_prop,
            "value_type": value_type,
            "field_name": type_casting_utils.value_type_to_field[value_type],
            "cast_fn": type_casting_utils.value_type_cast_map[value_type],
        }

    # TODO: Design a more robust and efficient approach for handling updates to compound details
    def build_details_records(
        self,
        model: Type,
        properties: Dict[str, Any],
        entity_ids: Dict[str, Any],
        entity_type: enums.EntityType,
        include_user_fields: bool = True,
        update_checker: Optional[Callable[[str, int, Any], Optional[Dict[str, Any]]]] = None,
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        records_to_insert, records_to_update = [], []

        for prop_name, value in properties.items():
            prop_info = self.get_property_info(prop_name, entity_type)
            prop = prop_info["property"]
            value_type = prop_info["value_type"]
            cast_fn = prop_info["cast_fn"]
            field_name = prop_info["field_name"]
            prop_id = getattr(prop, "id")

            if value in ("", "none", None):
                continue

            #  Detect and parse value qualifiers
            value_qualifier = 0
            if value_type in ("double", "int") and isinstance(value, str):
                qualifiers = {
                    "<": 1,
                    ">": 2,
                    "=": 0,
                }

                first_char = value[0] if value else ""
                if first_char in qualifiers:
                    value_qualifier = qualifiers[first_char]
                    value = value[1:].strip()

            try:
                casted_value = cast_fn(value)
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Error casting value for property {prop_name}: {e}")

            detail = {
                **entity_ids,
                "property_id": prop_id,
                field_name: casted_value,
            }

            # TODO: Refactor to generically handle all value_* fields without hardcoding model-specific attributes
            mapper = inspect(model)
            value_columns = [
                col.key for col in mapper.columns if col.key.startswith("value") and col.key not in field_name
            ]

            for col_name in value_columns:
                if col_name == "value_qualifier":
                    detail[col_name] = value_qualifier
                    continue

                column = mapper.columns[col_name]
                default = None
                if column.default is not None:
                    if not callable(column.default.arg):
                        default = column.default.arg
                detail[col_name] = default

            if include_user_fields:
                detail["created_by"] = admin.admin_user_id
                detail["updated_by"] = admin.admin_user_id

            if update_checker:
                result = update_checker(entity_ids, detail, field_name, casted_value)
                if result.action == "skip":
                    continue
                elif result.action == "update":
                    records_to_update.append(result.update_data)
                    continue
                elif result.action == "insert":
                    pass

            records_to_insert.append(detail)

        return records_to_insert, records_to_update

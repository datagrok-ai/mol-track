from fastapi import HTTPException
from sqlalchemy import inspect
from utils import type_casting_utils
import main
from typing import Callable, Dict, Any, List, Optional, Tuple, Type


class PropertyService:
    def __init__(self, property_records_map: Dict[str, Any]):
        self.property_records_map = property_records_map

    def get_property_info(self, prop_name: str, scope: str) -> Dict[str, Any]:
        prop = self.property_records_map.get(prop_name)
        if prop is None:
            raise HTTPException(status_code=400, detail=f"Unknown property: {prop_name}")

        retieved_scope = getattr(prop, "scope", None)
        if retieved_scope != scope:
            raise HTTPException(
                status_code=400, detail=f"Property '{prop_name}' has scope '{retieved_scope}', expected '{scope}'"
            )

        value_type = getattr(prop, "value_type", None)
        if (
            value_type not in type_casting_utils.value_type_to_field
            or value_type not in type_casting_utils.value_type_cast_map
        ):
            raise HTTPException(status_code=400, detail=f"Unsupported or unknown value type for property: {prop_name}")

        return {
            "property": prop,
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
        scope: str,
        include_user_fields: bool = True,
        update_checker: Optional[Callable[[str, int, Any], Optional[Dict[str, Any]]]] = None,
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        records_to_insert, records_to_update = [], []

        for prop_name, value in properties.items():
            prop_info = self.get_property_info(prop_name, scope)
            prop = prop_info["property"]
            cast_fn = prop_info["cast_fn"]
            field_name = prop_info["field_name"]
            prop_id = getattr(prop, "id")

            try:
                casted_value = None if value in ("", None) else cast_fn(value)
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Error casting value for property {prop_name}: {e}")

            detail = {
                **entity_ids,
                "property_id": prop_id,
                field_name: casted_value,
            }

            # TODO: Refactor to generically handle all value_* fields without hardcoding model-specific attributes
            mapper = inspect(model)
            value_columns = [col.key for col in mapper.columns if col.key.startswith("value") and col.key != field_name]

            for col_name in value_columns:
                column = mapper.columns[col_name]
                default = None
                if column.default is not None:
                    if not callable(column.default.arg):
                        default = column.default.arg
                detail[col_name] = default

            if include_user_fields:
                detail["created_by"] = main.admin_user_id
                detail["updated_by"] = main.admin_user_id

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

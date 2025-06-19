from fastapi import HTTPException
from datetime import datetime
import utils
import main
from typing import Dict, Any, List


class PropertyService:
    def __init__(self, property_records_map: Dict[str, Any]):
        self.property_records_map = property_records_map

    def get_property_info(self, prop_name: str) -> Dict[str, Any]:
        prop = self.property_records_map.get(prop_name)
        if prop is None:
            raise HTTPException(status_code=400, detail=f"Unknown property: {prop_name}")

        value_type = getattr(prop, "value_type", None)
        if value_type not in utils.value_type_to_field or value_type not in utils.value_type_cast_map:
            raise HTTPException(status_code=400, detail=f"Unsupported or unknown value type for property: {prop_name}")

        return {
            "property": prop,
            "value_type": value_type,
            "field_name": utils.value_type_to_field[value_type],
            "cast_fn": utils.value_type_cast_map[value_type],
        }

    def build_details_records(
        self, properties: Dict[str, Any], entity_ids: Dict[str, Any], include_user_fields: bool = True
    ) -> List[Dict[str, Any]]:
        records = []
        for prop_name, value in properties.items():
            prop_info = self.get_property_info(prop_name)
            detail = {
                **entity_ids,
                "property_id": getattr(prop_info["property"], "id"),
                "value_datetime": datetime.now(),
                "value_num": 0,
                "value_string": None,
            }
            if include_user_fields:
                detail["created_by"] = main.admin_user_id
                detail["updated_by"] = main.admin_user_id

            try:
                detail[prop_info["field_name"]] = prop_info["cast_fn"](value)
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Error casting value for property {prop_name}: {e}")

            records.append(detail)
        return records

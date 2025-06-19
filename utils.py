from datetime import datetime
from typing import Any, Callable, Dict
import uuid
from dateutil.parser import parse


def cast_datetime(value: Any) -> datetime:
    if isinstance(value, datetime):
        return value
    try:
        return datetime.fromisoformat(value)
    except Exception:
        try:
            # Fall back to dateutil for flexible parsing
            return parse(value)
        except Exception:
            raise ValueError(f"Invalid datetime format: {value}")


def cast_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        val = value.strip().lower()
        if val in ("true", "1"):
            return True
        elif val in ("false", "0"):
            return False
    if isinstance(value, (int, float)):
        return bool(value)
    raise ValueError(f"Invalid boolean format: {value}")


def cast_uuid(value: Any) -> uuid.UUID:
    if isinstance(value, uuid.UUID):
        return value
    try:
        return uuid.UUID(str(value))
    except Exception:
        raise ValueError(f"Invalid UUID format: {value}")


value_type_to_field: Dict[str, str] = {
    "datetime": "value_datetime",
    "int": "value_num",
    "double": "value_num",
    "string": "value_string",
    "uuid": "value_uuid",
    "bool": "value_bool",
}

value_type_cast_map: Dict[str, Callable[[Any], Any]] = {
    "datetime": cast_datetime,
    "int": int,
    "double": float,
    "string": str,
    "uuid": cast_uuid,
    "bool": cast_bool,
}

import csv
from datetime import datetime
from io import BytesIO, StringIO
import json
import re
from typing import Any, Dict, List
import uuid
from fastapi import Response
import pandas as pd
from sqlalchemy.orm import Session
from sqlalchemy import text
from app.models import Level, Aggregation
from app.utils import enums


def create_alias_mapping(columns: List[str], aggregations: List[Aggregation]) -> Dict[str, str]:
    output_aliases = {sanitize_field_name(field): (field, (field, None)) for field in columns}
    output_aliases.update(
        {
            sanitize_field_name(agg.field, agg.operation): (
                f"{agg.operation.value}({agg.field})",
                (agg.field, agg.operation.value),
            )
            for agg in aggregations
        }
    )
    return output_aliases


def create_alias(table: Level) -> str:
    if table != enums.SearchEntityType.ASSAY_RUNS.value:
        name_parts = table.split("_")
        if len(name_parts) == 2:
            # Since the table name is in the format "table_name", we can use the first letter of the first part and the
            # first letter of the second part
            return f"{name_parts[0][0]}{name_parts[1][0]}"
        else:
            return f"{name_parts[0][0]}"
    else:
        return "rn"


def singularize(word: str) -> str:
    if word.endswith("es"):
        return word[:-2]
    elif word.endswith("s"):
        return word[:-1]
    return word


def get_table_columns(table_name: str, session: Session) -> list:
    """
    Get the columns of a table in the database.
    """

    result = session.execute(
        text("SELECT column_name FROM information_schema.columns WHERE table_name = :table_name"),
        {"table_name": table_name},
    )
    columns = [row[0] for row in result.fetchall()]
    return columns


def sanitize_field_name(field_name: str, agg_op: str = None) -> str:
    """Sanitize field name for use in SQL aliases"""
    # Replace dots and special characters with underscores
    sanitized = re.sub(r"[^a-zA-Z0-9_]", "_", field_name)
    # Ensure it doesn't start with a number
    if sanitized and sanitized[0].isdigit():
        sanitized = f"field_{sanitized}"
    if agg_op:
        agg_op = agg_op.replace(" ", "_")
        sanitized = f"{agg_op.lower()}_{sanitized}"
        sanitized = sanitized.lower()
    return sanitized


def convert_datetime(obj):
    if isinstance(obj, datetime):
        return obj.isoformat()
    return obj


def prepare_search_output(results: List[Any], headers: List[str], output_format: enums.SearchOutputFormat):
    match output_format:
        # UUID objects are not JSON serializable, so convert to str to prevent circular reference errors
        case enums.SearchOutputFormat.json:
            return_obj = {
                "status": "success",
                "total_count": len(results),
                "columns": headers,
                "data": [
                    dict(zip(headers, (str(cell) if isinstance(cell, uuid.UUID) else cell for cell in row)))
                    for row in results
                ],
            }
            json_output = json.dumps(return_obj, default=convert_datetime)
            return Response(
                content=json_output,
                media_type="application/json",
                headers={"Content-Disposition": "attachment; filename=result.json"},
            )
        case enums.SearchOutputFormat.csv:
            output = StringIO()
            writer = csv.writer(output)
            writer.writerow(headers)
            writer.writerows(results)
            csv_content = output.getvalue()
            return Response(
                content=csv_content,
                media_type="text/csv",
                headers={"Content-Disposition": "attachment; filename=result.csv"},
            )
        case enums.SearchOutputFormat.parquet:
            df = pd.DataFrame(results, columns=headers)
            buffer = BytesIO()
            df.to_parquet(buffer, engine="pyarrow", index=False)
            return Response(
                content=buffer.getvalue(),
                media_type="application/octet-stream",
                headers={"Content-Disposition": "attachment; filename=result.parquet"},
            )

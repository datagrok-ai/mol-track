import re
from typing import List
from sqlalchemy.orm import Session
from sqlalchemy import text


def create_alias(table: str) -> str:
    name_parts = table.split("_")
    if len(name_parts) == 2:
        # Since the table name is in the format "table_name", we can use the first letter of the first part and the
        # first letter of the second part
        return f"{name_parts[0][0]}{name_parts[1][0]}"
    else:
        return f"{name_parts[0][0]}"


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


def sanitize_field_name(field_name: str) -> str:
    """Sanitize field name for use in SQL aliases"""
    # Replace dots and special characters with underscores
    sanitized = re.sub(r"[^a-zA-Z0-9_]", "_", field_name)
    # Ensure it doesn't start with a number
    if sanitized and sanitized[0].isdigit():
        sanitized = f"field_{sanitized}"
    return sanitized


class JoinOrderingTool:
    """
    Utility class to manage and organize SQL JOIN clauses for complex search queries.

    This tool helps in constructing and ordering JOIN statements for different
    entities (such as compounds, batches, details, assay results, and properties)
    when building dynamic SQL queries. It maintains a dictionary of join lists,
    categorized by entity type, and provides methods to expand and update these
    lists as needed.
    """

    def __init__(self):
        self.joins_dict = {
            "compounds": set(),
            "batches": set(),
            "assay_results": set(),
            "compound_details": set(),
            "batch_details": set(),
            "assay_result_details": set(),
            "properties": set(),
        }

    def add(self, joins: List[str], keys: List[str]) -> bool:
        for i, join in enumerate(joins):
            self.joins_dict[keys[i]].add(join)

    def getListOfJoins(self) -> List[str]:
        keys = [
            "compounds",
            "batches",
            "compound_details",
            "batch_details",
            "assay_results",
            "assay_result_details",
            "properties",
        ]
        joins = [join for key in keys for join in self.joins_dict.get(key, set())]
        return " ".join(joins) if joins else ""

    def getJoinSQL(self, reversed: bool = False) -> str:
        if not reversed:
            keys = ["compounds", "batches", "assay_results"]
        else:
            keys = ["assay_results", "batches", "compounds"]
        keys.extend(["compound_details", "batch_details", "assay_result_details", "properties"])
        joins = [join for key in keys for join in self.joins_dict.get(key, set())]
        return " ".join(joins) if joins else ""

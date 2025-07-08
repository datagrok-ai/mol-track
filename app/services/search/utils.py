import re
from typing import List


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
            "compound_details": set(),
            "batch_details": set(),
            "assay_results": set(),
            "properties": set(),
        }

    def add(self, joins: List[str], keys: List[str]) -> bool:
        for i, join in enumerate(joins):
            self.joins_dict[keys[i]].add(join)

    def getListOfJoins(self) -> List[str]:
        joins = []
        joins.extend(list(self.joins_dict.get("compounds", set())))
        joins.extend(list(self.joins_dict.get("batches", set())))
        joins.extend(list(self.joins_dict.get("compound_details", set())))
        joins.extend(list(self.joins_dict.get("batch_details", set())))
        joins.extend(list(self.joins_dict.get("assay_results", set())))
        joins.extend(list(self.joins_dict.get("properties", set())))

        return joins

    def getJoinSQL(self) -> str:
        joins = []
        joins.extend(list(self.joins_dict.get("compounds", set())))
        joins.extend(list(self.joins_dict.get("batches", set())))
        joins.extend(list(self.joins_dict.get("compound_details", set())))
        joins.extend(list(self.joins_dict.get("batch_details", set())))
        joins.extend(list(self.joins_dict.get("assay_results", set())))
        joins.extend(list(self.joins_dict.get("properties", set())))

        return " ".join(joins) if joins else ""

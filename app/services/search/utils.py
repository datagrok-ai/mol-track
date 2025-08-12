import re
from typing import List
from sqlalchemy.orm import Session
from sqlalchemy import text
from collections import deque
from app.models import Level


def create_alias(table: Level) -> str:
    if table != "assay_runs":
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
    when building dynamic SQL queries.
    """

    def __init__(self):
        self.keys = []
        self.joins = []

    def add(self, joins: List[str], keys: List[str]) -> bool:
        for i in range(len(joins)):
            if keys[i] not in self.keys:
                self.keys.append(keys[i])
                self.joins.append(joins[i])

    def getJoinSQL(self) -> str:
        return " ".join(self.joins) if self.joins else ""

    def getLastTableAlias(self) -> str:
        if len(self.keys):
            return create_alias(self.keys[-1])

    def checkLastJoin(self, table: Level) -> bool:
        if len(self.keys):
            return table == self.keys[-1]
        return False

    def joinCount(self):
        return len(self.joins)


class JoinResolutionError(Exception):
    """Custom exception for join resolution errors"""

    pass


class JoinResolver:
    def __init__(self, schema):
        self.relationships = {
            ("compounds", "batches"): f"INNER JOIN {schema}.batches b ON b.compound_id = c.id",
            ("batches", "compounds"): f"INNER JOIN {schema}.compounds c ON c.id = b.compound_id",
            ("batches", "assay_results"): f"INNER JOIN {schema}.assay_results ar ON ar.batch_id = b.id",
            ("assay_results", "batches"): f"INNER JOIN {schema}.batches b ON b.id = ar.batch_id",
            ("assay_results", "assay_runs"): f"INNER JOIN {schema}.assay_runs rn ON rn.id = ar.assay_run_id",
            ("assay_runs", "assay_results"): f"INNER JOIN {schema}.assay_results ar ON rn.id = ar.assay_run_id",
            ("assays", "assay_runs"): f"INNER JOIN {schema}.assay_runs rn ON rn.assay_id = a.id",
            ("assay_runs", "assays"): f"INNER JOIN {schema}.assays a ON a.id = rn.assay_id",
        }
        self.graph = self._build_graph()

    def _build_graph(self):
        graph = {}
        for (t1, t2), join in self.relationships.items():
            graph.setdefault(t1, []).append((t2, join))
        return graph

    def _find_path(self, start, end):
        queue = deque([(start, [])])
        visited = set()

        while queue:
            current, path = queue.popleft()
            if current == end:
                return path
            visited.add(current)
            for neighbor, join_clause in self.graph.get(current, []):
                if neighbor not in visited:
                    queue.append((neighbor, path + [(neighbor, join_clause)]))
        return None

    def resolve_join_components(self, from_level: str, to_level, subquery: bool = False, details: bool = False):
        result = self._find_path(from_level, to_level)
        if not result:
            raise JoinResolutionError(f"No relationship defined from {from_level} to {to_level}")
        from_table = None
        if subquery:
            from_table = from_level
            if len(result) >= 1 and result[0][1].find(f"{singularize(from_level)}_id") != -1:
                from_table = result.pop(0)[0]
            if details and len(result) >= 1 and result[-1][1].find(f"{singularize(to_level)}_id") != -1:
                result.pop(-1)
        joins, tables = [], []
        for table, join in result:
            joins.append(join)
            tables.append(table)
        return joins, tables, from_table

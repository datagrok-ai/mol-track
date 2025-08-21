from collections import deque
from typing import List
from app.models import Level
from app.services.search.utils.helper_functions import create_alias, singularize


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
            if keys[i] not in self.keys or (keys[i] == "properties" and joins[i] not in self.joins):
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
    def __init__(self, schema, table_configs):
        self._generate_relationships(table_configs, schema)
        self.graph = self._build_graph()

    def _generate_relationships(self, table_configs, schema):
        self.relationships = {}
        for from_name, from_table in table_configs.items():
            for to_name, to_table in table_configs.items():
                if from_name == to_name:
                    continue
                from_fk, from_alias = from_table["details_fk"], from_table["alias"]
                to_fk, to_alias = to_table["details_fk"], to_table["alias"]
                if from_fk in to_table["direct_fields"]:
                    self.relationships[(from_name, to_name)] = (
                        f"INNER JOIN {schema}.{to_name} {to_alias} ON {to_alias}.{from_fk} = {from_alias}.id"
                    )
                if to_fk in from_table["direct_fields"]:
                    self.relationships[(from_name, to_name)] = (
                        f"INNER JOIN {schema}.{to_name} {to_alias} ON {to_alias}.id = {from_alias}.{to_fk}"
                    )

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

    def resolve_join_components(
        self, from_level: Level, to_level: Level, subquery: bool = False, details: bool = False
    ):
        join_path = self._find_path(from_level, to_level)
        if not join_path:
            raise JoinResolutionError(f"No relationship defined from {from_level} to {to_level}")
        from_table = None
        if subquery:
            from_table = from_level
            if len(join_path) >= 1 and join_path[0][1].find(f"{singularize(from_level)}_id") != -1:
                from_table = join_path.pop(0)[0]
            if details and len(join_path) >= 1 and join_path[-1][1].find(f"{singularize(to_level)}_id") != -1:
                join_path.pop(-1)
        joins, tables = [], []
        for table, join in join_path:
            joins.append(join)
            tables.append(table)
        return joins, tables, from_table

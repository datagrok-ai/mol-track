import re
from typing import List, Optional, Union


class SimpleValidator:
    # Operators
    NONE = "none"
    EQUALS = "="
    NOT_EQUALS = "!="
    GT = ">"
    GTOE = ">="
    LT = "<"
    LTOE = "<="
    IN = "in"
    NOT_IN = "not in"
    RANGE = "-"
    IS_NULL = "is null"
    IS_NOT_NULL = "is not null"

    # Regex patterns
    _num = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
    numRegex = re.compile(f"^{_num}")
    unaryRegex = re.compile(r"(?:[$]{.*})?\s*(!=|>=|<=|=|>|<|\bnot in\b|\bin\b)")
    rangeRegex = re.compile(rf"({_num})\s?(?:-|\.{2})\s?({_num})")
    inRegex = re.compile(r"\((.*)\)")
    isNullRegex = re.compile(r"(^is null$|^is not null$)")

    unaryOperators = [EQUALS, NOT_EQUALS, GT, GTOE, LT, LTOE, IN, NOT_IN]

    def __init__(
        self, op: str, v1: Optional[float] = None, v2: Optional[float] = None, values: Optional[List[float]] = None
    ):
        self.op = op
        self.v1 = v1
        self.v2 = v2
        self.values = values or []

    @staticmethod
    def parse(query: str) -> Optional["SimpleValidator"]:
        query = query.strip()
        if not query:
            return SimpleValidator(SimpleValidator.NONE)

        # Simple numeric → equals
        try:
            x = float(query)
            return SimpleValidator(SimpleValidator.EQUALS, x)
        except ValueError:
            pass

        # Unary operators (=, !=, >, <, etc.)
        match = SimpleValidator.unaryRegex.search(query)
        if match:
            op = match.group(1)
            rest = query[match.end() :].strip()

            if op in (SimpleValidator.IN, SimpleValidator.NOT_IN):
                match_values = SimpleValidator.inRegex.search(rest)
                if match_values:
                    values = [float(s.strip()) for s in match_values.group(1).split(",")]
                    return SimpleValidator(op, values=values)
            else:
                match_value = SimpleValidator.numRegex.match(rest)
                if match_value:
                    v1 = float(match_value.group(0))
                    return SimpleValidator(op, v1=v1)

        # Range like 5-10
        range_match = SimpleValidator.rangeRegex.search(query)
        if range_match:
            return SimpleValidator(SimpleValidator.RANGE, float(range_match.group(1)), float(range_match.group(2)))

        # Null check
        null_match = SimpleValidator.isNullRegex.match(query.lower())
        if null_match:
            return SimpleValidator(null_match.group(1))

        return None

    def match(self, x: Optional[Union[int, float]]) -> bool:
        if x is None:
            return self.op in (self.NONE, self.IS_NULL)

        if self.op == self.NONE:
            return False
        if self.op == self.GT:
            return x > self.v1
        if self.op == self.GTOE:
            return x >= self.v1
        if self.op == self.LT:
            return x < self.v1
        if self.op == self.LTOE:
            return x <= self.v1
        if self.op == self.EQUALS:
            return x == self.v1
        if self.op == self.NOT_EQUALS:
            return x != self.v1
        if self.op == self.IN:
            return x in self.values
        if self.op == self.NOT_IN:
            return x not in self.values
        if self.op == self.RANGE:
            return self.v1 <= x <= self.v2
        if self.op == self.IS_NULL:
            return False
        if self.op == self.IS_NOT_NULL:
            return True
        raise ValueError(f"Unknown operation {self.op}")

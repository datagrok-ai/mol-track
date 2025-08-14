from typing import List, Dict, Any
import re
from sqlalchemy.orm import Session
from app.models import Token
from app.utils.enums import CompareOp, LogicOp
from app.services.search.parser import Parser
from app.services.search.engine import SearchEngine


class SearchFilterBuilder:
    def __init__(self, db: Session):
        self.COMPARE_OPS = [op.value for op in CompareOp]
        self.COMPARE_OPS = sorted(self.COMPARE_OPS, key=lambda s: len(s.split()), reverse=True)
        self.LOGICAL_OPS = [op.value for op in LogicOp]
        self.db = db

    def build_filter(self, expression: str) -> Dict[str, Any]:
        expression = expression.replace('"', "'")
        tokens = self._tokenize(expression=expression)

        parser = Parser(tokens, self.LOGICAL_OPS, self.COMPARE_OPS)
        filter_tree = parser.parse()

        search_engine = SearchEngine(self.db)
        err = search_engine._validate_filter(filter_tree, "")
        return {"filter": filter_tree.model_dump() if err == [] else None, "errors": err}

    def _tokenize(self, expression: str) -> List[Token]:
        tokens: List[Token] = []
        i = 0
        length = len(expression)

        while i < length:
            c = expression[i]

            # Skip whitespace
            if c.isspace():
                i += 1
                continue

            # Parentheses
            if c in "()":
                tokens.append(Token("LPAREN" if c == "(" else "RPAREN", c))
                i += 1
                continue

            # Array literal
            if c == "[":
                end = i
                bracket_level = 1
                while end + 1 < length and bracket_level > 0:
                    end += 1
                    if expression[end] == "[":
                        bracket_level += 1
                    elif expression[end] == "]":
                        bracket_level -= 1
                if bracket_level != 0:
                    raise ValueError("Unclosed array literal")

                array_content = expression[i + 1 : end].strip()
                if not array_content:
                    values = []
                else:
                    raw_items = re.split(r",\s*", array_content)
                    values = []
                    for item in raw_items:
                        item = item.strip()
                        if item.startswith("'") and item.endswith("'"):
                            values.append(item[1:-1])
                        elif re.match(r"^-?\d+\.\d+$", item):
                            values.append(float(item))
                        elif re.match(r"^-?\d+$", item):
                            values.append(int(item))
                        elif item.lower() == "true":
                            values.append(True)
                        elif item.lower() == "false":
                            values.append(False)
                        else:
                            raise ValueError(f"Unsupported array item: {item}")

                tokens.append(Token("ARRAY_LITERAL", values))
                i = end + 1
                continue

            # String literal (single quotes)
            if c == "'":
                end = i + 1
                end = expression.find("'", i + 1)
                if end == -1:
                    raise ValueError("Unclosed string literal")
                tokens.append(Token("STRING_LITERAL", expression[i + 1 : end]))
                i = end + 1
                continue

            # Number literal
            number_match = re.match(r"-?\d+(\.\d+)?", expression[i:])
            if number_match:
                num_str = number_match.group(0)
                value = float(num_str) if "." in num_str else int(num_str)
                tokens.append(Token("NUMBER_LITERAL", value))
                i += len(num_str)
                continue

            # Boolean literal
            match = re.match(r"(true|false)", expression[i:], re.IGNORECASE)
            if match:
                tokens.append(Token("BOOLEAN_LITERAL", match.group(1).lower() == "true"))
                i += len(match.group(0))
                continue

            # Multi-word operators
            matched_op = next((op for op in self.COMPARE_OPS + self.LOGICAL_OPS if expression[i:].startswith(op)), None)
            if matched_op:
                token_type = "COMPARE_OP" if matched_op in self.COMPARE_OPS else "LOGICAL_OP"
                tokens.append(Token(token_type, matched_op))
                i += len(matched_op)
                continue

            # Field name (dot notation, alphanumeric + underscores)
            field_match = re.match(r"[a-zA-Z_][a-zA-Z0-9_.]*", expression[i:])
            if field_match:
                field = field_match.group(0)
                tokens.append(Token("FIELD", field))
                i += len(field)
                continue

            raise ValueError(f"Unexpected character at position {i}: '{c}'")
        return tokens

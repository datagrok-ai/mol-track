from typing import List, Union
from pydantic import ValidationError
from app.models import AtomicCondition, LogicalNode, Token, Filter
from app.utils.enums import CompareOp, LogicOp


class ParseError(Exception):
    """Raised on invalid syntax or validation issues during parsing."""

    pass


class Parser:
    def __init__(self, tokens: List[Token], LOGICAL_OPS: List[str], COMPARE_OPS: List[str]):
        self.tokens = tokens
        self.pos = 0
        self.ops = {
            "LOGICAL_OP": LOGICAL_OPS,
            "COMPARE_OP": COMPARE_OPS,
        }

    def peek(self) -> Union[Token, None]:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def advance(self) -> Token:
        tok = self.peek()
        self.pos += 1
        return tok

    def expect(self, expected_type: str) -> Token:
        tok = self.peek()
        if not tok or tok.type != expected_type:
            raise ParseError(
                f"Expected token of type '{expected_type}', but got {tok}. Possible values: {self.ops[expected_type]}"
            )
        return self.advance()

    def parse(self) -> Filter:
        result = self.parse_expression()
        if self.peek() is not None:
            raise ParseError(f"Unexpected trailing token: {self.peek()}")
        return result

    def parse_expression(self) -> Filter:
        left = self.parse_term()

        while self.peek() and self.peek().type == "LOGICAL_OP":
            op_token = self.advance()
            logic_op_str = op_token.value

            try:
                logic_op = LogicOp(logic_op_str)
            except ValueError:
                raise ParseError(f"Invalid logical operator: {logic_op_str}")

            right = self.parse_term()

            if isinstance(left, LogicalNode) and left.operator == logic_op:
                left.conditions.append(right)
            else:
                left = LogicalNode(operator=logic_op, conditions=[left, right])

        return left

    def parse_term(self) -> Filter:
        tok = self.peek()

        if tok and tok.type == "LPAREN":
            self.advance()
            node = self.parse_expression()
            if not self.peek() or self.peek().type != "RPAREN":
                raise ParseError("Expected closing parenthesis ')'")
            self.advance()
            return node
        else:
            return self.parse_atomic()

    def parse_atomic(self) -> AtomicCondition:
        field_token = self.expect("FIELD")
        operator_token = self.expect("COMPARE_OP")
        value_token = self.advance()

        if not value_token or value_token.type not in {
            "STRING_LITERAL",
            "NUMBER_LITERAL",
            "BOOLEAN_LITERAL",
            "ARRAY_LITERAL",
        }:
            raise ParseError(f"Expected a value literal, got {value_token}")

        field = field_token.value
        operator_str = operator_token.value

        try:
            operator = CompareOp(operator_str)
        except ValueError:
            raise ParseError(f"Unsupported comparison operator: {operator_str}")

        value = value_token.value
        threshold = None

        # Optional threshold for IS SIMILAR
        if operator == CompareOp.IS_SIMILAR:
            if self.peek() and self.peek().type == "NUMBER_LITERAL":
                threshold_token = self.advance()
                threshold = float(threshold_token.value)
            else:
                raise ParseError("IS SIMILAR operator requires a numeric threshold after the value")

        try:
            return AtomicCondition(field=field, operator=operator, value=value, threshold=threshold)
        except ValidationError as e:
            raise ParseError(f"Validation failed: {e}")

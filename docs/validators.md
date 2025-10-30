# Data Validation System Documentation

## Overview

This system provides a flexible and multi-layered validation framework for database entities.
Validation happens at two levels:

1. **Single-property validation** – field-level constraints such as numeric ranges, allowed string values, and specific numeric conditions.
2. **Complex validation** – record-level rules that involve multiple fields, defined using the **CEL (Common Expression Language)**.

The goal is to ensure data integrity, enforce business rules, and provide extensibility for future validation logic.

---

## 1. Single-Property Validation

Every property (database field) can have **validation metadata** stored alongside it. The metadata fields are:

- `min` (numeric, optional) - Lower bound for numeric values
- `max` (numeric, optional) - Upper bound for numeric values
- `choices` (string, optional) - A JSON-encoded list of valid string values
- `validators` (list, optional) - A JSON-encoded list of advanced numeric validation expressions

If a field’s metadata entry is populated, the value must satisfy **all conditions** associated with that field.

### 1.1 Min/Max Constraints

For numeric properties:

- `min` - Value must be greater than or equal to this number
- `max` - Value must be less than or equal to this number

**Example:**

```json
  min =  0
  max = 100
```

**Validation:**
The value must lie between **0** and **100** (inclusive).

### 1.2 Choices (for Strings)

For string properties, the system can restrict allowed values using a JSON list of choices.

**Example:**

```json
  choices = '["draft", "submitted", "approved", "rejected"]'
```

**Validation:**
The value must be one of `"draft"`, `"submitted"`, `"approved"`, or `"rejected"`.
### 1.3 Validators (Advanced Numeric Conditions)

Validators are stored as a JSON-encoded list of numeric expressions, parsed using the `NumericConstraint` class.

#### Supported Operators

| Operator      | Example        | Meaning                                   |
|---------------|----------------|-------------------------------------------|
| `=`           | `= 10`         | Value must equal 10                       |
| `!=`          | `!= 5`         | Value must not equal 5                    |
| `>`           | `> 3`          | Value must be greater than 3              |
| `>=`          | `>= 7`         | Value must be greater than or equal to 7  |
| `<`           | `< 20`         | Value must be less than 20                |
| `<=`          | `<= 100`       | Value must be less than or equal to 100   |
| `Range (-)`   | `5-10`         | Value must be between 5 and 10 (inclusive)|
| `Range (..)`  | `5..10`        | Same as above                             |
| `in (...)`    | `in (1,2,3)`   | Value must be one of 1, 2, or 3           |
| `not in (...)`| `not in (4,5)` | Value must not be one of 4 or 5           |
| `is null`     | `is null`      | Field must be NULL                        |
| `is not null` | `is not null`  | Field must not be NULL                    |

#### Example JSON

```json
  validators = ["> 10", "<= 50", "!= 25", "in (5, 15, 20, 30)"]
```

### 1.4 Single-Property Validation Flow

1. Check if the field has validation metadata.
2. If **min** or **max** exist - enforce numeric range.
3. If **choices** exist - check if string value is in the list.
4. If **validators** exist - parse each validator expression and test against the value.
5. Value passes only if **all conditions** are satisfied.

## 2. Complex Validation (Record-Level)

Single-property validation is not always sufficient. Some rules span multiple fields
(e.g., `end_date` must be greater than `start_date`).

For this, the system supports **CEL-based rules** with custom functions and preprocessing.

### 2.1 CEL (Common Expression Language) with Extended Functionality

Rules are written in CEL expressions and executed against the entire record (dictionary of property values).

CEL supports:

- Arithmetic operations
- Comparison and logical operators
- String and collection operations
- Date/time operations
- Control flow (ternary expressions)
- The above custom functions

#### Custom Functions Available in CEL Context

- `size(value)` - Returns the length of a string or list.
- `matches(value, pattern)` - Returns `true` if the string `value` matches the regex `pattern`.
- `today()` - Returns the current date as an ISO string.
- `date(str)` - Converts a string in ISO format to a date object.

#### Preprocessing Features

- `${field}` - Translated to `field`.
- `${field}.length` - Translated to `size(field)`.
- `is null` - Translated to `== null`.
- `is not null` - Translated to `!= null`.

### 2.2 Example Rules

**Numeric Relationship**

```cel
age >= 18 && age <= 120
```

***Meaning:*** Age must be between 18 and 120.

**Cross-Field Validation**

```cel
end_date > start_date
```

***Meaning***: The end date must be later than the start date.

**String Conditions**

```cel
status in ["draft", "submitted", "approved"]
```

***Meaning***: Status must be one of the allowed states.

For more information on expression format, check out the [CEL library documentation](https://python-common-expression-language.readthedocs.io/en/latest/tutorials/cel-language-basics/).

### 2.3 Validation Flow

1. Collect record values as a map `{ field: value }`.
2. Load all CEL rules applicable to the entity.
3. Evaluate each rule against the record.
4. Validation passes if **all rules return true**.
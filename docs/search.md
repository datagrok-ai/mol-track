# Search
The search functionality allows users to query compounds, batches, and assay results using flexible filters and field selection. Users can define which fields to return, apply precise conditions using logical operators, and build simple or complex queries. This functionality and be accessed through emdpoints ```/v1/search/compounds```, ```/v1/search/batches``` and ```/v1/search/assay_results```. This document explains how to structure search requests and use the available options effectively. 
# Search request
Format of the search filter is as follows:
```json
{
  "level": <level>,
  "output": <output_list>,
  "filter": <filter>
}
```
* `<level>` - Specifies the main entity to search over.

  Accepted values:

  * *compounds*

  * *batches*

  * *assay_results*
  
  This field is automatically populated based on the called endpoind. 

* `<output>` - A list of wanted output fields (all fields in this list must be related to the specified `<level>` and must follow the format rules described in [next section](#fields))
* `<filter>` - Filter represents the filter criteria and it's exact format is described in [this section](#filter).
## Fields
Field names follow a specific notation, depending on whether a standard field or a dynamic property is referenced:
* `<table>.<field>`

  Refers to a direct column of the table.

  *Examples*: `compounds.canonical_smiles`, `batches.batch_regno`

* `<table>.details.<propertyName>`

  Refers to a dynamic property stored in the corresponding details table:

  * *compound_details* for compounds

  * *batch_details* for batches

  * *assay_results* for assays - WIP

  *Example*: `compounds.details.MolLogP`

This notation applies when writing filter expressions and stating output fields.
## Filter
The filter defines the search criteria based on fields, values, and logical operators.

Filter supports:
* A simple condition (e.g., a field equals a value)
* A complex condition created by combining multiple conditions using *AND* / *OR* logic
### Supported operators
The search functionality supports the following operators based on the data type:
* String: `=`, `!=`, `IN`, `STARTS WITH`, `ENDS WITH`, `LIKE`, `CONTAINS`
* Number: 
  * `=`, `!=`, `<`, `>`, `<=`, `>=`
  * `RANGE (min, max)` 
* Datetime: `BEFORE`, `AFTER`, `ON`
* Bool: `=`, `!=`
* Molecule (uses RDKit cartridge operator)
  * `HAS SUBSTRUCTURE (value: string)`
  * `IS SUBSTRUCTURE OF (value: string)` 
  * `IS SIMILAR (value: string, threshold: number)`

**Note:** Molecule operations can be applied only to field ```compounds.structure```
### Simple Condition Format
  Each simple condition in the filter object follows a common structure:
  ```json
  {
    "field": <field_name>,
    "operator": <operator>,
    "value": <value>,
    "threshold": <number or null>
  }
```
Field descriptions:
- ***field*** - the name of the field to filter on (as describes in [this section](#fields))
- ***operator*** - the comparison operator to use (as described in [this section](#supported-operators))
- ***value*** - the target value or list of values to compare against; can be a string, number, or array
- ***threshold*** - only used with the `IS SIMILAR` operator; for all other operators, this should be left out or `null`
### Complex Condition
Complex condition follows a common structure:
```json
{
  "operator": <"AND" or "OR">,
  "conditions": <list_of_conditions>
}
```
Field descriptions:
* *operator* - Logical operator "AND" or "OR"
* *conditions* - List of conditions (list can contain simple conditions, complex conditions or both)
### Example Expression
```json
{
   "level":"compounds",
   "output":[
      "compounds.canonical_smiles",
      "compounds.details.chembl",
      "compounds.details.polarSurface"
   ],
   "filter":{
      "operator":"AND",
      "conditions":[
         {
            "field":"compounds.structure",
            "operator":"IS SIMILAR",
            "value":"Cc1ccc",
            "threshold":0.95
         },
         {
            "operator":"OR",
            "conditions":[
               {
                  "field":"compounds.details.chembl",
                  "operator":"IN",
                  "value":[
                     "CHEMBL123",
                     "CHEMBL123"
                  ]
               },
               {
                  "field":"compounds.details.project",
                  "operator":"=",
                  "value":"My project"
               }
            ]
         },
         {
            "operator":"AND",
            "conditions":[
               {
                  "field":"assay_results.ic50",
                  "operator":">",
                  "value":0.42
               },
               {
                  "field":"assay_results.ec50",
                  "operator":"<",
                  "value":100
               }
            ]
         }
      ]
   }
}
```
## Search query parser
### Overview
Users can use the `/vi/search/generate-filter` endpoint to create structured filters from a human-readable expression. 
This makes building complex search filters easier and more intuitive.
### Expression Format
Users can use:
- **Fields**: Any [searchable field](#fields) (like `compounds.details.chembl`)
- **Operators**: Any [search operator](#supported-operators) (like `IN` or `IS SIMILAR`)
- **Values**: strings (`'CHEMBL1'`), numbers (`3.5`), or arrays (`['A', 'B']`)
- **Logical operators**: `AND`, `OR`
- **Parentheses** To ensure conditions are evaluated in the intended order.
### Example Expressions
#### 1. Simple AND
```
(compounds.details.chembl IN ['CHEMBL1', 'CHEMBL2'] AND compounds.details.MolLogP > 1.5)
```
#### 2. Similarity
```
compounds.canonical_smiles IS SIMILAR 'CCCO' 0.9
```
**Note**: Operator `IS SIMILAR` requires two values smiles and threshold.

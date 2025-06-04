# Search

We want to be able to search on any level (compounds | batches | assay_results), with dynamically constructed filtering 
criteria on any level.

When defining filters for the level lower than the query level, it means that the condition has to be satisfied for
any of the child entities, not for all. So if we search on a compounds level for "assay_results.ic50 > 0.5", we 
would get all compounds for which we have at least one recorded ic50 > 0.5

Within one condition group, you can't mix multiple levels. This can be checked in advance before executing the query.

Fields notation (applies to compounds.xxx, batches.xxx, assay.results):
* <tablename>.property defines immediate table field, such as compounds.canonical_smiles, or batches.description
* <tablename>.details.<propertyName> defines a dynamic property linked via the corresponding details table (compound_details,  batch_details, or assay_results)

## Supported criteria

Most of the criteria are written in a form <field> <operator> <value>.

* String: `=`, `!=`, `IN`, `STARTS WITH`, `ENDS WITH`, `LIKE`, `CONTAINS`
* Number: 
  * `=`, `!=`, `<`, `>`, `<=`, `>=`
  * `RANGE (min, max)` 
* Datetime: `BEFORE`, `AFTER`, `ON`
* Bool: `=`, `!=`
* Molecule (uses RDKit cartridge operator)
  * `CONTAINS` (value: string)
  * `IS CONTAINED` (value: string)
  * `IS SIMILAR (value: string, threshold: number)`


```json
{
  "level": "compounds",                // compounds | batches | assay_results
  "output": [
    "compounds.canonical_smiles",      // columns to return
  	"compounds.details.chembl",
	  "compounds.details.polarSurface"
  ],
  "filter": {
    "operator": "AND",
    "conditions": [
      {
        "field": "compounds.structure",
        "operator": "is similar",
        "value": "Cc1ccc",
        "threshold": 0.95
      },
      
      {
        "operator": "OR",
        "conditions": [
          {
            "field": "compounds.details.chembl",
            "operator": "in",
            "value": ["CHEMBL123", "CHEMBL123"]
          },
          {
            "field": "compounds.details.project",
            "operator": "=",
            "value": "My project"
          }
        ]
      },

      {
        "operator": "AND",
        "conditions": [
          {
            "field": "assay_results.ic50",
            "operator": ">",
            "value": 0.42
          },
          {
            "field": "assay_results.ec50",
            "operator": "<",
            "value": 100
          }
        ]
      }
    ]
  }
}
```
# Example 2 Dataset #

This is a synthetic dataset with a reasonably coherent set of molecule and assay data at both single point and curve level.

## Compounds ##

- [schema](compounds_schema.json)
- [mapping](compounds_mapping.json)
- [data](compounds.csv)

## Batches ##

- [schema](batches_schema.json)
- [mapping](batches_mapping.json)
- [data](batches.csv)

## Assay Data ##

- Assay Types
  - [schema](assay_types_schema.json)
- Assays
  - [schema](assays_schema.json)
  - [mapping](assays_mapping.json)
  - [data](assays.csv)
- Assay Results
  - [mapping](.)
  - [data](.)

```json
{
    "Assay Date": "assays.assay_date",
    "Assay Name": "assay_types.name",
    "Assay Type": "assay_type -> details -> (cell-based, biochemical) -> (assay_format, something else)",
    "End Point Unit (nM)": "assay_type_properties missing some details in unit. this is an IC50",
    "Experiment ID": "assays -> details -> experiment id", 
    "Hill Coefficient": "assay_type_properties -> Hill coefficient",
    "IC50 Absolute": "",
    "Inflection": "",
    "Max": "curve max",
    "Min": "curve min",
    "Protein Batch": "assays -> details",
    "Protocol ID": "assays -> details - this is combo of assay_type and target",
    "Target": "assays -> details - target is experimental condition of an assay run",
    "Qualified Reported IC50", "",
    "rSquared": "",
    "Signal to Background": "",
    "Therapeutic Area": "",
    "Z Prime": "",
    "Cell Line": "",
    "Number results": "",
    "All results": "",
    "Avg IC50": "",
    "StdDev": ""
}
```
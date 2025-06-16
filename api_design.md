# API Design Thinking #

## **WORK-IN-PROGRESS** ##

We assemble some thinking about the general structure of APIs for [Moltrack](https://github.com/datagrok-ai/mol-track)

What are the strategies?  Typical RESTful API design where there is a close connection between the underlying database design and the endpoint design.  For each database table, there are corresponding CRUD needs fulfilled by CRUD endpoints. (CRUD - create, read, update, delete).  The corresponding REST method are POST, GET, PUT or PATCH, and DELETE.  There is an open standard, [json:api](https://jsonapi.org/), that Revvity's API is modeled after, that we may want to consider.

We have a number of intents.

1. set up system
2. set up properties, synonyms, associations
3. register virtual compounds
4. register batches
5. set up assay_types and corresponding assay_type_details and assay_type_properties
6. register assays and assay_results.
7. search within a domain (compounds, batches, assay_results) or in combination

**IMPORTANT: All paths should be prefixed by the appropriate version.**  Right now the version would be `/v1`.  In this document, the version prefix is assumed for clarity.  FASTAPI has a [router](https://fastapi.tiangolo.com/reference/apirouter/) concept that may allow us to declare versioning.

## Schema - WIP ##

When a user wants to register data into [Moltrack](https://github.com/datagrok-ai/mol-track), they typically want to register main entities like batches, compounds, or assay data.  The data set will typically have suppementary data like properties and synonyms that add context.  Moltrack was designed to handle this data in a controlled way.  To do that we need to allow the user to register the definitions of that context data.  We call this definition a schema.  The contract of this schema definition is to ensure that the appropriate definitions exist, to reuse previously defined definitions where possible, and to reject definitions that are incomplete.  The 'schema' is used in the definition of a set of resources but **not** a resource in and of itself.

*We will be developing this by example to understand its breath and depth.*  We know that we need to be able to define basic entities like properties and synonyms and associations.  The immediate need is for creations of these definitions.  Later we will think about the design for updates and (logical) deletions.

- Semantic Types: name, description
- Properties: name, value_type, semantic_type, property_class, unit, scope
- Synonym Types: synonym_level, name, pattern, description
- Additions: name, description, code, is_active, formula, molecular_weight, smiles, molfile, role

- `POST /schema` - used to define/ensure that a set of context definitions exist.

   For the [Black demo data set](./demo-data/black)  to register compounds and batches we need the following:

   ```json
   {
       "properties": [
         {"name": "MolLogP", "scope": "compound", "property_class": "calculated", "value_type": "double", "unit": ""}
       ],
       "synonym_types": [
            {"name": "batch_corporate_id", "level": "batch", "pattern": ""},
            {"name": "corporate_id", "level": "compound", "pattern": ""},
            {"name": "CAS number", "level": "compound", "pattern": "r'\b[1-9]{1}[0-9]{1,6}-\d{2}-\d\b'"},
            {"name": "common name", "level": "compound", "pattern": ""},
            {"name": "USAN", "level": "compound", "pattern": ""}
            ]
   }
   ```

  - Output
  The output should indicate the success of the overall request (`200 OK` or `201 Created`). *I can't identify the appropriate response for mixed suceess.* There should also be an item-by-item success status and if necessary a failure message:
    - Success - new property or synonym type is successfully registered.
    - Skipped - property or synonym type exists in MolTrack already.
    - Failed - the insert of the property or synonym type failed.  A 'failure message' attribute should be added to the item.

### Getter endpoints for schema ###

The GET endpoints allow us to see the configuration of MolTrack from perspective of the properties, synonyms, etc. associated with entities in the system.  At the moment we are provide resources specific to the entities but could also make it more dynamic with a query parameter.  We could also consider adding equivalent paths like `/compounds/schema` and `/batches/schema`

- `GET /schema/compounds` - returns properties and synonym_types associated with the compounds entity_type using the json form identified above.
- `GET /schema/batches` - returns properties, synonym_types and additions associated with the batches entity_type using the json form identified above.

## Additions ##

Additions are addtional chemical entities present in the batch/lot of a material produced.  Examples include HCl salts or hydrates.  For these to be available as a part of the batch registration they need to be pre-registered into MolTrack.  One of smiles and molfile must exist.  Any missing values should be computed if possible.  The return will be an array of additions added including the addition_id. `Is_archived` is archived by default; `Is_active` is true by default.  The MolTrack administator may update an addition's `is_active` flag (via `PUT`) to `false` to disallow subsequent use of the addtion.  A work-in-progress set of additions can be found in [additions.csv](./demo-data/additions.csv).

- `POST /additions` - used to register one or more additions into the MolTrack system.  The endpoint accepts an array of addition definitions in csv format (selected as the format for the MVP).  The return will show the status of the overal request. There should also be an item-by-item success status and if necessary a failure message:
  - Success - new addition is successfully registered.
  - Skipped - addition exists in MolTrack already.
  - Failed - the insert of the addition failed.  A 'failure message' attribute should be added to the item.

- `PUT /additions` – updates one or more existing additions in the MolTrack system. Accepts an array of complete addition definitions, each including a valid addition_id. Returns an array of the updated additions. The return will show the status of the overal request. There should also be an item-by-item success status and if necessary a failure message:
  - Success - new addition is successfully updated.
  - Skipped - no addition update required; addition exists and all attributes match in MolTrack already.
  - Failed - the addition update failed.  A 'failure message' attribute should be added to the item.

- `GET /additions` - retrieve all the salts and solvates that have been registered.  Retrieve in a csv format.
- `GET /additions/salts` - retrieve all additions with role of *salts*.
- `GET /additions/solvates` - retrieve all additions with role of *solvates*.
- `GET /additions/{addition_id}` - retrieve all information for a specific addition.

### Update endpoints ###

- `PUT /additions/{addition_id}` – update information for the specified addition_id.
- `DELETE /additions/{addition_id}` - soft delete from a given addition but only if there are no existing dependent batches.

## Register Batches ##

Batch registration will be performed as singletons or in bulk, by the chemist or registrar.  Use cases include supporting individual batch synthesis, library synthesis generated internally or by a contract research organization (CRO) and acquistion from a vendor.  Each registrations entry will contain information about the batch including batch details, batch synonyms, compound information (including synonyms and properties) and batch additions.  Properties, synonym_types, and additons must be defined ahead of time and if missing will result in a registration failure.  The inputs will be a collection of batch records that will be presented in a CSV-style format, an SDfile format or a parquet file .  A single registration will be a collection with a single record.  

Registrations could executed synchronously or asynchronously.  In general, registration should be performed asynchronously to accommmodate multiple scenarios where processing would be require a lot of time, including molecules with complex ring structure or a collection that contains many entries.  *For the MVP, we have decided to go with synchronous registration.*

Mapping is optional, but assumes that the field names directly map to the database concepts.

- `POST /batches` - synchronous registration endpoint for CSV input
  - Input
    - query parameters
      - ~~`input_format=csv`  [default csv, optionally mol-v3000]~~
      - `output_format=csv` [default csv, optional mol-v3000]
    - example in csv

      ```text
      "batch_corporate_id", structure, common_name, cas, usan
      "EPA-001-001","OC(=O)c1cc2c(cc1)c(=O)O","1,3-benezenedicarboxylic acid","121-91-5"
      "EPA-120-001","c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N",,"169590-42-5","celecoxib"
      ```

    - csv data in the json body

      ```json
      {
         "data": "\"batch_corporate_id\", smiles, common_name, cas, usan\n\"EPA-001-001\",\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",     \"121-91-5\"\n\"EPA-120-001\",\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\"\n",
         "mapping": {
            "smiles": "compounds.smiles",
            "batch_corporate_id": "batches_synonyms.batch_corporate_id",
            "common_name": "compounds_synonyms.common_name",
            "cas": "compounds_synonyms.cas",
            "usan": "compounds_synonyms.usan"
         }
      }
      ```

    - An example set of batches with addtions and batch properties to be registered can be found [2_batches_with_additions.csv](./demo-data/black/2_batches_with_additions.csv).  The additions are listed in columns with the cells showing the equivalents.  Schemas posts should be performed with [compounds_schema.json](./demo-data/black/compounds_schema.json) and [batches_schema.json](./demo-data/black/batches_schema.json) first to prepare the appropriate schema.

  - Output
   The batches, batches_detail, batches_synonym, batch_additions, additions, compounds, compound_details, and compound_synonyms tables will be joined appropriately, pivoted and concatenated appropriately to produce a single record.  The batch_regno and id should also be returned.  The row-by-row registration status and if necessary error message should be added as attributes.
    - 200

      ```json
      {
         "status_message": "Success",
         "data": "\"batch_corporate_id\",smiles,common_name,cas,usan,batch_id,compound_id,compound_corporate_id\n\"EPA-001-001\",\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",\"121-91-5\",1,2,\"EPA-001\"\n\"EPA-120-001\",\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\",2,4,\"EPA-120\"\n"
      }
      ```

- `GET /batches/` - returns a representation of batches including details, synonyms, additions and compounds including details and synonyms

   __Future potential capabilities__
   1. pagination must be supported
   2. output format must be supported. `output-format=[json|csv|mol-v3000]` with a default of `csv`.  If `json` format is selected, the data is returned in a nested dictionary structure.  If the format is `csv` or `mol-v3000`, then the data is pivoted/concatenated/flattened
   3. there should probably be a hard limit

- `GET /batches/{id}` - Returns a representation of batches including details, synonyms, additions and compounds including details and synonyms

   __Future potential capability__
  - Query parameter `output-format=[json|csv|mol-v3000]` with a default of `csv`.  If `json` format is selected, the data is returned in a nested dictionary structure.  If the format is `csv` or `mol-v3000`, then the data is pivoted/concatenated/flattened

- `GET /batches/{batch_id}/properties`
- `GET /batches/{batch_id}/synonyms`
- `GET /batches/{batch_id}/additions` - should returns additions, equivalents, and summary

- `PUT /batches/{batch_id}` - Used for updating information (details, compound, synonym, additions) for a given batch.  `batch_regno` is not updateable.

## Register Virtual Compounds ##

Virtual compound registration will typically be done in bulk, by a cheminformatician or registrar.  These registrations will be contain the structure, synonyms, and possibly properties.  Properties and synonym_types will need to be defined ahead of time.  The inputs will be a collection of compound records that will be presented in a CSV-style format, an SDfile format or a parquet file. A single registration will be a array with a single record. 

Registrations could be executed synchronously or asynchronously.  In general, registration should be performed asynchronously to accommmodate multiple scenarios where processing would be require a lot of time, including molecules with complex ring structure or a collection that contains many entries.  

*For the MVP, we have decided to go with*:

   1. *Synchronous registration.*
   2. *csv as the only input format with the structure encoded as either smiles or molfile (V2000 or V3000) format.*
   3. *csv as the only output format.*
   4. *Two failure-handling strategies are possible: (a) Reject all – the entire request fails if any entry is invalid; (b) Reject row – invalid entries are skipped, while valid ones are processed. By default, the system uses the Reject all strategy.*

- `POST /compounds` - synchronous registration endpoint for CSV input
  - Input
    - *Deferred* query parameters
      - `input_format=csv`  [default csv, optionally mol-v3000] - *No longer required.  The data format should be determined by the content.*
      - `output_format=csv` [default csv, optional mol-v3000]
    - example in csv

      ```text
      "structure, common_name, cas, usan
      "OC(=O)c1cc2c(cc1)c(=O)O","1,3-benezenedicarboxylic acid","121-91-5"
      "c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N",,"169590-42-5","celecoxib"
      ```

    - Mapping data is optional included in the input body. If not provided, it will be inferred from column names.  If provided, it will follow the format below:

      ```json
      {
         "mapping": {
            "<column name>": "<field name>",
            "<column_name>": "<table_name>.[<field_name | dynamic attribute name>]"
         }
      }
      ```

    - Example csv data in the json body with mapping

      ```json
      {
         "data": "structure, common_name, cas, usan\n\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",     \"121-91-5\"\n\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\"\n",
         "mapping": {
            "structure": "smiles",
            "common_name": "compounds_synonyms.common_name",
            "cas": "compounds_synonyms.cas",
            "usan": "compounds_synonyms.usan"
         }
      }
      ```

  - Output
   The compounds, compound_details, and compound_synonyms tables will be joined appropriately, pivoted and concatenated appropriately to produce a single record.
    - 200

      ```json
      {
         "status_messge": "Success",
         "data": "smiles,common_name,cas,usan,batch_id,compound_id,compound_corporate_id,registration_status,registration_error_message\n\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",\"121-91-5\",2,\"EPA-001\","success",""\n\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\",4,\"EPA-120\","success",""\n"
      }
      ```

    - 422 Unprocessable Content  - *For the MVP, any process error will for entire transaction to fail* 

      ```json
      {
         "status_messge": "Failed",
         "data": "smiles,common_name,cas,usan,batch_id,compound_id,compound_corporate_id,registration_status,registration_error_message\n\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",\"121-91-5\",2,\"EPA-001\","success",""\n\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\",4,\"EPA-120\","failed","invalid structure"\n"
      }
      ```

### Getter endpoints ###

These are included for symmetry of user experience.  They may be deprecated in favor of a search experience for become synonyms for that search endpoint.

- `GET /compounds/` - returns an array of compounds with each entry having the comound, the properties and the synonyms.  An optional query parameter will allow the user to select the output format: json, csv-style, sd-file, parquet.  Properties and synonyms will need to be pivoted (and concatenated as necessary) to accomodate the csv, sd and parquet format expectations.  ~~We will need to have query parameters to support pagination of the resulting output (`start`,`stop`,`max-per-page`)~~
- `GET /compounds/{compound_id}`
- `GET /compounds/{compound_id}/properties` - I think these will be infrequently used but included for completeness
- `GET /compounds/{compound_id}/synonyms` - I think these will be infrequently used but included for completeness

### Update endpoints ###

- `PUT /compounds/{compound_id}` - Used to update information (structure, properties, synonyms) for the provided compound_id.  Based on the business rule configuration, structure changes for compounds with batches attached are not permitted but rather may only be performed via the `PUT /batches/{batch_id}` endpoint.
- `DELETE /compounds/{compound_id}` - Used to perform a soft delete.  Only allowed from compounds that have no dependent batches.

## Assay Data Domain ##

A key capability for moltrack is to capture assay data related to a sample (batch).  Typically this will be measured biological activity (potency, selectivity, toxicity).  It can be used to capture measured physical/chemical attributes as well.

- `POST /schema/assay`  will define allowed/expected properties for assay_type_details, assay_details, assay_type_properties.  The properties represent categorization of the assay type whose values would be stored in *assay_type_details*, experimental conditions that would be stored as the assays level in the *assay_details* table, and result types and experimental conditions that would be stored at the *assay_results* level.  *in vivo*, *in vitro*, *in celluo* are examples of an **assay format** property that would likely be declared at the *assay type* level.

N.B. I am not happy with the titling of the major sections of the following schema

   ```json
   {
      "assay_type_details properties": [
         {
            "name": "assay format",
            "scope": "assay_types",
            "property_class": "asserted",
            "value_type": "string"
         },
         {
            "name": "biological system",
            "scope": "assay_types",
            "property_class": "asserted",
            "value_type": "string"
         }
      ],
      "assay_details properties":[
            {
            "name": "Cell Species",
            "value_type": "string",
            "required": true,
            "scope": "assay",
            "property_class": "asserted"
        },
        {
            "name": "Cell Lot",
            "value_type": "string",
            "required": true,
            "scope": "assay",
            "property_class": "measured"
        },
        {
            "name": "Cell Concentration",
            "value_type": "double",
            "required": true,
            "scope": "assay",
            "property_class": "measured",
            "units": "10^6 cells/mL"
        },
        {
            "name": "Assay Run Date",
            "value_type": "datetime",
            "required": true,
            "scope": "assay",
            "property_class": "measured"
        },
        {
            "name": "Assayer",
            "value_type": "uuid",
            "required": true,
            "scope": "assay",
            "property_class": "measured"
        }
      ],
       "assay_type_properties": [
            {
                "name": "Reported CLint",
                "value_type": "double",
                "unit": "uL/min/10^6 cells",
                "required": true,
                "property_class": "measured",
                "scope": "assay_results"
            },
            {
                "name": "Mean HTC recovery",
                "unit": "%",
                "required": false,
               "value_type": "double",
               "property_class": "measured",
                "scope": "assay_results"
            },
            {
                "name": "SD HTC recovery",
                "unit": "%",
                "required": false,
               "value_type": "double",
               "property_class": "measured",
                "scope": "assay_results"
            },
            {
                "name": "Dosed Concentration",
                "unit": "uM",
                "required": true,
               "value_type": "double",
               "property_class": "measured",
                "scope": "assay_results"
            }
       ]
   }
   ```

- `POST /assay_types` will create an instance of an assay type and values in assay_type_details and assay_type_properties. An example is presented below.

   ```json
   {
    "assay_type": {
        "name": "Hepatocyte Stability",
        "details": {
            "assay format": "in cellulo",
            "biological_system": "hepatocyte"
        }
      }
   }
   ```

- `POST /assays` will create an instance of the 'assays' or assay run and populate property values in *assay_details* table.  Example properties to be populated would include assay date, assayer, and might include cell lot number for a *in cellulo* assay.  An example is presented below:

   ```json
   {
      "assay_details": {
         "cell species": "Human",
         "cell lot": "H1",
         "assayer": "Jane Doe",
         "assay_date": "2024-02-01"
      }
   }
   ```

- `POST /assay_results` This endpoint will be used populate data in assays, assay_details, and assay_results with input from a csv file and mapping.  The properties that are populated here will mostly be result types like IC50, SD, % inhibtion, ...  Certain result level experimental conditions may also be populated here, like dosed concentration for a stability study.  The mostly like input will be a csv file with a row per sample (batch) and columns per property from assays and assay_results levels.  There will need to be a mapping.  Since there were be multiple results rows per assay run, the assays instance will be to be determined by matching appropriate properties.
  - See example [assay data](./demo-data/black/assay_results.csv)
  - See example [mapping](./demo-data/black/assay_results_mapping.json)

### Assay Data Getter endpoint ###

These are included for symmetrical thinking of user experience.  They may be deprecated in favor of a search experience or become alternative paths for that search endpoint.

- `GET /assay_types` - Return a list of all assay_types with included assay_type_details and assay_type_properties.
- `GET /assay_types/{assay_type_id}` - Return the details for a specific assay_type.
- `GET /assays` - Return a list of all assays including assay_type information, assay_details information.
- `GET /assays/{assay_id}` - Return the details for a specific assay.
- `GET /assay_results`
   Query parameters include possible filters for {assay_type_id} and/or {assay_id}

### Assay Data Putter endpoints ###

### Assay Data Deleter endpoints ###

An example CSV can be found here [assay_results.csv](./demo-data/black/assay_results.csv).  An example mapping can be seen here [assay_results_mapping.json](./demo-data/black/assay_results_mapping.json).  Note that all fields are not mapped.

## Search - WIP ##

We provide a generalized search capabilities across all the various entities in our model: compounds, batches, assay_results.  A simple or complex set of search criteria are presented in a json format.  ~~Query parameters provide pagination support and output format style.~~ For the MVP, pagination is not required and only the CSV style output will be provided.  Because our primary output is csv style output, the output content needs to be coherent and hence limited by the entity type row orientation.  If we were using json output we could imagine batches and compounds in one result set.

User stories:

   1. Find all batches made for project X in the past week.
   2. Find all stereo-similar compounds to my molecule, so I can determine whether I have expressed the stereo nature of my new molecule correctly.
   3. Find all batches made by WuXi for project Y.
   4. Show me all details for batch with a specific corporate id.
   5. Find all batches made for project Y including details and assay results
   6. find all assay-data where
      1. substructure criteria
      2. similarity criteria
      3. batch property of source = WuXi
      4. kinase IC50 < 100 uM
   7. Show me all the *in cellulo* assay_types with details
   8. Shom me all the assay runs for a specific assay_type

The endpoint path indicates the shape of the returned data: compound per row; batch per row; assay result value per row.

- `/search/compounds` - compound per row, batch data aggregated or pivoted to the compound row.
- `/search/batches` - batch per row, compound data denormalized to each batch row, assay results aggreated or pivoted to the batch row
- `/search/assay-data` -

- `POST /search`
- `POST /search/compounds/structure`
  - exact  -->  smiles + standardized by Moltrack settings.  optional pattern of standardization  [array of standardization steps].  uses all layers of the registration mol hash
    - returns - 0 or 1 compound entries, exact match might return multiple tautomers
  - less precise  -- substructure, tautomer, no-stereo, similar to a key compound
    - returns multiple rows

   Caller does not need to how we call the fields and parameters

    - exact( query_molecule: smiles)
    - substructure( query_molecule: smarts)
    - tautomer( query_molecule: smiles)
    - no-stereo( query_molecule: smiles )
    - similar( query_molecule: smiles, similarity_threshold: float = 0.9)

   How to structure the json?

   ```json
    {
      "query_mol": "smiles | smarts",
      "structure_search": "exact (default) | substructure | tautomer | no-stereo | similar "
    },
    { "structure_search_type": ""}
    ```
   If we can decide on the proper structure search condition expression then these can just be part of the `/search/compounds` body and we don't need the `/search/compounds/structure` path

- `POST /search/complex` - should this be the generalized search?  Krzysztof has implement body like

```json
{
  "conditions": [
    {
      "table": "compounds",
      "field": "hash_tautomer",
      "operator": "=",
      "query_smiles": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
      "columns": ["id", "canonical_smiles"]
    },
    {
      "table": "assays",
      "field": "IC50",
      "operator": "<",
      "unit": "nM"
      "value": 100,
      "columns": ["assay_id", "IC50", "target"]
    }
  ],
  "logic": "AND"
}
```

- `POST /search/batches`
- `POST /search/assay-data`

logical conditions -> postgres sql condition phrases

- operators
  - AND
  - OR
  - NOT
  - in
  - "="
  - "is similar"
  - ">"
  - "<"
  - range ?
  - exists
  - like
  - contains

- conditions
  - field
  - operator
  - value - required for all operators except exists
  - threshold - required only for the "is similar" operator

- level - is currently one of
  - compounds
  - batches
  - assay_results

  We may need
  - assay_types
  - assays

- output
  - "*" - means all fields

- join - Joins may be reequired when including output fields from different levels

____________________________________________________________________________________________________________________
## *DEFERRED* ##

Properties and synonym types can probably follow the standard approaches:

## *Deferred* - Properties ##

- `POST /properties` takes a json structure in the body of the call with fields corresponding to a single property definition, registers the single property and returns the json structure for that property and including the id.
- `GET /properties` returns an array of json structures detailing each property.  A query string, defining the scope and limiting the return set, is also allowed.
- `GET /properties/{property_id}` returns the json structure for a single property.
- do we need updates.
- register a set of properties at one time through a "schema"

## *Deferred* - Synonym Types ##

- `POST /synonym-types/`
- `GET /synonym-types/`
- `GET /synonym-types/{synonym_id}`

- *Deferred* `POST /batches-collection`  To be implemented post-MVP.

   A mapping between input fields and database structure will need to be provided or inferred.  The POST endpoint creates a job that needs to be performed and returns the `jobId`.
   For this endpoint we will need to consider the ability to include all the data in the json body, gzipping the content, or provide a file upload mechanism.  The backend job should probably separate the content into individual entries so that we can run in parallel threads and still capture progress and errors.

   Returns:

   202

   ```json
   {
      "status_path": "/batches-collection/123/status"
   }
   ```

- *Deferred* `GET /batches-collection/{jobId}/status`. To be implemented post-MVP.
   Returns

   200

   ```json
   {
      "batches-collection-job-id": 123456,
      "global-status": "In Progress",
      "individual-statuses": [
         { "batch1": "Success"},
         { "batch2": "Failed", "error-message": "Business rules failed: project not recognized" },
         { "batch3": "In Progress"}
      ]
   }
   ```

- *Deferred* `POST /compounds-collection/`
A simplified json example is presented here to demonstrate the input structure:

   ```json
   [
      {
         "smiles": "c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N",
         "USAN": "celecoxib",
         "CAS Number": "169590-42-5"
      }
   ]
   ```

   or

   ```json
   [
      {
         "smiles": "c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N",
         "synonyms": [
            {"USAN": "celecoxib"},
            {"CAS Number": "169590-42-5"}
         ],
         "properties": []
      }
   ]
   ```

   A mapping between input fields and database structure will need to be provided or inferred.  The POST endpoint creates a job that needs to be performed and returns the `jobId`.
   For this endpoint we will need to consider the ability to include all the data in the json body, gzipping the content, or provide a file upload mechanism.  The backend job should probably separate the content into individual entries so that we can run in parallel threads and still capture progress and errors.

   Returns:

   202

   ```json
   {
      "status_path": "/compounds-collection/123/status"
   }
   ```

- *Deferred* `GET /compounds-collection/{jobId}/status` returns

   200

   ```json
   {
      "compounds-collection-job-id": 123456,
      "global-status": "In Progress",
      "individual-statuses": [
         { "molecule1": "Success"},
         { "molecule2": "Failed", "error-message": "Business rules failed: chemist not recognized" },
         { "molecule3": "In Progress"}
      ]
   }
   ```


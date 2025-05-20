# API Design Thinking #

## **WORK-IN-PROGRESS** ##

We assemble some thinking about the general structure of APIs for [Moltrack](https://github.com/datagrok-ai/mol-track)

What are the strategies?  Typical RESTful API design where there is a close connection between the underlying database design and the endpoint design.  For each database table, there are corresponding CRUD needs fulfilled by CRUD endpoints. (CRUD - create, read, update, delete).  The corresponding REST method are POST, GET, PUT or PATCH, and DELETE.  There is an open standards, that Revvity's API is modeled after.

We have a number of intents.

1. set up system
2. set up properties, synonyms, associations
3. register virtual compounds
4. register batches
5. set up assay_types and corresponding assay_type_details and assay_type_properties
6. register assays and assay_results.
7. search within a domain (compounds, batches, assay_results) or in combination

## Schema - WIP ##

- `POST /schema`

## Register Batches ##

Batch registration will be performed as singletons or in bulk, by the chemist or registrar.  Use cases include supporting individual batch synthesis, library synthesis generated internally or by a contract research organization (CRO) and acquistion from a vendor.  Each registrations entry will contain information about the batch including batch details, batch synonyms, compound information (including synonyms and properties) and batch additions.  Properties, synonym_types, and additons must be defined ahead of time and if missing will result in a registration failure.  The inputs will be a collection of batch records that will be presented in a CSV-style format, an SDfile format or a parquet file .  A single registration will be a collection with a single record.  

Registrations could executed synchronously or asynchronously.  In general, registration should be performed asynchronously to accommmodate multiple scenarios where processing would be require a lot of time, including molecules with complex ring structure or a collection that contains many entries.  *For the MVP, we have decided to go with synchronous registration.*

- `POST /batches` - synchronous registration endpoint for CSV input
  - Input
    - query parameters
      - `input_format=csv`  [default csv, optionally mol-v3000]
      - `output_format=csv` [default csv, optional mol-v3000]
    - example in csv

      ```text
      "batch_corporate_id", smiles, common_name, cas, usan
      "EPA-001-001","OC(=O)c1cc2c(cc1)c(=O)O","1,3-benezenedicarboxylic acid","121-91-5"
      "EPA-120-001","c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N",,"169590-42-5","celecoxib"
      ```

    - csv data in the json body

      ```json
      {
         "data": "\"batch_corporate_id\", smiles, common_name, cas, usan\n\"EPA-001-001\",\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",     \"121-91-5\"\n\"EPA-120-001\",\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\"\n",
         "mapping": {
            "smiles": "compounds.smiles",
            "batch_corporate_id": "batches.synonyms.batch_corporate_id",
            "common_name": "compounds.synonyms.common_name",
            "cas": "compounds.synonyms.cas",
            "usan": "compounds.synonyms.usan"
         }
      }
      ```

  - Output
   The batches, batches_detail, batches_synonym, batch_additions, additions, compounds, compound_details, and compound_synonyms tables will be joined appropriately, pivoted and concatenated appropriately to produce a single record.
    - 200

      ```json
      {
         "status_messge" = "Success",
         "data": "\"batch_corporate_id\",smiles,common_name,cas,usan,batch_id,compound_id,compound_corporate_id\n\"EPA-001-001\",\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",\"121-91-5\",1,2,\"EPA-001\"\n\"EPA-120-001\",\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\",2,4,\"EPA-120\"\n"
      }
      ```

- `GET /batches/` - returns a representation of batches including details, synonyms, additions and compounds including details and synonyms

   1. pagination must be supported
   2. output format must be supported. `output-format=[json|csv|mol-v3000]` with a default of `csv`.  If `json` format is selected, the data is returned in a nested dictionary structure.  If the format is `csv` or `mol-v3000`, then the data is pivoted/concatenated/flattened
   3. there should probably be a hard limit

- `GET /batches/{id}` - .returns a representation of batches including details, synonyms, additions and compounds including details and synonyms
   - Query parameter `output-format=[json|csv|mol-v3000]` with a default of `csv`.  If `json` format is selected, the data is returned in a nested dictionary structure.  If the format is `csv` or `mol-v3000`, then the data is pivoted/concatenated/flattened

- `GET /batches/{id}/properties`
- `GET /batches/{id}/synonyms`
- `GET /batches/{id}/additions` - should returns additions, equivalents, and summary


## Register Virtual Compounds ##

Virtual compound registration will typically be done in bulk, by a cheminformatician or registrar.  These registrations will be contain the structure, synonyms, and possibly properties.  Properties and synonym_types will need to be defined ahead of time.  The inputs will be a collection of compound records that will be presented in a CSV-style format, an SDfile format or a parquet file .  A single registration will be a array with a single record. 

Registrations could executed synchronously or asynchronously.  In general, registration should be performed asynchronously to accommmodate multiple scenarios where processing would be require a lot of time, including molecules with complex ring structure or a collection that contains many entries.  *For the MVP, we have decided to go with synchronous registration.*

- `POST /compounds` - synchronous registration endpoint for CSV input
  - Input
    - query parameters
      - `input_format=csv`  [default csv, optionally mol-v3000]
      - `output_format=csv` [default csv, optional mol-v3000]
    - example in csv

      ```text
      "smiles, common_name, cas, usan
      "OC(=O)c1cc2c(cc1)c(=O)O","1,3-benezenedicarboxylic acid","121-91-5"
      "c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N",,"169590-42-5","celecoxib"
      ```

    - csv data in the json body

      ```json
      {
         "data": "smiles, common_name, cas, usan\n\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",     \"121-91-5\"\n\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\"\n",
         "mapping": {
            "smiles": "compounds.smiles",
            "common_name": "compounds.synonyms.common_name",
            "cas": "compounds.synonyms.cas",
            "usan": "compounds.synonyms.usan"
         }
      }
      ```

  - Output
   The compounds, compound_details, and compound_synonyms tables will be joined appropriately, pivoted and concatenated appropriately to produce a single record.
    - 200

      ```json
      {
         "status_messge" = "Success",
         "data": "smiles,common_name,cas,usan,batch_id,compound_id,compound_corporate_id\n\"OC(=O)c1cc2c(cc1)c(=O)O\",\"1,3-benezenedicarboxylic acid\",\"121-91-5\",2,\"EPA-001\"\n\"c1cc(C)ccc1c2cc(C(F)(F)F)nn2c3ccc(cc3)S(=O)(=O)N\",,\"169590-42-5\",\"celecoxib\",4,\"EPA-120\"\n"
      }
      ```

- `GET /compounds/` - returns an array of compounds with each entry having the comound, the properties and the synonyms.  An optional query parameter will allow the user to select the output format: json, csv-style, sd-file, parquet.  Properties and synonyms will need to be pivoted (and concatenated as necessary) to accomodate the csv,
sd and parquet format expectations.  We will need to have query parameters to support pagination of the resulting output (`start`,`stop`,`max-per-page`)
- `GET /compounds/{compound_id}`
- `GET /compounds/{compound_id}/properties` - I think these will be less used but included for completeness
- `GET /compounds/{compound_id}/synonyms` - I think these will be less used but included for completeness

## Assay Data - WIP ##

`/assay_results`

## Search - WIP ##

Can we provide a generalized search capabilities across all the various entities in our model: compounds, batches, assay_results.  

- `/search/compounds`
- `/search/batches`
- `/search/assay-data`

- `POST /search`
A simple or complex set of search criteria are presented in a json format.  Query parameters provide pagination support and output format style.  User stories:

   1. Find all batches made for project X in the past week.
   2. Find all stereo-similar compounds to I can determine whether I have expressed the stereo nature of my new molecule correctly.
   3. Find all batches made by WuXi for project Y.
   4. Show me all details for batch with a specific corporate id.
   5. Find all batches made for project Y including details and assay results
   6. find all assay-data where
      1. substructure criteria
      2. similarity criteria
      3. batch property of source = WuXi
      4. kinase IC50 < 100 uM


`/search/compounds/structure`
- exact  -->  smiles + standardized by Moltrack settings.  optional pattern of standardization  [array of standardization steps].  uses all layers of the registration mol hash
   - returns - 0 or 1 compound entries, exact match can return multiple tautomers
- less precise  -- substructure, tautomer, no-stereo
substructure + similar to a key compound

`/search/batches`
`/search/assay-data`

logical conditions -> postgres sql condition phrases

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

   A mapping between input fields and database structure will need to be provided or inferred.  The POST method creates a job that needs to be performed and returns the `jobId`.
   For this method we will need to consider the ability to include all the data in the json body, gzipping the content, or provide a file upload mechanism.  The backend job should probably separate the content into individual entries so that we can run in parallel threads and still capture progress and errors.

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
      "batches-collection-job-id": id,
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
            {"CAS Number": "169590-42-5}"
         ],
         "properties": []
      }
   ]
   ```

   A mapping between input fields and database structure will need to be provided or inferred.  The POST method creates a job that needs to be performed and returns the `jobId`.
   For this method we will need to consider the ability to include all the data in the json body, gzipping the content, or provide a file upload mechanism.  The backend job should probably separate the content into individual entries so that we can run in parallel threads and still capture progress and errors.

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
      "compounds-collection-job-id": id,
      "global-status": "In Progress",
      "individual-statuses": [
         { "molecule1": "Success"},
         { "molecule2": "Failed", "error-message": "Business rules failed: chemist not recognized" },
         { "molecule3": "In Progress"}
      ]
   }
   ```


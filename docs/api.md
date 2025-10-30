## API Reference

### Additions

| Method | Path | Description |
|--------|------|-------------|
| GET | /v1/additions/ | Read Additions V1 |
| POST | /v1/additions/ | Create Additions |
| GET | /v1/additions/salts | Read Additions Salts V1 |
| GET | /v1/additions/solvates | Read Additions Solvates V1 |
| GET | /v1/additions/{addition_id} | Read Addition V1 |
| PUT | /v1/additions/{addition_id} | Update Addition V1 |
| DELETE | /v1/additions/{addition_id} | Delete Addition |

### Admin

| Method | Path | Description |
|--------|------|-------------|
| PATCH | /v1/admin/settings | Update Settings |
| PATCH | /v1/admin/update-standardization-config | Update Standardization Config |

### Assay results

| Method | Path | Description |
|--------|------|-------------|
| POST | /v1/assay_results/ | Create Assay Results |

### Assay runs

| Method | Path | Description |
|--------|------|-------------|
| POST | /v1/assay_runs/ | Create Assay Runs |
| GET | /v1/assay_runs/ | Get Assay Runs |
| GET | /v1/assay_runs/{assay_run_id} | Get Assay Run By Id |

### Assays

| Method | Path | Description |
|--------|------|-------------|
| POST | /v1/assays | Create Assays |
| GET | /v1/assays/ | Get Assays |
| GET | /v1/assays/{assay_id} | Get Assay By Id |

### Auto-map-columns

| Method | Path | Description |
|--------|------|-------------|
| POST | /v1/auto-map-columns | Auto Map Columns |

### Batches

| Method | Path | Description |
|--------|------|-------------|
| GET | /v1/batches | Get Batch By Any Synonym |
| POST | /v1/batches/ | Register Batches |
| GET | /v1/batches/ | Get Batches |
| GET | /v1/batches/additions | Get Batch Additions |
| GET | /v1/batches/properties | Get Batch Properties |
| GET | /v1/batches/synonyms | Get Batch Synonyms |
| DELETE | /v1/batches/{corporate_batch_id} | Delete Batch By Any Synonym |

### Compounds

| Method | Path | Description |
|--------|------|-------------|
| GET | /v1/compounds | Get Compound By Any Synonym |
| POST | /v1/compounds/ | Register Compounds |
| GET | /v1/compounds/ | Get Compounds |
| GET | /v1/compounds/properties | Get Compound Properties |
| GET | /v1/compounds/synonyms | Get Compound Synonyms |
| PUT | /v1/compounds/{corporate_compound_id} | Update Compound By Id |
| DELETE | /v1/compounds/{corporate_compound_id} | Delete Compound By Id |

### Schema

| Method | Path | Description |
|--------|------|-------------|
| GET | /v1/schema/ | Get Schema |
| POST | /v1/schema/ | Preload Schema |
| GET | /v1/schema/batches | Get Schema Batches |
| GET | /v1/schema/batches/synonyms | Get Schema Batch Synonyms |
| GET | /v1/schema/compounds | Get Schema Compounds |
| GET | /v1/schema/compounds/synonyms | Get Schema Compound Synonyms |

### Schema-direct

| Method | Path | Description |
|--------|------|-------------|
| GET | /v1/schema-direct/ | Get Schema Direct |

### Search

| Method | Path | Description |
|--------|------|-------------|
| POST | /v1/search/assay-results | Search Assay Results Advanced |
| POST | /v1/search/assay-runs | Search Assay Runs Advanced |
| POST | /v1/search/assays | Search Assays Advanced |
| POST | /v1/search/batches | Search Batches Advanced |
| POST | /v1/search/compounds | Search Compounds Advanced |
| POST | /v1/search/generate-filter | Generate Search Filter |

### Validators

| Method | Path | Description |
|--------|------|-------------|
| GET | /v1/validators/ | Get Validators |
| POST | /v1/validators/ | Register Validators |
| DELETE | /v1/validators/{validator_name} | Delete Validator By Name |


# Database Readme

## ~~Moltrack schema generation~~

~~The design is evolved with DBML in [schema.dbml](schema.dbml).  Install `npm install -g @dbml/cli`~~

~~`dbml2sql schema.dbml -o schema.sql`~~

## Instantiation

### Prerequisites

1. Available postgress engine
2. user account to the engine with necessary permissions to create a database, database objects and grant permissions
3. Appropriate database client like psql, pgadmin, dbeaver.

### Initialization

You can directly instatiate with 

1. Execute [db.sql](./db.sql)
1. Execute [schema.sql](./schema.sql)
1. Execute [schema_moltrack_privileges.sql](./schema_moltrack_privileges.sql)
1. Execute [schema_rdkit.sql](./schema_rdkit.sql)
1. Execute [schema_rdkit_privileges.sql](./schema_rdk_privileges.sql)


Or you can deploy via the provided [Dockerfile](../Dockerfile)

You can additionally deploy the [moltrack_indexes](moltrack_indexes.sql) to add performance when following foreign key relationships.
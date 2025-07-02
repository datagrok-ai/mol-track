FROM informaticsmatters/rdkit-cartridge-debian:Release_2024_03_6

ENV POSTGRES_USER=postgres
ENV POSTGRES_PASSWORD=postgres
ENV POSTGRES_DB=moltrack

ADD db/db.sql /docker-entrypoint-initdb.d/01_db.sql
ADD db/schema.sql /docker-entrypoint-initdb.d/02_schema.sql
Add db/setup.sql /docker-entrypoint-initdb.d/03_setup.sql
ADD db/schema_rdkit.sql /docker-entrypoint-initdb.d/04_schema_rdkit.sql
ADD db/schema_moltrack_privileges.sql /docker-entrypoint-initdb.d/05_schema_moltrack_privileges.sql
ADD db/schema_rdk_privileges.sql /docker-entrypoint-initdb.d/06_schema_rdk_privileges.sql

RUN sed -i 's/\^/\;/g' /docker-entrypoint-initdb.d/04_schema_rdkit.sql

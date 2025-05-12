FROM mcs07/postgres-rdkit

ENV POSTGRES_USER=postgres
ENV POSTGRES_PASSWORD=postgres
ENV POSTGRES_DB=moltrack

ADD db/db.sql /docker-entrypoint-initdb.d/
ADD db/schema.sql /docker-entrypoint-initdb.d/
ADD db/schema_moltrack_privileges.sql /docker-entrypoint-initdb.d/
ADD db/schema_rdkit.sql /docker-entrypoint-initdb.d/
ADD db/schema_rdk_privileges.sql /docker-entrypoint-initdb.d/

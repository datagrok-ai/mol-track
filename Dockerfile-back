FROM mcs07/postgres-rdkit

ENV POSTGRES_USER=postgres
ENV POSTGRES_PASSWORD=postgres
ENV POSTGRES_DB=moltrack

ADD db/schema.sql /docker-entrypoint-initdb.d/ 
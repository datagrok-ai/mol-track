FROM mcs07/postgres-rdkit

ENV POSTGRES_USER=postgres
ENV POSTGRES_PASSWORD=postgres
ENV POSTGRES_DB=moltrack

ADD db/db.sql /docker-entrypoint-initdb.d/db.sql
ADD db/schema.sql /docker-entrypoint-initdb.d/schema.sql

RUN sed -i 's/\^/\;/g' /docker-entrypoint-initdb.d/schema.sql
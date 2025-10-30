# Architecture

MolTrack server component consists of the FastAPI service that uses RDKit, and 
a Postgres database with the RDKit cartridge for chemical intelligence.

## Property calculations

Basic properties are calculated via RDKit. More advanced properties can be calculated
via webhooks (events triggered when compounds are registered).

## Security

For the role and privilege management system, MolTrack uses the same tables as Datagrok.
This allows to either host MolTrack on its own (in this case, the necessary tables get deployed to the
MolTrack schema), or use tables already present in the Datagrok installation. In the latter case,
the UI for role management comes out of the box.
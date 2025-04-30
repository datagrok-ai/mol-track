# MolTrack Minimal Viable Product

We describe here the capabilities that we see as appropriate for a minimal viable product that can be used to test the chemistry, alias, property, and assay features and provide feedback on evolution and hackability.  We use the terminology of a role (implementer, admin, chemist, biologist).  *In this MVP, this will be a single global trusted user.*

1. Implementor can
    1. deploy postgres database with RDkit
    2. configure business rules
        1. compound standardization rules
        2. compound uniqueness rules
        3. compound identification rules and aliases / synonym
        4. batch identification rules and aliases / synonym
    3. configure
        1. required compound properties
        2. optional compound properties
        3. required batch properties
        4. optional batch properties
        5. compound level alias types
        6. batch level alias types
    4. deploy mol-track service

2. Users can be managed by the system or related systems
    1. Run in trusted manner with just 1 global user.
    2. ~~Users can be exclusively managed within MolTrack~~
    3. ~~Users managed by Datagrok can be used within MolTrack~~

3. Permissions can be managed by the system or related systems
    1. Global trusted user will have all permissions.
    2. ~~Role based permissions can be managed within MolTrack.~~
    3. ~~Permissions managed by Datagrok can be used within MolTrack.~~

4. Admin can *(for MVP, global user)*
    1. manage controlled vocabularies - add allowed values to database tables.

5. Users can *(for MVP, global user)*
    1. Search for compounds by exact search, substructure search, no stereo, tautomer, similarity
    2. Search for compounds by aliases.
    3. Search for batches based on compounds using above criteria.
    4. Search for batches based on aliases.
    5. Search for batches based on properties.
6. Chemists can *(for MVP, global user)*
    1. Register compounds by smiles+json or sd-file (ctab + properties)
    2. Register batches  with compound, compound alias, batch, batch aliases, batch properties with smiles+json or sd-file with v3000-ctab and properties
    3. Update compound based on business rules and permissions
    4. Update batches based on business rules and permissions
7. Biologists can *(for MVP, global user)*
    1. configure result_types (properties).
    2. configure assay_types and property associations.
    3. register assays
    4. register assay_results
8. Biologists and chemists can *(for MVP, global user)*
    1. search, list assay_types and associated properties
    2. retrieve assay_results based on compound criteria and/or batch criteria and/or assay_type criteria
9. UX
    1. Swagger based API calling
    2. python command line client

# MolTrack Minimal Viable Product

We describe here the capabilities that we see as appropriate for a minimal viable product that can be used to test the chemistry, synonym, property, and assay features and provide feedback on evolution and hackability.  We use the terminology of a role (implementer, admin, chemist, biologist).  *In this MVP, this will be a single global trusted user.*

1. Implementor can
    1. deploy postgres database with RDkit
    2. configure business rules
        - [x] compound standardization rules
        - [ ] compound uniqueness rules
        - [ ] compound identification rules and synonym
        - [ ] batch identification rules and synonym
    3. configure
        - [ ] required compound properties
        - [x] optional compound properties
        - [ ] required batch properties
        - [x] optional batch properties
        - [x] compound level synonym types
        - [x] batch level synonym types
    4. deploy mol-track service

2. Users can be managed by the system or related systems
    - [ ] Run in trusted manner with just 1 global user.
    - ~~Users can be exclusively managed within MolTrack~~
    - ~~Users managed by Datagrok can be used within MolTrack~~

3. Permissions can be managed by the system or related systems
    1. Global trusted user will have all permissions.
    2. ~~Role based permissions can be managed within MolTrack.~~
    3. ~~Permissions managed by Datagrok can be used within MolTrack.~~

4. Admin can *(for MVP, global user)*
    1. manage controlled vocabularies - add allowed values to database tables.

5. Users can *(for MVP, global user)*
    - [] Search for compounds by exact search, substructure search, no stereo, tautomer, similarity
    - [] Search for compounds by synonyms.
    - [] Search for batches based on compounds using above criteria.
    - [] Search for batches based on synonyms.
    - [] Search for batches based on properties.
    - [ ] Search for compounds or batches using a combination of search conditions
6. Chemists can *(for MVP, global user)*
    1. [] Register compounds by csv with smiles or ctab ~~or sd-file (ctab + properties)~~
    2. [] Register batches  with compound, compound synonyms, batch, batch synonyms, batch properties with smiles+json or sd-file with v3000-ctab and properties
    3. Update compound based on business rules and permissions
    4. Update batches based on business rules and permissions
7. Biologists can *(for MVP, global user)*
    1. [] configure result_types (properties).
    2. [] configure assay_types and property associations.
    3. [] register assays
    4. [] register assay_results
8. Biologists and chemists can *(for MVP, global user)*
    1. search, list assay_types and associated properties
    2. retrieve assay_results based on compound criteria and/or batch criteria and/or assay_type criteria
9. UX
    1. [] Swagger based API calling
    2. python command line client

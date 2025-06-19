# MolTrack Minimal Viable Product

We describe here the capabilities that we see as appropriate for a minimal viable product that can be used to test the chemistry, synonym, property, and assay features and provide feedback on evolution and hackability.  We use the terminology of a role (implementer, admin, chemist, biologist).  *In this MVP, this will be a single global trusted user.*

1. Implementor can
    - [x] deploy postgres database with RDkit
    - [ ] configure business rules
        - [x] compound standardization rules
        - [ ] compound uniqueness rules
        - [ ] compound identification rules and synonym
        - [ ] batch identification rules and synonym
    - [ ] configure
        - [ ] required compound properties
        - [x] optional compound properties
        - [ ] required batch properties
        - [x] optional batch properties
        - [x] compound level synonym types
        - [x] batch level synonym types
        - [x] additions (salts and solvates)
    - [ ] deploy mol-track service

2. Users can be managed by the system or related systems
    - [x] Run in trusted manner with just 1 global user.
    - ~~Users can be exclusively managed within MolTrack~~
    - ~~Users managed by Datagrok can be used within MolTrack~~

3. Permissions can be managed by the system or related systems
    - [x] Global trusted user will have all permissions.
    - ~~Role based permissions can be managed within MolTrack.~~
    - ~~Permissions managed by Datagrok can be used within MolTrack.~~

4. Admin can *(for MVP, global user)*
    - [ ] manage controlled vocabularies - add allowed values to database tables.

5. Chemists can *(for MVP, global user)*
    - [x] Register compounds by csv with smiles or ctab including compounds, properties, synonyms. ~~or sd-file (ctab + properties)~~
    - [x] Register batches  with compound, compound synonyms, batch, batch synonyms, batch properties, batch additions with csv, ~~smiles+json or sd-file with v3000-ctab and properties~~
    - [x] Update compound based on business rules and permissions
      - Hard coded business rule to specify that only compounds without batches can be updated.
    - [x] Delete compounds based on business rules
    - [ ] Update batches based on business rules and permissions
    - [ ] Delete batches

6. Biologists can *(for MVP, global user)*
    - [x] configure result_types (properties).
    - [x] register assays and property associations.
    - [x] register assay_runs
    - [x] register assay_results

7. Users can *(for MVP, global user)*
    - [ ] Search for compounds by exact search, substructure search, no stereo, tautomer, similarity
    - [ ] Search for compounds by synonyms.
    - [ ] Search for batches based on compounds using above criteria.
    - [ ] Search for batches based on synonyms.
    - [ ] Search for batches based on properties.
    - [ ] Search for compounds or batches using a combination of search conditions
    - [ ] Search, list assays and associated properties
    - [ ] Retrieve assay_results based on compound criteria and/or batch criteria and/or assays, assay_runs criteria

8. UX
    - [x] Swagger based API calling
    - [ ] python command line client

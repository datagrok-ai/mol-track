# Features

Work in progress:
1. Compound Registration
    * [x] **Unique Identifiers**: Automatically assigns unique identifiers (e.g., registration numbers, UUIDs) to new compounds.
    * [x] **Duplicate Detection**: Prevents registration of duplicates.
    * [x] **Structure Validation**: Valence checking, standardization of stereochemistry, etc.
    * [ ] **Structure Standardization**: Converts entered structures into a consistent format, handling tautomerization, salts, and stereo conventions.
2. Metadata
    * [ ] **Custom Attributes**: Supports capturing custom metadata (e.g., biological data, physicochemical properties, origin information) and ties it to the appropriate entity (compound, batch/lot).
    * [ ] **Attachment Management**: Allows attaching documents (NMR spectra, mass spectrometry data, analytical certificates).
3. Batches and Lots
    * [x] **Batch Registration**: Manages registration of multiple batches or lots for a single compound.
    * [x] **Duplicate Detection**: Prevents the registration of duplicates
    * [ ] **Purity and Inventory Tracking**: Tracks batch-specific details such as purity, quantity, storage location, supplier, and expiration dates.
4. Protocols and Assay Results
    * [ ] **Protocols**: Define assay types used to measure batches.
    * [ ] **Assay Results**: Register and query assay results.
5. Search
    * [x] **Structure-based Search**: Supports exact, substructure, similarity, and Markush searches.
    * [x] **Metadata Search**: Enables querying by metadata fields such as IDs, names, properties, and batch information.
6. Audit and Compliance
    * [ ] **Audit Trails**: Records detailed logs of registration, editing, and deletion activities for compliance and traceability.
    * [ ] **Role-based Access Control**: Implements security controls to ensure sensitive data is accessible only by authorized users.
7. Integration and APIs
    * [x] **API Access**: Provides RESTful APIs to facilitate integration with other lab informatics systems (ELNs, LIMS, inventory management systems).
9. User Interface
    * [ ] **Chemical Drawing Integration**: Allows users to input structures directly using chemical drawing tools (e.g., MarvinJS, ChemDraw, Ketcher).
    * [ ] **Custom Reports**: Generates reports on compound libraries, registration statistics, and inventory statuses.
    * [ ] **Visualization Tools**: Includes dashboards and data visualization features for quick analysis and decision-making.
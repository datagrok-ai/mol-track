## Concepts
- As a compound registration system, meaning registering and tracking the parent structure and obtaining a unique corporate id for the structure.  
- As a batch registration system, where the batch entity is the primary entity and compound represents a structural annotation for the batch.
- Salts and/or solvents are an additional annotation on the batch, but all serve as a point of aggregation - 'find all salt/solvent form for batches of this structure.'
- The inital focus is small molecules and no-structs, but the system should be evolvable to handle other modalities: peptides, siRNA, ASO, biologics.


## Persona / Roles
- Administrator has all the appropriate rights to administer the system
- Chemist (or editor) has the ability to add batches and compounds and update those that created.
- Registar has the ability to add batches and compounds in bulk and edit all compounds and batches.
- Viewer has the ability to search and retrieve compounds, batches/lots and appropriate properties


## User Stories

Here are the user stories that MolTrack should support:

- As an administrator, I can add a chemist to the user management system, so they may register their batches.
- As a chemist, I want to add a batch of a compound that I have synthesized and receive an appropriate batch identifier and compound identifier so that the biologist can generate identifiable unique assay results.
- As a chemist, I want to edit a batch that I have registered to correct or add additional information or update the structure.
- As a chemist, I should only be allowed to update the compound structure if I own all the associated batches, so that I do not perturb the registration number for batches/compounds made by others.
- As a chemist, the molecule coordinates that I use to register should be saved to preserve an appropriate orientation.
- 
- As a registrar, I want to add batches (and compounds) in bulk to support the case when we might acquire a set of batches and compounds for research purposes.
- As a registrar, I want to edit the batches (and compounds) that others have registered to correct errors.
- As a biologist, I want to view all the information for a specific batch so that I might prepare my experiment correctly.
- As an administrator, I want to configure the business rules for batch uniquenss verification.
- As an administrator, I want to configure the business rules for compound uniqueness verification.
- As an administrator, I want to configure the business rules for no-structure uniqueness.
- As a chemist, I want to register a batch with a no-structure compound and receive appropriate batch and compound identifiers.
- As a chemist or registrar, I want to annotate my registration with aliases, so that I can capture external identifiers.
- As chemist, I want to query compounds using substructure, similarity, exact, no-stereo exact or tautomeric methods, so that I might retrieve the exact subset that is scientifically interesting.
- As a scientist, I want to query batches by properties and/or aliases to find the relevant batches.
- As a scientist, I want to query compounds by properties and/or aliases to find the relevant compounds.
- As a peptide chemist, I want to register my compounds using a helm representation, 
- As an administrator, I want to configure my system so that certain properties must adhere to controlled vocabularies, in order to maintain data cleanliness.
- As the system, I should be able to provide a preferred rendering of the compound as an image or svg, so that I can used in downstream systems that don't have native molecular rendering.
- As a chemist or registrar, I should be able to associate a preferred rendering for compound, so that consumers will be able to retrieve and benefit from the preferred rendering.
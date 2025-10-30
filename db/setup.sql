BEGIN;

INSERT INTO moltrack.users (
  id, email, first_name, last_name,
  has_password, is_active, is_service_account,
  created_by, updated_by
) VALUES (
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d',
  'admin@datagrok.ai', 'Admin', 'Admin',
  true, true, true,
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d',
  '3f5b8c3e-1a72-4c09-9aeb-2f12a7a81e8d'
);

INSERT INTO moltrack.semantic_types (name, description) 
VALUES
    ('Synonym', 'A semantic type representing a synonym or alternative identifier'),
    ('Molecule', 'A semantic type that represents the chemical structure of a compound')
ON CONFLICT (name) DO NOTHING;

INSERT INTO moltrack.settings (name, value, description)
VALUES
    ('corporate_compound_id_pattern', 'DG-{:06d}', 'Pattern for corporate compound IDs'),
    ('corporate_batch_id_pattern', 'DGB-{:06d}', 'Pattern for corporate batch IDs');

INSERT INTO moltrack.settings (name, value, description)
VALUES
    ('corporate_compound_id_friendly_name', 'Grok ID', 'Friendly name for corporate compound IDs'),
    ('corporate_batch_id_friendly_name', 'Grok Batch ID', 'Friendly name for corporate batch IDs');

INSERT INTO moltrack.settings (name, value, description)
VALUES
    ('compound_sequence_start', '1', 'Starting value for the molregno sequence'),
    ('batch_sequence_start', '1', 'Starting value for the batchregno sequence');

with ADMIN AS (
  SELECT id FROM moltrack.users WHERE email = 'admin@datagrok.ai'
),
STYPE AS (
  SELECT id FROM moltrack.semantic_types WHERE name = 'Synonym'
)
INSERT INTO moltrack.properties (created_by, updated_by, name, description, value_type, semantic_type_id, property_class, entity_type, pattern, friendly_name)
VALUES (
  (SELECT id FROM ADMIN),
  (SELECT id FROM ADMIN),
  'corporate_compound_id', 'Official institution synonym for compounds',
  'string', 
  (SELECT id FROM STYPE), 
  'DECLARED', 'COMPOUND',
  (SELECT value FROM moltrack.settings WHERE name = 'corporate_compound_id_pattern'),
  (SELECT value FROM moltrack.settings WHERE name = 'corporate_compound_id_friendly_name')
), (
  (SELECT id FROM ADMIN),
  (SELECT id FROM ADMIN),
  'corporate_batch_id', 'Official institution synonym for batches',
  'string', 
  (SELECT id FROM STYPE), 
  'DECLARED', 'BATCH',
  (SELECT value FROM moltrack.settings WHERE name = 'corporate_batch_id_pattern'),
  (SELECT value FROM moltrack.settings WHERE name = 'corporate_batch_id_friendly_name')
);

INSERT INTO moltrack.settings (name, value, description)
VALUES ('Compound Matching Rule',
        'ALL_LAYERS',
        'Defines the rule for matching compounds. Possible values: ALL_LAYERS (default), STEREO_INSENSITIVE_LAYERS, TAUTOMER_INSENSITIVE_LAYERS');

INSERT INTO moltrack.settings (name, value, description)
VALUES (
  'Molecule standardization rules',
  'operations:
  - type: "Cleanup"
    description: "Basic cleanup of the molecule, including removal of explicit hydrogens"
    enable: true

  - type: "FragmentParent"
    description: "Retains the largest parent fragment of the molecule."
    enable: true

  - type: "Uncharger"
    description: "Neutralizes charges on the molecule."
    enable: true
  ',
  'Defines the molecule standardization pipeline'
);

INSERT INTO moltrack.api_keys (
  owner_id, prefix, secret_hash,
  privileges, status, ip_allowlist
) VALUES (
  (SELECT id FROM moltrack.users WHERE email = 'admin@datagrok.ai'),
  'Ib-vwjI',
  'y4XqY3MbIzGb1/aXL8LAKWI4tImIKUOhhn/7bxvf+C4=',
  ARRAY['admin'],
  'active',
  ARRAY[]::cidr[]
);

COMMIT;
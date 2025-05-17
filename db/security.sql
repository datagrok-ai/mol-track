-- Groups of users. They come in two flavors:
-- 1. "wrapper" groups for each individual users. They are created automatically when a user is created (privileges work on a group level). See users.group_id.
-- 2. "regular" groups for grouping users together. 
--This is a subset of the Datagrok "groups" table (and in the case of Datagrok deployment would be used directly).
CREATE TABLE moltrack.groups
(
  id uuid NOT NULL primary key,
  name varchar(1024),           -- name to be used in the sctipting and URLs
  description varchar(16384),   -- description of the group
  personal boolean,             -- if true, this is just a "wrapper" for the user
  author_id uuid,               -- id of the user who created the group
  created_on timestamp without time zone,
  updated_on timestamp without time zone
);

-- Users table. 
-- This is a subset of the Datagrok "users" table (and in the case of Datagrok deployment would be used directly).
CREATE TABLE moltrack.users
(
  id uuid NOT NULL primary key,
  email varchar(64),
  first_name varchar(64),
  last_name varchar(64),
  login varchar(64),
  group_id uuid references moltrack.groups(id),  -- "wrapper" group for the user
  status varchar(255),
  joined timestamp without time zone,
  has_password boolean,
  email_confirmed boolean,
  is_service bool default false
);

-- Groups could be nested.
-- User privileges is a sum of all privileges from all groups a user belongs to.
-- This is a subset of the Datagrok "groups_relations" table (and in the case of Datagrok deployment would be used directly).
CREATE TABLE moltrack.groups_relations
(
  id uuid primary key,
  parent_id uuid references moltrack.groups(id),
  child_id uuid references moltrack.groups(id),
  is_admin boolean  -- if true, child group (usually a person) is an admin of the group
);

-- Used for grouping entities together
-- This is a subset of the Datagrok "projects" table (and in the case of Datagrok deployment would be used directly).
CREATE TABLE moltrack.projects
(
  id uuid primary key,
  name varchar(1024),
  description varchar(16384),
  author_id uuid references moltrack.users(id),
  created_on timestamp without time zone,
  updated_on timestamp without time zone
);

-- Types of entities (such as compounds, batches, assays, etc.).
CREATE TABLE moltrack.entities_types
(
  id uuid NOT NULL primary key,
  name varchar(512) NOT NULL
);

-- Entities (such as compounds, batches, assays, etc.) that are subject to checking roles and permissions.
CREATE TABLE moltrack.entities
(
  id uuid NOT NULL primary key,
  entity_type_id uuid references moltrack.entities_types(id),
  name varchar(1024),
  is_deleted bool default false
);

-- Connects entities to projects. Projects could form a hierarchy.
-- This is a subset of the Datagrok "project_relations" table (and in the case of Datagrok deployment would be used directly).
create table moltrack.project_relations (
  id uuid primary key,
  project_id uuid references moltrack.projects(id),
  entity_id uuid references moltrack.entities(id)
);

-- Permissions for entities types (such as "change structure" for "compounds")
CREATE TABLE moltrack.entities_types_permissions
(
  id uuid NOT NULL primary key,
  entity_type_id uuid references moltrack.entities_types(id),
  name varchar(1024),
  description varchar(65536),
  project_permission_id uuid references moltrack.entities_types_permissions(id)
);


-- Granted permissions for entities to groups
CREATE TABLE moltrack.permissions
(
  id uuid NOT NULL primary key,
  entity_id uuid,
  user_group_id uuid references moltrack.groups(id),
  permission_id uuid references moltrack.entities_types_permissions(id)
);
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_database WHERE datname = 'moltrack') THEN
        CREATE DATABASE moltrack;
    END IF;
END $$;

GRANT ALL PRIVILEGES ON DATABASE moltrack TO CURRENT_USER;

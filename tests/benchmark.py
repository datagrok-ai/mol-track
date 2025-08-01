import csv
import random
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterator

from sqlalchemy import text
from sqlalchemy.engine import Connection

from app.main import get_db


FIRST_NAMES = ["Alice", "Bob", "Charlie", "Diana", "Edward", "Fiona"]
LAST_NAMES = ["Smith", "Johnson", "Williams", "Brown", "Jones", "Garcia"]
DOMAINS = ["example.com", "test.org", "sample.net"]


def generate_random_name() -> str:
    first = random.choice(FIRST_NAMES)
    last = random.choice(LAST_NAMES)
    return f"{first} {last}"


def generate_email(name: str) -> str:
    username = name.lower().replace(" ", ".")
    return f"{username}@{random.choice(DOMAINS)}"


def generate_timestamp() -> str:
    now = datetime.now()
    delta = timedelta(seconds=random.randint(-1_000_000, 1_000_000))
    return (now + delta).isoformat(sep=" ")


def generate_data(record_count: int, output_csv_path: Path) -> None:
    print(f"Generating {record_count} records...")
    with output_csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "age", "email", "created_at"])
        for _ in range(record_count):
            name = generate_random_name()
            age = random.randint(18, 90)
            email = generate_email(name)
            timestamp = generate_timestamp()
            writer.writerow([name, age, email, timestamp])
    print(f"Data written to {output_csv_path}")


def get_create_table_sql() -> str:
    return """
CREATE TABLE IF NOT EXISTS benchmark_data (
    id SERIAL PRIMARY KEY,
    name TEXT,
    age INTEGER,
    email TEXT,
    created_at TIMESTAMP
);
"""


def truncate_table(db: Connection, table_name: str = "benchmark_data") -> None:
    db.execute(text(f"TRUNCATE TABLE {table_name}"))
    db.commit()


def insert_data_from_csv(db, csv_path):
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    print(f"Preparing {len(rows)} records for bulk insert...")

    columns = ["name", "age", "email", "created_at"]

    values_sql = ",\n".join(
        [
            f"('{row['name'].replace("'", "''")}', {int(row['age'])}, '{row['email'].replace("'", "''")}', '{row['created_at']}')"
            for row in rows
        ]
    )

    insert_sql = f"""
    INSERT INTO benchmark_data ({", ".join(columns)})
    VALUES {values_sql}
    """

    start = time.perf_counter()
    db.execute(text(insert_sql))
    db.commit()
    end = time.perf_counter()
    duration = end - start
    print(f"Inserted {len(rows)} records in {end - start:.3f} seconds")
    return duration


def copy_data_from_csv(db: Connection, csv_path: Path) -> float:
    abs_path = f"/data/{csv_path}"
    copy_sql = f"""
    COPY benchmark_data(name, age, email, created_at)
    FROM '{abs_path}'
    WITH CSV HEADER
    """
    print(f"Using COPY command to import from {abs_path} ...")
    start = time.perf_counter()
    db.execute(text(copy_sql))
    db.commit()
    end = time.perf_counter()
    duration = end - start
    print(f"Inserted using COPY in {duration:.3f} seconds")
    return duration


def run_benchmarks():
    db_gen: Iterator[Connection] = get_db()
    db = next(db_gen)

    try:
        print("Creating table...")
        db.execute(text(get_create_table_sql()))
        db.commit()
        print("Table created successfully.")

        powers = range(1, 8)
        csv_path = Path("benchmark_data.csv")

        for power in powers:
            record_count = 10**power
            print(f"\n=== Benchmark for 10^{power} records ===")
            generate_data(record_count, csv_path)

            truncate_table(db)
            insert_time = insert_data_from_csv(db, csv_path)

            truncate_table(db)
            copy_time = copy_data_from_csv(db, csv_path)

            faster = "COPY" if copy_time < insert_time else "INSERT"
            slower = "INSERT" if faster == "COPY" else "COPY"
            percent = (abs(copy_time - insert_time) / max(copy_time, insert_time)) * 100

            print(f"{faster} was faster than {slower} by {percent:.2f}%")
    finally:
        db.close()


if __name__ == "__main__":
    run_benchmarks()

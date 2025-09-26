import csv
import typer
from dateutil.parser import parse as date_parse


csv_app = typer.Typer()


# TODO: Ask whether we really want to keep this command
@csv_app.command("analyze")
def analyze_csv_columns(
    input_csv: str = typer.Argument(..., help="Path to the input CSV file"),
    output_csv: str = typer.Argument(..., help="Path to the output CSV file"),
):
    """
    Analyze a CSV file's columns and output a CSV with:
    - column name
    - detected type (double, int, datetime, or string)
    - number of null values
    - number of unique values
    """

    def detect_type(values):
        non_nulls = [v for v in values if v.strip() != ""]
        if not non_nulls:
            return "string"
        # Try int
        try:
            for v in non_nulls:
                int(v)
            return "int"
        except Exception:
            pass
        # Try float
        try:
            for v in non_nulls:
                float(v)
            return "double"
        except Exception:
            pass
        # Try datetime, but only if most values parse and values are not all short strings
        parsed_count = 0
        for v in non_nulls:
            try:
                date_parse(v)
                parsed_count += 1
            except Exception:
                pass
        # Only call it datetime if >80% parse and the average value length is at least 8
        if parsed_count / len(non_nulls) > 0.8 and (sum(len(v) for v in non_nulls) / len(non_nulls)) > 8:
            return "datetime"
        return "string"

    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = list(csv.DictReader(f))
        if not reader:
            typer.echo("No data rows found.")
            raise typer.Exit(1)
        columns = reader[0].keys()
        results = []
        for col in columns:
            values = [row[col] for row in reader]
            col_type = detect_type(values)
            null_count = sum(1 for v in values if v.strip() == "")
            unique_count = len(set(v for v in values if v.strip() != ""))
            results.append([col, col_type, null_count, unique_count])
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["column", "type", "null_count", "unique_count"])
        writer.writerows(results)
    typer.echo(f"✅ Analysis complete. Results written to {output_csv}")

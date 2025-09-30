from sqlalchemy import text
from sqlmodel import SQLModel
import typer

from client.config import settings
from client.utils.api_helpers import get_table_row_counts
from client.utils.display import display_data_table, display_database_stats_table

try:
    from app.setup.database import engine, DB_SCHEMA
except ImportError:
    engine = None
    DB_SCHEMA = None

database_app = typer.Typer()


# Database Commands
@database_app.command("stats")
def database_stats():
    """
    Show database table statistics with row counts.
    """
    try:
        # Get row counts for all tables using shared utility
        row_counts = get_table_row_counts()

        # Convert to the format expected by display function
        table_stats = [(DB_SCHEMA, table_name, count) for table_name, count in row_counts.items()]

        # Display results in a table
        display_database_stats_table(table_stats, DB_SCHEMA)

    except Exception as e:
        typer.echo(f"❌ Error connecting to database: {e}", err=True)
        raise typer.Exit(1)


@database_app.command("clean")
def database_clean(
    force: bool = typer.Option(False, "--force", "-f", help="Skip confirmation prompt"),
    url: str = settings.API_BASE_URL,
):
    """
    Clean the database by deleting rows in dependency order.

    All assay, batch, compound, and property data are removed,
    while semantic_types, settings, and users are preserved.
    """

    try:
        # Tables to preserve
        tables_to_preserve = ["semantic_types", "settings", "users"]

        # Tables to clean in dependency order (child tables first)
        tables_in_dependency_order = reversed(SQLModel.metadata.sorted_tables)
        tables_to_clean = [t.name for t in tables_in_dependency_order if t.name not in tables_to_preserve]
        tables_to_clean.append("rdk.mols")  # Include the rdkit mols table

        # Show what will be cleaned
        typer.echo("🗑️  Database Clean Operation")
        typer.echo("=" * 50)
        typer.echo("Tables to be cleaned (in dependency order):")
        for i, table in enumerate(tables_to_clean, 1):
            typer.echo(f"  {i:2d}. {table}")

        typer.echo("\nTables to be preserved:")
        for table in tables_to_preserve:
            typer.echo(f"  ✓ {table}")

        # Get row counts before deletion
        typer.echo("\n📊 Getting current row counts...")
        row_counts = get_table_row_counts(tables_to_clean)

        # Display current row counts in a rich table
        typer.echo("\nCurrent row counts:")
        table_data = [
            {"table": table, "count": row_counts[table]} for table in tables_to_clean if row_counts[table] > 0
        ]

        if table_data:
            # TODO: look for optimizations
            # Create rich table for row counts
            columns = [("Table Name", "cyan", {"no_wrap": True}), ("Row Count", "red", {"justify": "right"})]

            def extract_count_row(item):
                return [item["table"], f"{item['count']:,}"]

            display_data_table(
                data=table_data,
                title="Tables to be cleaned",
                columns=columns,
                row_extractor=extract_count_row,
                show_total=True,
                total_label="TOTAL ROWS TO DELETE",
            )
        else:
            typer.echo("  No tables have data to delete.")

        total_rows = sum(row_counts.values())
        if total_rows == 0:
            typer.echo("\n✅ Database is already clean - no rows to delete.")
            return

        # Confirm unless --force is used
        if not force:
            typer.echo("\n⚠️  WARNING: This operation will permanently delete all data from the specified tables!")
            confirm = typer.confirm("Are you sure you want to proceed?")
            if not confirm:
                typer.echo("❌ Operation cancelled.")
                raise typer.Exit(0)

        # Perform deletion
        typer.echo("\n🧹 Starting database cleanup...")
        deleted_counts = {}

        with engine.connect() as connection:
            for table in tables_to_clean:
                count = row_counts[table]
                if count > 0:
                    try:
                        typer.echo(f"  Deleting from {table}...", nl=False)

                        # Special handling for properties table to preserve corporate IDs
                        if table == "properties":
                            # Delete all properties except corporate_compound_id and corporate_batch_id
                            delete_query = text("""
                                DELETE FROM properties
                                WHERE name NOT IN ('corporate_compound_id', 'corporate_batch_id')
                            """)
                            result = connection.execute(delete_query)
                            deleted_count = result.rowcount
                        else:
                            # Use TRUNCATE for better performance on large tables
                            # TRUNCATE is faster than DELETE as it doesn't generate individual row events
                            delete_query = text(f"TRUNCATE TABLE {table} RESTART IDENTITY CASCADE")
                            connection.execute(delete_query)
                            # TRUNCATE doesn't return rowcount, so use the pre-operation count
                            deleted_count = count

                        deleted_counts[table] = deleted_count
                        typer.echo(f" ✅ {deleted_count:,} rows deleted")
                    except Exception as e:
                        typer.echo(f" ❌ Error: {e}")
                        deleted_counts[table] = 0
                else:
                    typer.echo(f"  Skipping {table} (already empty)")
                    deleted_counts[table] = 0

            # Commit the transaction after all deletions are complete
            connection.commit()

        # Summary
        typer.echo("\n✅ Database cleanup completed!")

        # Create summary table for deleted rows
        summary_data = [
            {"table": table, "deleted": deleted_counts[table]} for table in tables_to_clean if deleted_counts[table] > 0
        ]

        if summary_data:
            # TODO: look for optimizations
            typer.echo("\nSummary of deleted rows:")
            columns = [("Table Name", "cyan", {"no_wrap": True}), ("Rows Deleted", "green", {"justify": "right"})]

            def extract_summary_row(item):
                return [item["table"], f"{item['deleted']:,}"]

            display_data_table(
                data=summary_data,
                title="Deletion Summary",
                columns=columns,
                row_extractor=extract_summary_row,
                show_total=True,
                total_label="TOTAL ROWS DELETED",
            )
        else:
            typer.echo("\nNo rows were deleted.")

    except Exception as e:
        typer.echo(f"❌ Error during database cleanup: {e}", err=True)
        raise typer.Exit(1)

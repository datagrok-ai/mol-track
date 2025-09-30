# client.py
import sys
from pathlib import Path
import typer


from client.cli.batches import batch_app
from client.cli.schema import schema_app
from client.cli.compounds import compound_app
from client.cli.additions import additions_app
from client.cli.assays import assays_app
from client.cli.database import database_app
from client.cli.directory import directory_app
from client.cli.search import search_app
from client.cli.admin import admin_app


# Add parent directory to Python path for imports
# sys.path.append("..")
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))


app = typer.Typer()
app.add_typer(schema_app, name="schema", help="Schema management commands")
app.add_typer(compound_app, name="compounds", help="Compound management commands")
app.add_typer(batch_app, name="batches", help="Batch management commands")
app.add_typer(additions_app, name="additions", help="Addition management commands")
app.add_typer(assays_app, name="assays", help="Assays management commands")
app.add_typer(database_app, name="database", help="Database management commands")
app.add_typer(directory_app, name="directory", help="Directory loading commands")
app.add_typer(search_app, name="search", help="Search functionality")
app.add_typer(admin_app, name="admin", help="Administrative functions")

if __name__ == "__main__":
    app()

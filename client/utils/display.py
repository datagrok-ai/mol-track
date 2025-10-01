import csv
from datetime import datetime
import sys
import typer
from rich.console import Console
from rich.table import Table
from rich.color import ANSI_COLOR_NAMES


def create_rich_table(title: str, columns: list[tuple[str, str, dict]]) -> Table:
    """
    Create a rich table with the specified columns.

    Args:
        title: Table title
        columns: List of (header, style, kwargs) tuples for each column

    Returns:
        Configured Table object
    """
    table = Table(title=title)

    for header, style, kwargs in columns:
        table.add_column(header, style=style, **kwargs)

    return table


def display_data_table(
    data: list[dict],
    title: str,
    columns: list[tuple[str, str, dict]],
    row_extractor: callable,
    show_total: bool = False,
    total_label: str = "TOTAL",
    max_rows: int = None,
) -> None:
    """
    Display data in a rich table format.

    Args:
        data: List of data dictionaries
        title: Table title
        columns: List of (header, style, kwargs) tuples for each column
        row_extractor: Function that extracts row values from a data item
        show_total: Whether to show a total row
        total_label: Label for the total row
        max_rows: Maximum number of rows to display (None = all)
    """
    console = Console()
    table = create_rich_table(title, columns)

    total_value = 0
    total_rows = len(data)
    if max_rows is not None:
        data = data[:max_rows]
    typer.echo(f"Rows: {len(data)} (total available: {total_rows})")

    for item in data:
        row_values = row_extractor(item)
        table.add_row(*row_values)

        # Update total if needed and if the last value is numeric
        if show_total and len(row_values) > 0:
            try:
                total_value += float(str(row_values[-1]).replace(",", ""))
            except (ValueError, TypeError):
                pass

    # Add total row if requested
    if show_total:
        table.add_row("", "")
        table.add_row(total_label, f"{total_value:,.0f}", style="bold")

    console.print(table)


# Table display utility functions
def format_timestamp(timestamp_str: str) -> str:
    """
    Format a timestamp string to a readable format.

    Args:
        timestamp_str: ISO format timestamp string

    Returns:
        Formatted timestamp string or original if parsing fails
    """
    if not timestamp_str:
        return ""

    try:
        # Parse and format the timestamp
        dt = datetime.fromisoformat(timestamp_str.replace("Z", "+00:00"))
        return dt.strftime("%Y-%m-%d %H:%M:%S")
    except (ValueError, AttributeError):
        # If parsing fails, use the original value
        return timestamp_str


def display_compounds_table(compounds_data):
    """Display compounds data in a rich table format."""
    columns = [
        # ("ID", "cyan", {"no_wrap": True}),
        ("Friendly Name", "blue", {"no_wrap": True}),
        ("Common Name", "green", {"no_wrap": True}),
        ("CAS", "yellow", {"no_wrap": True}),
        ("SMILES", "magenta", {"no_wrap": True}),
        ("Synonyms", "white", {"no_wrap": True}),
        ("Created At", "red", {}),
    ]

    def extract_compound_row(compound):
        # Extract common name, CAS, and synonyms from properties
        common_name = ""
        cas_number = ""
        synonyms = []
        friendly_name = ""

        for prop in compound.get("properties", []):
            prop_name = prop.get("name", "").lower()
            prop_value = prop.get("value_string", "")
            semantic_type_id = prop.get("semantic_type_id")

            if prop_name == "common name":
                common_name = prop_value
            elif prop_name == "cas":
                cas_number = prop_value
            elif semantic_type_id == 1 and prop_value:  # Synonyms have semantic_type_id=1
                synonyms.append(f"{prop.get('name', '')}: {prop_value}")
                if prop_name == "corporate_compound_id":
                    friendly_name = prop_value

        # Join synonyms with semicolon
        synonyms_str = "; ".join(synonyms[:3])  # Limit to first 3 synonyms
        if len(synonyms) > 3:
            synonyms_str += f" (+{len(synonyms) - 3} more)"

        # Truncate long fields for better display
        smiles = compound.get("canonical_smiles", "")
        if len(smiles) > 35:
            smiles = smiles[:32] + "..."

        if len(common_name) > 25:
            common_name = common_name[:22] + "..."

        if len(synonyms_str) > 43:
            synonyms_str = synonyms_str[:40] + "..."

        return [
            friendly_name,
            common_name,
            cas_number,
            smiles,
            synonyms_str,
            format_timestamp(compound.get("created_at", "")),
        ]

    display_data_table(data=compounds_data, title="Compounds", columns=columns, row_extractor=extract_compound_row)


def display_batches_table(batches_data):
    """Display batches data in a rich table format."""
    columns = [
        # ("ID", "cyan", {"no_wrap": True}),
        ("Friendly Name", "blue", {"no_wrap": True}),
        ("Compound ID", "green", {"no_wrap": True}),
        ("Synonyms", "white", {"no_wrap": True}),
        ("Notes", "yellow", {"no_wrap": True}),
        ("Created At", "red", {}),
    ]

    def extract_batch_row(batch):
        synonyms = []
        friendly_name = ""
        for prop in batch.get("properties", []):
            prop_name = prop.get("name", "").lower()
            prop_value = prop.get("value_string", "")
            semantic_type_id = prop.get("semantic_type_id")
            if semantic_type_id == 1 and prop_value:  # Synonyms have semantic_type_id=1
                synonyms.append(f"{prop_name}: {prop_value}")
                if prop_name == "corporate_batch_id":
                    friendly_name = prop_value

        synonyms_str = "; ".join(synonyms[:3])  # Limit to first 3 synonyms
        if len(synonyms) > 3:
            synonyms_str += f" (+{len(synonyms) - 3} more)"
        if len(synonyms_str) > 43:
            synonyms_str = synonyms_str[:40] + "..."
        return [
            friendly_name,
            str(batch.get("compound_id", "")),
            synonyms_str,
            batch.get("notes", "")[:30] + "..."
            if batch.get("notes") and len(batch.get("notes", "")) > 30
            else batch.get("notes", ""),
            format_timestamp(batch.get("created_at", "")),
        ]

    display_data_table(data=batches_data, title="Batches", columns=columns, row_extractor=extract_batch_row)


def display_assays_table(assays_data, assay_entity: str = None):
    """Display assays data in a rich table format."""
    columns = [
        ("ID", "cyan", {"no_wrap": True}),
        ("Assay ID", "yellow", {"no_wrap": True}) if assay_entity else None,
        ("Batch ID", "yellow", {"no_wrap": True}) if assay_entity == "result" else None,
        ("Name", "blue", {"no_wrap": True}),
        ("Description", "green", {"no_wrap": True}) if assay_entity != "result" else None,
        ("Created At", "red", {}),
    ]
    columns = [item for item in columns if item]

    def extract_assay_row(assay):
        values = [str(assay.get("id", ""))]
        if assay_entity:
            values.append(str(assay.get("assay", "")["id"]))
        if assay_entity == "result":
            values.append("batch iD")
        values.append(assay.get("name", ""))
        if assay_entity != "result":
            values.append(
                assay.get("description", "")[:30] + "..."
                if assay.get("description") and len(assay.get("description", "")) > 30
                else assay.get("description", "")
            )

        values.append(format_timestamp(assay.get("created_at", "")))

        return values

    title = "Assays" if not assay_entity == "runs" else "Assay Runs"
    display_data_table(data=assays_data, title=title, columns=columns, row_extractor=extract_assay_row)


def extract_value_from_property(prop):
    value_type = prop.get("value_type", "")
    if value_type == "string":
        return prop.get("value_string", "")
    elif value_type in ("int", "double"):
        return prop.get("value_num", "")
    elif value_type == "uuid":
        return prop.get("value_uuid", "")
    elif value_type == "datetime":
        return prop.get("value_datetime", "")
    elif value_type == "bool":
        return prop.get("value_bool")


def display_properties_table(properties_data, max_rows=None, display_value=False):
    """Display properties data in a rich table format."""
    # Define the entity_type order for sorting (matching the API's uppercase format)
    entity_type_order = ["COMPOUND", "BATCH", "ASSAY", "ASSAY_RUN", "ASSAY_RESULT"]

    # Create a function to get entity_type priority for sorting
    def get_entity_type_priority(entity_type):
        try:
            return entity_type_order.index(entity_type.upper())
        except ValueError:
            # If entity_type is not in the predefined order, put it at the end
            return len(entity_type_order)

    # Sort properties by entity_type first, then by name alphabetically
    sorted_properties = sorted(
        properties_data,
        key=lambda prop: (get_entity_type_priority(prop.get("entity_type", "")), prop.get("name", "").lower()),
    )

    columns = [("Name", "cyan", {"no_wrap": True})]

    if not display_value:
        columns.extend(
            [
                ("Entity Type", "magenta", {}),
                ("Value Type", "green", {}),
                ("Semantic Type ID", "yellow", {}),
                ("Created At", "blue", {}),
            ]
        )
    if display_value:
        columns.append(("Value", "white", {}))

    def extract_property_row(prop):
        table_values = [prop.get("name", "")]
        if not display_value:
            table_values.extend(
                [
                    prop.get("entity_type", ""),
                    prop.get("value_type", ""),
                    str(prop.get("semantic_type_id", "")),
                    format_timestamp(prop.get("created_at", "")),
                ]
            )
        if display_value:
            table_values.append(str(extract_value_from_property(prop)))
        return table_values

    display_data_table(
        data=sorted_properties,
        title="Properties",
        columns=columns,
        row_extractor=extract_property_row,
        max_rows=max_rows,
    )


def display_additions_table(resp, max_rows=None):
    """
    Display additions data in a rich table format.
    """
    columns = [
        ("ID", typer.colors.CYAN, {"no_wrap": True}),
        ("Name", typer.colors.BLUE, {"no_wrap": True}),
        ("Role", typer.colors.GREEN, {"no_wrap": True}),
        ("SMILES", typer.colors.MAGENTA, {"no_wrap": True}),
        ("Molecular Formula", typer.colors.YELLOW, {"no_wrap": True}),
        ("Molecular Weight", typer.colors.WHITE, {"no_wrap": True}),
        ("Description", typer.colors.BRIGHT_BLUE, {"no_wrap": True}),
        ("Created At", typer.colors.RED, {"no_wrap": True}),
    ]

    def extract_addition_row(addition):
        smiles = addition.get("smiles", "")
        if smiles and len(smiles) > 35:
            smiles = smiles[:32] + "..."

        desc = addition.get("description", "")
        desc = "" if not desc else desc[:30] + "..."

        return [
            str(addition.get("id", "")),
            addition.get("name", ""),
            addition.get("role", ""),
            smiles,
            addition.get("formula", ""),
            f"{addition.get('molecular_weight', ''):,.2f}" if addition.get("molecular_weight") else "",
            desc,
            format_timestamp(addition.get("created_at", "")),
        ]

    display_data_table(
        data=resp,
        title="Additions",
        columns=columns,
        row_extractor=extract_addition_row,
        max_rows=max_rows,
    )


def display_search_table(resp, max_rows=None):
    """
    Display search results in a table format using rich.
    """
    output = resp["columns"]
    level = resp["columns"][0].split(".")[0]
    all_colors = list(ANSI_COLOR_NAMES)
    columns = [
        (
            col.replace(f"{level}.details.", "").replace(f"{level}.", ""),
            all_colors[i % len(all_colors)],
            {"no_wrap": True},
        )
        for i, col in enumerate(output)
    ]  # Keep original for headers

    def extract_row(item):
        return [str(item[col]) for col in output]

    display_data_table(
        data=resp["data"],
        title=f"Search Results ({level.capitalize()})",
        columns=columns,
        row_extractor=extract_row,
        max_rows=max_rows,
    )


def display_search_csv(data, max_rows=None):
    """
    Display search results in CSV format.
    """
    # Write CSV header
    writer = csv.writer(sys.stdout)
    data = data.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    data = [item for item in data if item.strip() != ""]

    # Write data rows
    total_rows = len(data) - 1  # Exclude header from count
    if max_rows is not None:
        data = data[:max_rows]
    typer.echo(f"Rows: {len(data) - 1} (total available: {total_rows})")

    for row in data:
        writer.writerow(row.split(","))


def display_database_stats_table(table_stats, schema_name):
    """Display database statistics in a rich table format."""
    columns = [("Table Name", "cyan", {"no_wrap": True}), ("Row Count", "green", {"justify": "right"})]

    def extract_stats_row(stat):
        table_name = stat[1]
        row_count = stat[2] or 0
        return [table_name, f"{row_count:,}"]

    display_data_table(
        data=table_stats,
        title=f"Database Statistics - Schema: {schema_name}",
        columns=columns,
        row_extractor=extract_stats_row,
        show_total=True,
    )

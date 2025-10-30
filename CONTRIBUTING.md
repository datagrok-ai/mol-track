# Developer introduction to MolTrack

This document provides an introduction for developers new to the MolTrack project, covering documentation, resources, tooling, setup instructions, and testing.

## MolTrack documentation

For core project information, see:

- [README](./README.md)
- [MVP](./docs/mvp.md)
- [API Design](./docs/api-design.md)
- [Search](./docs/search.md)
- [User Stories](./docs/user-stories.md)

## Resources

### Project links

- [MolTrack repository](https://github.com/datagrok-ai/mol-track)
- [MolTrack GitHub project](https://github.com/orgs/datagrok-ai/projects/14)

### Communication

- Slack channel with bookmarks of relevant information

### Background reading

- [Cheminformatics Beginners Guide](https://neovarsity.org/blogs/cheminformatics-beginners-guide)
- [Cheminformatics: Principles and Applications](https://www.drugdesign.org/chapters/cheminformatics/)
- [Datagrok Cheminformatics](https://datagrok.ai/cheminformatics)
- [RDKit Overview](https://rdkit.org/docs/Overview.html)


## Tooling

- [Commit linting guide](./docs/commit-linting.md)

### Manual setup

Follow these steps to set up MolTrack manually:

#### Create and activate a virtual environment

```bash
python3 -m venv .venv
```

Activate the environment:

* **Windows (CMD):**
  ```cmd
  .venv\Scripts\activate
  ```
* **macOS/Linux:**
  ```bash
  source .venv/bin/activate
  ```

#### Install `uv`

Install `uv` package using pip:

```bash
pip install uv
```

*For alternative installation options, see the [official docs](https://docs.astral.sh/uv/guides/install-python/#getting-started).*

#### Initialize the project environment with `uv`

Set up the virtual environment and dependencies using `uv` commands:

```bash
uv venv
uv sync
```

* `uv venv` creates the virtual environment.
* `uv sync` installs all required dependencies.

#### Configure the database connection

Edit the `database.py` file and update the `SQLALCHEMY_DATABASE_URL` variable with your PostgreSQL connection string:

```python
SQLALCHEMY_DATABASE_URL = "postgresql://user:password@host:port/database"
```

Make sure your database server is running and accessible.

#### Run the server

Start the FastAPI server with:

```bash
uv run --active uvicorn app.main:app --reload
```

You can now access the API at [http://localhost:8000](http://localhost:8000).

### Setting up pytest in VS Code

To configure pytest in VS Code, follow these steps:
1. Install the **Python** extension
   * Open the **Extensions** view (`Ctrl+Shift+X` on Windows/Linux or `Cmd+Shift+X` on macOS).
   * Search for **Python** and install the official extension by Microsoft.
2. Click the **Testing** icon (beaker icon) in the **Activity bar**.
3. Configure python tests
   * Click on **Configure Python Tests** button.
   * When prompted, select:
     * **Test framework**: `pytest`
     * **Test directory**: folder containing the tests (important: ensure it contains an `__init__.py` file — this is required for test discovery to work properly)
Your tests should now be detected and listed in the **Testing panel**.

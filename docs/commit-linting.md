# Commit linting

It's important to follow consistent code conventions to maintain clean, readable, and maintainable code. That’s why we now use [ruff](https://github.com/astral-sh/ruff) — a faster alternative to `flake8`, `black`, and other tools for style checking and formatting.

### Setup

1. Install `pre-commit`:

   ```bash
   pip install pre-commit
   ```

2. Install hooks:

   ```bash
   pre-commit install
   ```

On the first run, it may take some time to install dependencies and set up the hooks. However, on subsequent runs, the hooks will run much faster since they are cached.

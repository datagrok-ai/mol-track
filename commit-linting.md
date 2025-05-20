# Commit linting

It's important to follow the same code conventions to maintain clean, consistent, and readable code. That's why we added tools like `flake8` for style checks, `isort` for sorting imports, and `black` for code formatting.

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

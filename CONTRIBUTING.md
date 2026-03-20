# Contributing

Thanks for considering improvements to SODA.

## Development Setup

```bash
pip install -e .[dev]
```

## Recommended Workflow

1. Create a feature branch.
2. Keep changes focused and small when possible.
3. Run a quick syntax check before opening a pull request.

```bash
python -m compileall src examples
```

4. If tests are added, run them locally before submitting.

## Code Style

- Prefer clear, small public APIs.
- Keep package code separate from paper-specific scripts and outputs.
- Avoid committing generated data, figures, or large result files into the core package tree.

## Pull Request Scope

Good pull requests for this repository include:

- algorithm improvements to the SODA implementation
- packaging and documentation improvements
- reproducible examples for real datasets
- bug fixes in AnnData or CLI integration


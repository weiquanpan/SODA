# GitHub Release Checklist

Use this checklist before making the repository public.

## Repository Metadata

- Choose the final GitHub repository name.
- Add the final repository URL to `pyproject.toml`.
- Add maintainer or author metadata if desired.
- Decide whether the package name should remain `soda-sc`.

## Legal and Citation

- Add a license file.
- Add `CITATION.cff` if you want GitHub citation support.
- Update the citation section in `README.md` when the manuscript metadata is final.

## Packaging

- Confirm `pip install -e .` works in a fresh environment.
- Confirm `soda-run --help` works.
- Confirm the example scripts still match the public API.

## Documentation

- Review `README.md` for repository-specific URLs.
- Make sure the real-data demo instructions match the latest example script.
- Remove any local-only paths before publishing.

## Final Sanity Check

- Do not include datasets unless you explicitly want to distribute them.
- Do not include result folders, figures, caches, or temporary files.
- Confirm `.gitignore` covers local outputs.


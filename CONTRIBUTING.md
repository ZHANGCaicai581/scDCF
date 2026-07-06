# Contributing to scDCF

## Scope

This repository contains both the reusable `scDCF` package and project-specific
research materials. Changes intended for package users should focus on:

- `scDCF/` for library code
- `README.md` and `docs/` for user documentation
- `examples/` for public usage examples
- `tests/` for lightweight package verification

## Development Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Before Opening a Pull Request

Run the lightweight verification suite:

```bash
python tests/run_all_tests.py
```

If you change CLI behavior or outputs, update:

- `README.md`
- `docs/output_structure.md`
- relevant examples in `examples/`

## Style Guidelines

- Keep user-facing defaults simple.
- Prefer deterministic behavior when randomness is involved.
- Avoid exporting intermediate artifacts unless the user explicitly requests them.
- Keep examples aligned with the current public API.

## Reporting Issues

Please include:

- your `scDCF` version
- Python version
- operating system
- the command or API call that failed
- the full traceback or error message

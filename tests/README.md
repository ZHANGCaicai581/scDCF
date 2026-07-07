# Tests

The maintained verification entrypoints are:

```bash
python tests/run_all_tests.py
pytest -q
```

Root-level `pytest -q` is configured to collect from `tests/` only, so local
scratch directories such as `local_tests/` do not affect package verification.

The current maintained tests are:

- `tests/test_pipeline_smoke.py`: end-to-end CLI smoke test on the bundled synthetic dataset
- `tests/test_organization.py`: regression test for the public `organize_output` helper

Files in this directory should be either:

- maintained automated tests wired into `run_all_tests.py`, or
- clearly non-collected utilities with a narrow manual purpose

Collected test files should avoid top-level side effects so that standard
`pytest` discovery remains safe.

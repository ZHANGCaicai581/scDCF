# Tests

The maintained public verification entrypoint is:

```bash
python tests/run_all_tests.py
```

At the moment this runs the lightweight end-to-end smoke test in
`tests/test_pipeline_smoke.py`.

Some other files in this directory are older exploratory or compatibility
checks kept for reference. They should not be treated as the primary release
gate unless they are explicitly updated and wired into `run_all_tests.py`.

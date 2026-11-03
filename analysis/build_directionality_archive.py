"""Build the S7–S8 source archive from an explicit, verified file list."""

from pathlib import Path
import hashlib
import json
import zipfile


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "analysis" / "simulation"
OUTPUT = ROOT / "analysis" / "scDCF_directionality_analysis_2026-11-03.zip"
PREFIX = "scDCF_directionality_analysis/"


def main():
    metadata = json.loads((WORKFLOW / "archive_metadata.json").read_text())
    provenance = json.loads((WORKFLOW / "provenance.json").read_text())
    source_files = provenance["scDCF_source_sha256"]
    for name, expected in source_files.items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit(f"Recorded analysis source changed: {name}")
    if metadata["license"] != "mit":
        raise SystemExit("Archive metadata must agree with the MIT LICENSE.")

    workflow_names = [
        ".gitignore", "README.md", "VALIDATION.txt", "audit_matched.py",
        "generate_scenarios.R", "plot_directionality.R", "provenance.json",
        "reference_metrics.csv", "requirements-scdrs.txt", "requirements.txt",
        "run_benchmark.py", "score_reference_methods.R", "test_audit_matched.py",
        "archive_metadata.json",
    ]
    names = ["LICENSE", "pyproject.toml", "setup.py",
             "analysis/build_directionality_archive.py", *source_files]
    names += [f"analysis/simulation/{name}" for name in workflow_names]
    contents = {name: (ROOT / name).read_bytes() for name in names}
    contents[".zenodo.json"] = (WORKFLOW / "archive_metadata.json").read_bytes()
    contents["VERSION"] = (metadata["version"] + "\n").encode()
    contents["README.md"] = (
        "# scDCF matched directionality analysis\n\n"
        f"Archive version: {metadata['version']}. License: MIT (LICENSE).\n\n"
        "This source snapshot covers Supplementary Tables S7–S8 only.\n"
        "Simulation is implemented in R/Splatter. The included scDCF source\n"
        "matches the hashes recorded for the saved analyses; its existing\n"
        "package version metadata is retained. This analysis version does\n"
        "not claim a new scDCF software release.\n\n"
        "Start with [the workflow and commands](analysis/simulation/README.md).\n"
        "Validation details: analysis/simulation/VALIDATION.txt.\n"
        "SHA256SUMS identifies every other file in this archive.\n"
    ).encode()
    contents["SHA256SUMS"] = "".join(
        f"{hashlib.sha256(data).hexdigest()}  {name}\n"
        for name, data in sorted(contents.items())
    ).encode()
    temporary = OUTPUT.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(contents.items()):
            item = zipfile.ZipInfo(PREFIX + name, (2026, 11, 3, 0, 0, 0))
            item.external_attr = 0o100644 << 16
            item.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(item, data)
    with zipfile.ZipFile(temporary) as archive:
        if archive.testzip() is not None:
            raise SystemExit("ZIP integrity check failed")
    temporary.replace(OUTPUT)
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    OUTPUT.with_suffix(".zip.sha256").write_text(f"{digest}  {OUTPUT.name}\n")
    print(f"{OUTPUT}\nVersion: {metadata['version']}\nFiles: {len(contents)}\nSHA256: {digest}")


if __name__ == "__main__":
    main()

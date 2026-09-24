from pathlib import Path
import json
import subprocess

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLES = PROJECT / "results" / "tables"
FIGURES = PROJECT / "results" / "figures"


def require(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path.relative_to(PROJECT)}")


def main():
    required = [
        PROJECT / "README.md",
        PROJECT / "docs" / "METHODS.md",
        PROJECT / "docs" / "PORTFOLIO_BRIEF.md",
        PROJECT / "scripts" / "01_prepare_data.py",
        PROJECT / "scripts" / "02_run_models_leakage_safe.R",
        PROJECT / "scripts" / "03_make_figures.R",
        TABLES / "model_performance.csv",
        TABLES / "modeling_summary.json",
        TABLES / "train_validation_test_split.csv",
        FIGURES / "held_out_roc_curves.png",
        FIGURES / "standalone_vs_integrated_auc.png",
    ]
    for path in required:
        require(path)

    perf = pd.read_csv(TABLES / "model_performance.csv")
    expected_models = {
        "Diet-only Elastic Net",
        "Metabolomics-only Elastic Net",
        "Late stacked Elastic Net",
        "Block PLS + Elastic Net",
    }
    if set(perf["model"]) != expected_models:
        raise ValueError("Performance table does not contain the expected model set.")
    if set(perf["split"]) != {"validation", "test"}:
        raise ValueError("Performance table must contain validation and test rows.")

    summary = json.loads((TABLES / "modeling_summary.json").read_text())
    leakage_note = summary.get("leakage_fix", "").lower()
    for token in ["imputation", "scaling", "pls", "training split"]:
        if token not in leakage_note:
            raise ValueError(f"Modeling summary leakage note is missing: {token}")

    split = pd.read_csv(TABLES / "train_validation_test_split.csv")
    duplicated_samples = split["Sample"][split["Sample"].duplicated()].tolist()
    if duplicated_samples:
        raise ValueError(f"Samples appear in multiple splits: {duplicated_samples[:5]}")

    raw_tracked = subprocess.run(
        ["git", "ls-files", "data/raw", "data/processed/model_matrix_raw_missing.csv"],
        cwd=PROJECT,
        text=True,
        capture_output=True,
        check=False,
    )
    tracked = [line for line in raw_tracked.stdout.splitlines() if line and not line.endswith(".gitkeep")]
    if tracked:
        raise ValueError(f"Private/raw data files are tracked by Git: {tracked}")

    print("Project audit passed: required files, splits, leakage notes, and Git privacy checks look good.")


if __name__ == "__main__":
    main()

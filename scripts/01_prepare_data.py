from pathlib import Path
import json
import math
import re

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data" / "raw"
PROCESSED = PROJECT / "data" / "processed"
TABLES = PROJECT / "results" / "tables"
PROCESSED.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)


def clean_name(name: str) -> str:
    name = str(name).replace("\u00a0", " ")
    name = name.encode("ascii", "ignore").decode()
    name = re.sub(r"-(Diet|Metabolites)\s*\.0\.0$", "", name, flags=re.I)
    name = re.sub(r"\.(0)\.(0)-(Metabolites)$", "", name, flags=re.I)
    name = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
    return name


def make_unique(names):
    seen = {}
    out = []
    for name in names:
        if name not in seen:
            seen[name] = 0
            out.append(name)
        else:
            seen[name] += 1
            out.append(f"{name}_{seen[name]}")
    return out


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce")
    valid = p.dropna().sort_values()
    out = pd.Series(np.nan, index=p.index, dtype=float)
    m = len(valid)
    if m == 0:
        return out
    adjusted = valid * m / np.arange(1, m + 1)
    adjusted = adjusted.iloc[::-1].cummin().iloc[::-1].clip(upper=1)
    out.loc[adjusted.index] = adjusted
    return out


def normal_two_sided_from_z(z: float) -> float:
    if pd.isna(z):
        return np.nan
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def mann_whitney_normal_p(x0: pd.Series, x1: pd.Series):
    x0 = pd.to_numeric(x0, errors="coerce").dropna()
    x1 = pd.to_numeric(x1, errors="coerce").dropna()
    n0, n1 = len(x0), len(x1)
    if n0 == 0 or n1 == 0:
        return np.nan, np.nan
    combined = pd.concat([x0, x1], ignore_index=True)
    ranks = combined.rank(method="average")
    r1 = ranks.iloc[n0:].sum()
    u1 = r1 - n1 * (n1 + 1) / 2
    mean_u = n0 * n1 / 2
    sd_u = math.sqrt(n0 * n1 * (n0 + n1 + 1) / 12)
    z = (u1 - mean_u) / sd_u if sd_u else np.nan
    return u1, normal_two_sided_from_z(z)


def welch_z_p(x0: pd.Series, x1: pd.Series):
    x0 = pd.to_numeric(x0, errors="coerce").dropna()
    x1 = pd.to_numeric(x1, errors="coerce").dropna()
    n0, n1 = len(x0), len(x1)
    if n0 < 3 or n1 < 3:
        return np.nan, np.nan, np.nan
    m0, m1 = x0.mean(), x1.mean()
    v0, v1 = x0.var(ddof=1), x1.var(ddof=1)
    se = math.sqrt(v0 / n0 + v1 / n1)
    z = (m1 - m0) / se if se else np.nan
    pooled = math.sqrt(((n0 - 1) * v0 + (n1 - 1) * v1) / (n0 + n1 - 2))
    cohen_d = (m1 - m0) / pooled if pooled else np.nan
    return z, normal_two_sided_from_z(z), cohen_d


def main():
    diet_file = RAW / "Diet.xlsx"
    met_file = RAW / "Metabolomics.xlsx"
    if not diet_file.exists() or not met_file.exists():
        raise FileNotFoundError(
            "Place Diet.xlsx and Metabolomics.xlsx in data/raw/ before running this script."
        )

    diet_raw = pd.read_excel(diet_file)
    met_raw = pd.read_excel(met_file)
    diet_raw["casecontrol_clean"] = diet_raw["casecontrol"].astype(str).str.strip()
    met_raw["casecontrol_clean"] = met_raw["casecontrol"].astype(str).str.strip()

    diet = diet_raw.dropna(subset=["Sample"]).rename(columns=lambda c: clean_name(c)).copy()
    met = met_raw.dropna(subset=["Sample"]).rename(columns=lambda c: clean_name(c)).copy()
    diet.columns = make_unique(diet.columns)
    met.columns = make_unique(met.columns)

    diet_features = [c for c in diet.columns if c not in ["Sample", "casecontrol", "casecontrol_clean"]]
    met_features = [c for c in met.columns if c not in ["Sample", "casecontrol", "casecontrol_clean"]]

    label_overlap = met[["Sample", "casecontrol_clean"]].merge(
        diet[["Sample", "casecontrol_clean"]],
        on="Sample",
        how="inner",
        suffixes=("_metabolomics", "_diet"),
    )
    pd.crosstab(
        label_overlap["casecontrol_clean_metabolomics"],
        label_overlap["casecontrol_clean_diet"],
    ).to_csv(TABLES / "target_label_crosstab.csv")

    merged = diet[["Sample", "casecontrol_clean"] + diet_features].merge(
        met[["Sample"] + met_features], on="Sample", how="inner"
    )
    merged["y_crc"] = (merged["casecontrol_clean"] == "Colon Cancer").astype(int)

    x_diet = merged[diet_features].apply(pd.to_numeric, errors="coerce")
    negative = x_diet.lt(0).sum().rename("negative_code_count").reset_index()
    negative = negative.rename(columns={"index": "feature"})
    negative["negative_code_rate"] = negative["negative_code_count"] / len(x_diet)
    negative.sort_values("negative_code_count", ascending=False).to_csv(
        TABLES / "diet_negative_code_audit.csv", index=False
    )
    x_diet = x_diet.mask(x_diet < 0, np.nan)

    x_met = merged[met_features].apply(pd.to_numeric, errors="coerce")

    diet_missing = x_diet.isna().mean().rename("missing_rate").reset_index().rename(columns={"index": "feature"})
    met_missing = x_met.isna().mean().rename("missing_rate").reset_index().rename(columns={"index": "feature"})
    diet_missing.to_csv(TABLES / "diet_missingness.csv", index=False)
    met_missing.to_csv(TABLES / "metabolite_missingness.csv", index=False)

    per_sample = pd.DataFrame(
        {
            "Sample": merged["Sample"],
            "diet_missing_rate": x_diet.isna().mean(axis=1),
            "metabolomics_missing_rate": x_met.isna().mean(axis=1),
        }
    )
    per_sample.to_csv(TABLES / "per_sample_missingness.csv", index=False)

    y = merged["y_crc"]
    diet_tests = []
    for col in x_diet.columns:
        control, crc = x_diet.loc[y == 0, col], x_diet.loc[y == 1, col]
        u, p = mann_whitney_normal_p(control, crc)
        diet_tests.append(
            {
                "feature": col,
                "control_median": control.median(),
                "crc_median": crc.median(),
                "median_difference_crc_minus_control": crc.median() - control.median(),
                "mann_whitney_u": u,
                "p_value_approx": p,
            }
        )
    diet_tests = pd.DataFrame(diet_tests)
    diet_tests["fdr_bh"] = bh_fdr(diet_tests["p_value_approx"])
    diet_tests.sort_values(["fdr_bh", "p_value_approx"]).to_csv(TABLES / "diet_group_tests.csv", index=False)

    met_tests = []
    for col in x_met.columns:
        control, crc = x_met.loc[y == 0, col], x_met.loc[y == 1, col]
        z, p, d = welch_z_p(control, crc)
        met_tests.append(
            {
                "feature": col,
                "control_mean": control.mean(),
                "crc_mean": crc.mean(),
                "mean_difference_crc_minus_control": crc.mean() - control.mean(),
                "welch_z_approx": z,
                "p_value_approx": p,
                "cohen_d": d,
            }
        )
    met_tests = pd.DataFrame(met_tests)
    met_tests["fdr_bh"] = bh_fdr(met_tests["p_value_approx"])
    met_tests.sort_values(["fdr_bh", "p_value_approx"]).to_csv(
        TABLES / "metabolite_group_tests.csv", index=False
    )

    model_matrix = pd.concat(
        [
            merged[["Sample", "casecontrol_clean", "y_crc"]].rename(
                columns={"casecontrol_clean": "casecontrol"}
            ),
            x_diet.add_prefix("diet__"),
            x_met.add_prefix("metabolite__"),
        ],
        axis=1,
    )
    model_matrix.to_csv(PROCESSED / "model_matrix_raw_missing.csv", index=False)

    summary = {
        "matched_samples": int(len(model_matrix)),
        "target_counts": merged["casecontrol_clean"].value_counts().to_dict(),
        "diet_features_raw": int(x_diet.shape[1]),
        "metabolite_features_raw": int(x_met.shape[1]),
        "diet_missing_cell_rate": float(x_diet.isna().sum().sum() / x_diet.size),
        "metabolomics_missing_cell_rate": float(x_met.isna().sum().sum() / x_met.size),
        "diet_features_fdr_lt_0_05": int((diet_tests["fdr_bh"] < 0.05).sum()),
        "metabolite_features_fdr_lt_0_05": int((met_tests["fdr_bh"] < 0.05).sum()),
        "label_conflicts_after_stripping": int(
            (
                label_overlap["casecontrol_clean_metabolomics"]
                != label_overlap["casecontrol_clean_diet"]
            ).sum()
        ),
        "leakage_note": "No imputation is performed in this script. Modeling imputation is fit on training data only.",
    }
    (TABLES / "preprocessing_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

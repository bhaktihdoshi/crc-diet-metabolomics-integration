# CRC Diet + Metabolomics Integration

Portfolio-ready, reproducible analysis comparing late and intermediate integration strategies for colorectal cancer (CRC) classification using diet questionnaire data and NMR metabolomics.

## Project Question

Can diet features add complementary predictive and interpretive value to metabolomics for CRC case-control classification?

## Data Modalities

- **Target:** CRC status (`Colon Cancer` = 1, `Control` = 0)
- **Omics block:** NMR metabolomics features
- **Non-omics block:** diet questionnaire features

Raw data files are not tracked in Git. Place these files in `data/raw/`:

```text
Diet.xlsx
Metabolomics.xlsx
```

## Methods Summary

The current pipeline uses a leakage-safe train/validation/test design:

1. Raw data are cleaned and sample-matched.
2. Negative diet codes are recoded as missing.
3. Train/validation/test split is created after matching.
4. Feature missingness filtering, median imputation, and z-score scaling are fit on the **training set only**.
5. Validation data are used for alpha/component/threshold selection.
6. The held-out test set is used once for final evaluation.

### Main Models

**Late stacked Elastic Net**

- Fits separate diet-only and metabolomics-only Elastic Net models.
- Combines block-level CRC probabilities using an interpretable logistic stacking model.

**Intermediate Block PLS + Elastic Net**

- Fits PLS separately inside diet and metabolomics blocks.
- Combines block-level PLS scores.
- Uses Elastic Net as the final classifier on the combined scores.

## Reproducibility

From a fresh clone, create a Python environment and install the Python dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Install the R packages listed in `R_PACKAGES.md`, then run the full pipeline from the project root:

```bash
python3 scripts/01_prepare_data.py
Rscript scripts/02_run_models_leakage_safe.R
Rscript scripts/03_make_figures.R
```

Outputs are written to:

```text
data/processed/
results/tables/
results/figures/
```

## Latest Leakage-Safe Results

After correcting the earlier imputation-before-splitting issue, the held-out test results are:

| Model | Test AUC | Accuracy | Sensitivity | Specificity | F1 |
|---|---:|---:|---:|---:|---:|
| Diet-only Elastic Net | 0.551 | 0.559 | 0.400 | 0.715 | 0.473 |
| Metabolomics-only Elastic Net | 0.613 | 0.590 | 0.534 | 0.644 | 0.563 |
| Block PLS + Elastic Net | 0.612 | 0.578 | 0.495 | 0.659 | 0.536 |
| Late stacked Elastic Net | 0.617 | 0.594 | 0.571 | 0.616 | 0.582 |

Interpretation: the corrected pipeline still supports the same qualitative conclusion: metabolomics carries most of the signal, diet alone is weak, and late stacking gives a small complementary improvement. The corrected AUC values are lower than the earlier assignment draft because imputation/scaling are no longer performed before the split.

## Key Outputs

- `results/tables/model_performance.csv`
- `results/tables/model_thresholds.csv`
- `results/tables/delong_auc_tests.csv`
- `results/tables/selected_feature_overlap.csv`
- `results/figures/held_out_roc_curves.png`
- `results/figures/standalone_vs_integrated_auc.png`
- `results/figures/cross_model_feature_agreement.png`
- `results/figures/partial_spearman_correlation_network.png`
- `results/figures/exploratory_metabolite_theme_summary.png`

## Methodological Notes

- Median imputation is now fit on training data only, avoiding the earlier leakage risk from imputation before splitting.
- Elastic Net lambda tuning uses manual 5-fold CV inside the training set with fold-specific preprocessing.
- Validation data are used for model-selection choices and classification thresholds.
- Final AUCs are single held-out test estimates, not repeated-CV averages.
- No external validation or clinical covariate adjustment was available; this remains a limitation.

## Portfolio Framing

This project demonstrates:

- multimodal biomedical data integration,
- leakage-aware machine-learning workflow design,
- regularized modeling,
- intermediate vs late integration comparison,
- interpretable model reporting,
- clinically cautious biomarker interpretation.

# Methods

## Leakage-Safe Preprocessing

The earlier assignment workflow used a model-ready table where median imputation had already been applied before splitting. This project fixes that issue.

Current workflow:

1. Raw diet and metabolomics files are matched by `Sample`.
2. Negative diet questionnaire codes are recoded as missing.
3. A raw model matrix is written with missing values preserved.
4. Train/validation/test split is created.
5. Missingness filtering, median imputation, and z-score scaling are fit on the training set only.
6. During inner 5-fold CV, each fold gets its own fold-specific preprocessing fit.
7. For Block PLS + EN, each inner fold also refits PLS only on the inner-training samples before projecting the inner-validation samples.
8. Validation set is used for model selection and thresholds.
9. Held-out test set is used once for final reporting.

## Late Integration

Diet and metabolomics are modeled separately with Elastic Net logistic regression. The block-level predicted CRC probabilities are then combined by a logistic stacking model:

```text
final P(CRC) = logistic(beta0 + beta1 * diet_probability + beta2 * metabolomics_probability)
```

This tests whether diet adds complementary signal beyond metabolomics.

## Intermediate Integration

PLS is fit separately inside each block to generate supervised latent scores. Component count is selected by validation AUC, not by test AUC. The selected configuration was:

```text
diet PLS components = 1
metabolomics PLS components = 10
```

These PLS scores are concatenated and passed into an Elastic Net classifier.

During tuning, the PLS loadings are refit inside each inner CV training fold. This avoids using labels from an inner-validation fold to construct its own latent representation.

## Thresholds

Thresholds are selected on the validation set using Youden's rule:

```text
sensitivity + specificity - 1
```

AUC does not require a threshold and is reported on the held-out test set.

## Limitations

- Final AUCs are single held-out test estimates, not repeated/nested-CV averages.
- No external validation cohort was available.
- Clinical covariates such as age, BMI, stage, or medication were not available.
- Functional theme analysis is name-based and exploratory, not a formal database-mapped pathway enrichment.

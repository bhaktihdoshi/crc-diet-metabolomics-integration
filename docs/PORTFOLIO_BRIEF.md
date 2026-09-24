# Portfolio Brief

## Project Story

This project translates a multimodal cancer AI essay into a reproducible applied analysis. The essay argued that no single data modality fully captures cancer biology, and that integration strategies should preserve modality-specific structure while remaining interpretable. This repository applies that idea to CRC classification using diet questionnaire data and NMR metabolomics.

## Why These Models

Late stacked Elastic Net was used because it keeps diet and metabolomics as separate model blocks, then learns whether diet adds complementary predictive information beyond metabolomics. This is easy to explain clinically because each block produces a CRC probability before the final stack.

Block PLS + Elastic Net was used as the intermediate integration model because PLS learns supervised latent scores within each data block. Those scores reduce dimensionality while orienting the representation toward CRC separation. Elastic Net is then used on the combined latent scores for regularized classification.

## What The Results Show

The strongest standalone signal came from metabolomics. Diet alone was weak, but the late stacked model achieved the best held-out AUC, suggesting a small complementary diet contribution. Block PLS + EN performed similarly to metabolomics-only EN, supporting the idea that metabolite/lipid structure drives most of the classification signal.

## Why The Corrected Results Are Lower

An earlier assignment draft imputed data before splitting, which can leak information from validation/test samples into training. This repo fixes that: missingness filtering, median imputation, scaling, PLS projection, and inner-CV preprocessing are fit only within the correct training split. The resulting AUC values are more honest and therefore more suitable for a professional portfolio.

## Interview Talking Points

- I caught and corrected a preprocessing leakage issue.
- I compared late integration against intermediate integration rather than reporting only one model.
- I used validation-based model selection and kept the test set for final evaluation.
- I reported modest AUC honestly rather than overstating performance.
- I included XAI-style coefficient/loading summaries, feature overlap, partial correlation networks, and exploratory metabolite theme analysis to connect model performance with biological interpretation.
- I documented limitations: no external validation cohort, no clinical confounders, and exploratory pathway/theme analysis.

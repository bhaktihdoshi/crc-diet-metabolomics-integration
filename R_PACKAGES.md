# R Package Requirements

Install the R packages below before running the project:

```r
install.packages(c(
  "glmnet",
  "pROC",
  "pls",
  "ggplot2",
  "jsonlite"
))
```

The analysis was run with:

- manual 5-fold CV for Elastic Net lambda tuning,
- `glmnet` for regularized logistic regression,
- `pls` for block-level latent scores,
- `pROC` for AUC and DeLong tests,
- `ggplot2` for figures.

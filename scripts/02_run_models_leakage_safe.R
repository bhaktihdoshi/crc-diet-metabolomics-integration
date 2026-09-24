suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(pls)
  library(jsonlite)
})

set.seed(20260617)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grep("^--file=", args_all)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/02_run_models_leakage_safe.R"
project <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
processed_dir <- file.path(project, "data", "processed")
table_dir <- file.path(project, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(processed_dir, "model_matrix_raw_missing.csv"), check.names = FALSE)
df$casecontrol <- factor(df$casecontrol, levels = c("Control", "Colon Cancer"))
y <- df$y_crc
diet_cols <- grep("^diet__", names(df), value = TRUE)
met_cols <- grep("^metabolite__", names(df), value = TRUE)

stratified_split <- function(y, train_frac = 0.60, val_frac = 0.20) {
  train <- val <- test <- integer()
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    idx <- sample(idx)
    n_train <- floor(length(idx) * train_frac)
    n_val <- floor(length(idx) * val_frac)
    train <- c(train, idx[seq_len(n_train)])
    val <- c(val, idx[(n_train + 1):(n_train + n_val)])
    test <- c(test, idx[(n_train + n_val + 1):length(idx)])
  }
  list(train = sort(train), validation = sort(val), test = sort(test))
}

split <- stratified_split(y)
split_table <- rbind(
  data.frame(split = "train", Sample = df$Sample[split$train], casecontrol = df$casecontrol[split$train], y_crc = y[split$train]),
  data.frame(split = "validation", Sample = df$Sample[split$validation], casecontrol = df$casecontrol[split$validation], y_crc = y[split$validation]),
  data.frame(split = "test", Sample = df$Sample[split$test], casecontrol = df$casecontrol[split$test], y_crc = y[split$test])
)
write.csv(split_table, file.path(table_dir, "train_validation_test_split.csv"), row.names = FALSE)

auc_value <- function(y_true, prob) {
  as.numeric(pROC::auc(pROC::roc(y_true, prob, levels = c(0, 1), direction = "<", quiet = TRUE)))
}

metrics_at_threshold <- function(y_true, prob, threshold) {
  pred <- ifelse(prob >= threshold, 1, 0)
  tp <- sum(pred == 1 & y_true == 1)
  tn <- sum(pred == 0 & y_true == 0)
  fp <- sum(pred == 1 & y_true == 0)
  fn <- sum(pred == 0 & y_true == 1)
  sensitivity <- ifelse(tp + fn == 0, NA, tp / (tp + fn))
  specificity <- ifelse(tn + fp == 0, NA, tn / (tn + fp))
  precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  f1 <- ifelse(is.na(precision) || is.na(sensitivity) || precision + sensitivity == 0, NA, 2 * precision * sensitivity / (precision + sensitivity))
  data.frame(
    n = length(y_true), auc = auc_value(y_true, prob), threshold = threshold,
    accuracy = (tp + tn) / length(y_true), sensitivity = sensitivity,
    specificity = specificity, precision = precision, f1 = f1,
    tp = tp, fp = fp, tn = tn, fn = fn
  )
}

best_threshold_by_youden <- function(y_true, prob) {
  candidates <- sort(unique(prob))
  candidates <- candidates[candidates > 0 & candidates < 1]
  if (length(candidates) > 350) {
    candidates <- as.numeric(quantile(candidates, probs = seq(0.01, 0.99, length.out = 350), na.rm = TRUE))
  }
  scores <- vapply(candidates, function(t) {
    m <- metrics_at_threshold(y_true, prob, t)
    m$sensitivity + m$specificity - 1
  }, numeric(1))
  candidates[which.max(scores)]
}

prep_fit <- function(x_train, missing_threshold = 0.30) {
  missing_rate <- colMeans(is.na(x_train))
  keep <- names(missing_rate)[missing_rate <= missing_threshold]
  x_keep <- x_train[, keep, drop = FALSE]
  med <- apply(x_keep, 2, median, na.rm = TRUE)
  med[is.na(med)] <- 0
  x_imp <- x_keep
  for (j in seq_along(med)) x_imp[is.na(x_imp[, j]), j] <- med[j]
  nonconstant <- apply(x_imp, 2, function(v) length(unique(v)) > 1)
  keep2 <- names(nonconstant)[nonconstant]
  x_imp <- x_imp[, keep2, drop = FALSE]
  med <- med[keep2]
  mu <- colMeans(x_imp)
  sig <- apply(x_imp, 2, sd)
  sig[is.na(sig) | sig == 0] <- 1
  list(features = keep2, med = med, mu = mu, sig = sig, missing_rate = missing_rate)
}

prep_apply <- function(x, prep) {
  x <- x[, prep$features, drop = FALSE]
  for (j in seq_along(prep$med)) x[is.na(x[, j]), j] <- prep$med[j]
  sweep(sweep(x, 2, prep$mu, "-"), 2, prep$sig, "/")
}

make_folds <- function(y_train, k = 5) {
  folds <- integer(length(y_train))
  for (cls in sort(unique(y_train))) {
    idx <- sample(which(y_train == cls))
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  folds
}

fit_en_leakage_safe <- function(x_raw, y, train_idx, val_idx, alpha_grid = c(0, 0.25, 0.5, 0.75, 1),
                                lambda_grid = exp(seq(log(0.5), log(1e-4), length.out = 80)),
                                model_name = "Elastic Net") {
  folds <- make_folds(y[train_idx], k = 5)
  hpo <- list()
  best <- NULL
  for (alpha in alpha_grid) {
    fold_auc <- matrix(NA_real_, nrow = 5, ncol = length(lambda_grid))
    for (fold in 1:5) {
      inner_train <- train_idx[folds != fold]
      inner_val <- train_idx[folds == fold]
      prep <- prep_fit(x_raw[inner_train, , drop = FALSE])
      xtr <- prep_apply(x_raw[inner_train, , drop = FALSE], prep)
      xva <- prep_apply(x_raw[inner_val, , drop = FALSE], prep)
      fit <- glmnet(xtr, y[inner_train], family = "binomial", alpha = alpha,
                    lambda = lambda_grid, standardize = FALSE)
      prob <- predict(fit, newx = xva, s = lambda_grid, type = "response")
      fold_auc[fold, ] <- vapply(seq_along(lambda_grid), function(i) auc_value(y[inner_val], prob[, i]), numeric(1))
    }
    mean_auc <- colMeans(fold_auc, na.rm = TRUE)
    i <- which.max(mean_auc)
    prep_full <- prep_fit(x_raw[train_idx, , drop = FALSE])
    xtr_full <- prep_apply(x_raw[train_idx, , drop = FALSE], prep_full)
    xval <- prep_apply(x_raw[val_idx, , drop = FALSE], prep_full)
    fit_full <- glmnet(xtr_full, y[train_idx], family = "binomial", alpha = alpha,
                       lambda = lambda_grid, standardize = FALSE)
    p_val <- as.numeric(predict(fit_full, newx = xval, s = lambda_grid[i], type = "response"))
    row <- data.frame(
      model = model_name, alpha = alpha, lambda = lambda_grid[i],
      mean_inner_cv_auc = mean_auc[i], validation_auc = auc_value(y[val_idx], p_val)
    )
    hpo[[length(hpo) + 1]] <- row
    candidate <- list(fit = fit_full, prep = prep_full, lambda = lambda_grid[i], alpha = alpha, validation_auc = row$validation_auc)
    if (is.null(best) || candidate$validation_auc > best$validation_auc) best <- candidate
  }
  hpo_df <- do.call(rbind, hpo)
  hpo_df <- hpo_df[order(-hpo_df$validation_auc), ]
  list(best = best, hpo = hpo_df)
}

predict_en <- function(obj, x_raw) {
  x <- prep_apply(x_raw, obj$prep)
  as.numeric(predict(obj$fit, newx = x, s = obj$lambda, type = "response"))
}

project_block_pls <- function(x_diet_raw, x_met_raw, obj) {
  xd <- prep_apply(x_diet_raw, obj$prep_diet)
  xm <- prep_apply(x_met_raw, obj$prep_met)
  latent <- cbind(
    xd %*% obj$wd[, seq_len(obj$diet_components), drop = FALSE],
    xm %*% obj$wm[, seq_len(obj$metabolite_components), drop = FALSE]
  )
  colnames(latent) <- c(
    paste0("diet_PLS", seq_len(obj$diet_components)),
    paste0("met_PLS", seq_len(obj$metabolite_components))
  )
  latent
}

fit_block_pls_en_leakage_safe <- function(x_diet_raw, x_met_raw, y, train_idx, val_idx,
                                          diet_components, metabolite_components,
                                          alpha_grid = c(0, 0.25, 0.5, 0.75, 1),
                                          lambda_grid = exp(seq(log(0.5), log(1e-4), length.out = 80))) {
  folds <- make_folds(y[train_idx], k = 5)
  hpo <- list()
  best <- NULL

  for (alpha in alpha_grid) {
    fold_auc <- matrix(NA_real_, nrow = 5, ncol = length(lambda_grid))
    for (fold in 1:5) {
      inner_train <- train_idx[folds != fold]
      inner_val <- train_idx[folds == fold]

      prep_d <- prep_fit(x_diet_raw[inner_train, , drop = FALSE])
      prep_m <- prep_fit(x_met_raw[inner_train, , drop = FALSE])
      xd_tr <- prep_apply(x_diet_raw[inner_train, , drop = FALSE], prep_d)
      xd_va <- prep_apply(x_diet_raw[inner_val, , drop = FALSE], prep_d)
      xm_tr <- prep_apply(x_met_raw[inner_train, , drop = FALSE], prep_m)
      xm_va <- prep_apply(x_met_raw[inner_val, , drop = FALSE], prep_m)

      if (diet_components > ncol(xd_tr) || metabolite_components > ncol(xm_tr)) next

      pls_d <- plsr(y[inner_train] ~ xd_tr, ncomp = diet_components, center = FALSE, scale = FALSE, validation = "none")
      pls_m <- plsr(y[inner_train] ~ xm_tr, ncomp = metabolite_components, center = FALSE, scale = FALSE, validation = "none")
      wd_fold <- pls_d$loading.weights[, seq_len(diet_components), drop = FALSE]
      wm_fold <- pls_m$loading.weights[, seq_len(metabolite_components), drop = FALSE]
      latent_train <- cbind(xd_tr %*% wd_fold, xm_tr %*% wm_fold)
      latent_val <- cbind(xd_va %*% wd_fold, xm_va %*% wm_fold)
      colnames(latent_train) <- colnames(latent_val) <- c(
        paste0("diet_PLS", seq_len(diet_components)),
        paste0("met_PLS", seq_len(metabolite_components))
      )

      prep_latent <- prep_fit(latent_train)
      lat_tr <- prep_apply(latent_train, prep_latent)
      lat_va <- prep_apply(latent_val, prep_latent)
      fit <- glmnet(lat_tr, y[inner_train], family = "binomial", alpha = alpha,
                    lambda = lambda_grid, standardize = FALSE)
      prob <- predict(fit, newx = lat_va, s = lambda_grid, type = "response")
      fold_auc[fold, ] <- vapply(seq_along(lambda_grid), function(i) auc_value(y[inner_val], prob[, i]), numeric(1))
    }

    mean_auc <- colMeans(fold_auc, na.rm = TRUE)
    if (all(is.na(mean_auc))) next
    i <- which.max(mean_auc)

    prep_d_full <- prep_fit(x_diet_raw[train_idx, , drop = FALSE])
    prep_m_full <- prep_fit(x_met_raw[train_idx, , drop = FALSE])
    xd_train <- prep_apply(x_diet_raw[train_idx, , drop = FALSE], prep_d_full)
    xm_train <- prep_apply(x_met_raw[train_idx, , drop = FALSE], prep_m_full)
    xd_val <- prep_apply(x_diet_raw[val_idx, , drop = FALSE], prep_d_full)
    xm_val <- prep_apply(x_met_raw[val_idx, , drop = FALSE], prep_m_full)

    pls_d_full <- plsr(y[train_idx] ~ xd_train, ncomp = diet_components, center = FALSE, scale = FALSE, validation = "none")
    pls_m_full <- plsr(y[train_idx] ~ xm_train, ncomp = metabolite_components, center = FALSE, scale = FALSE, validation = "none")
    wd_full <- pls_d_full$loading.weights[, seq_len(diet_components), drop = FALSE]
    wm_full <- pls_m_full$loading.weights[, seq_len(metabolite_components), drop = FALSE]
    latent_train <- cbind(xd_train %*% wd_full, xm_train %*% wm_full)
    latent_val <- cbind(xd_val %*% wd_full, xm_val %*% wm_full)
    colnames(latent_train) <- colnames(latent_val) <- c(
      paste0("diet_PLS", seq_len(diet_components)),
      paste0("met_PLS", seq_len(metabolite_components))
    )

    prep_latent_full <- prep_fit(latent_train)
    lat_train <- prep_apply(latent_train, prep_latent_full)
    lat_val <- prep_apply(latent_val, prep_latent_full)
    fit_full <- glmnet(lat_train, y[train_idx], family = "binomial", alpha = alpha,
                       lambda = lambda_grid, standardize = FALSE)
    p_val <- as.numeric(predict(fit_full, newx = lat_val, s = lambda_grid[i], type = "response"))
    row <- data.frame(
      model = "Block PLS + Elastic Net", alpha = alpha, lambda = lambda_grid[i],
      mean_inner_cv_auc = mean_auc[i], validation_auc = auc_value(y[val_idx], p_val),
      diet_components = diet_components, metabolite_components = metabolite_components
    )
    hpo[[length(hpo) + 1]] <- row
    candidate <- list(
      fit = fit_full, lambda = lambda_grid[i], alpha = alpha,
      prep_diet = prep_d_full, prep_met = prep_m_full, prep_latent = prep_latent_full,
      wd = wd_full, wm = wm_full,
      diet_components = diet_components, metabolite_components = metabolite_components,
      validation_auc = row$validation_auc
    )
    if (is.null(best) || candidate$validation_auc > best$validation_auc) best <- candidate
  }

  hpo_df <- do.call(rbind, hpo)
  hpo_df <- hpo_df[order(-hpo_df$validation_auc), ]
  list(best = best, hpo = hpo_df)
}

predict_block_pls_en <- function(obj, x_diet_raw, x_met_raw) {
  latent <- project_block_pls(x_diet_raw, x_met_raw, obj)
  latent <- prep_apply(latent, obj$prep_latent)
  as.numeric(predict(obj$fit, newx = latent, s = obj$lambda, type = "response"))
}

x_diet <- as.matrix(df[, diet_cols])
x_met <- as.matrix(df[, met_cols])
storage.mode(x_diet) <- "double"
storage.mode(x_met) <- "double"

late_diet <- fit_en_leakage_safe(x_diet, y, split$train, split$validation, model_name = "Diet-only Elastic Net")
late_met <- fit_en_leakage_safe(x_met, y, split$train, split$validation, model_name = "Metabolomics-only Elastic Net")

p_diet_val <- predict_en(late_diet$best, x_diet[split$validation, , drop = FALSE])
p_diet_test <- predict_en(late_diet$best, x_diet[split$test, , drop = FALSE])
p_met_val <- predict_en(late_met$best, x_met[split$validation, , drop = FALSE])
p_met_test <- predict_en(late_met$best, x_met[split$test, , drop = FALSE])

stack_train <- data.frame(y = y[split$validation], diet_prob = p_diet_val, metabolomics_prob = p_met_val)
stack_fit <- glm(y ~ diet_prob + metabolomics_prob, data = stack_train, family = binomial())
p_stack_val <- as.numeric(predict(stack_fit, newdata = stack_train, type = "response"))
p_stack_test <- as.numeric(predict(stack_fit, newdata = data.frame(diet_prob = p_diet_test, metabolomics_prob = p_met_test), type = "response"))

# Intermediate Block PLS + EN.
train_prep_diet <- prep_fit(x_diet[split$train, , drop = FALSE])
train_prep_met <- prep_fit(x_met[split$train, , drop = FALSE])
max_d_pls <- min(8, length(train_prep_diet$features), length(split$train) - 1)
max_m_pls <- min(20, length(train_prep_met$features), length(split$train) - 1)
pls_grid <- expand.grid(
  diet_components = c(1, 2, 4, 8),
  metabolite_components = c(1, 2, 5, 10, 20)
)
pls_grid <- pls_grid[pls_grid$diet_components <= max_d_pls & pls_grid$metabolite_components <= max_m_pls, ]

block_results <- list()
block_hpo <- list()
block_preds <- list()
for (i in seq_len(nrow(pls_grid))) {
  dp <- pls_grid$diet_components[i]
  mp <- pls_grid$metabolite_components[i]
  fit <- fit_block_pls_en_leakage_safe(x_diet, x_met, y, split$train, split$validation, dp, mp)
  p_val <- predict_block_pls_en(fit$best, x_diet[split$validation, , drop = FALSE], x_met[split$validation, , drop = FALSE])
  p_test <- predict_block_pls_en(fit$best, x_diet[split$test, , drop = FALSE], x_met[split$test, , drop = FALSE])
  thr <- best_threshold_by_youden(y[split$validation], p_val)
  key <- paste("diet", dp, "met", mp, sep = "_")
  block_results[[key]] <- data.frame(
    key = key, diet_components = dp, metabolite_components = mp,
    validation_auc = auc_value(y[split$validation], p_val),
    test_auc = auc_value(y[split$test], p_test)
  )
  block_hpo[[key]] <- fit$hpo
  block_preds[[key]] <- list(fit = fit, p_val = p_val, p_test = p_test, threshold = thr, dp = dp, mp = mp)
}
block_summary <- do.call(rbind, block_results)
block_summary <- block_summary[order(-block_summary$validation_auc), ]
best_key <- block_summary$key[1]
best_block <- block_preds[[best_key]]
p_block_val <- best_block$p_val
p_block_test <- best_block$p_test

thr_diet <- best_threshold_by_youden(y[split$validation], p_diet_val)
thr_met <- best_threshold_by_youden(y[split$validation], p_met_val)
thr_stack <- best_threshold_by_youden(y[split$validation], p_stack_val)
thr_block <- best_threshold_by_youden(y[split$validation], p_block_val)

perf <- rbind(
  cbind(model = "Diet-only Elastic Net", split = "validation", metrics_at_threshold(y[split$validation], p_diet_val, thr_diet)),
  cbind(model = "Diet-only Elastic Net", split = "test", metrics_at_threshold(y[split$test], p_diet_test, thr_diet)),
  cbind(model = "Metabolomics-only Elastic Net", split = "validation", metrics_at_threshold(y[split$validation], p_met_val, thr_met)),
  cbind(model = "Metabolomics-only Elastic Net", split = "test", metrics_at_threshold(y[split$test], p_met_test, thr_met)),
  cbind(model = "Late stacked Elastic Net", split = "validation", metrics_at_threshold(y[split$validation], p_stack_val, thr_stack)),
  cbind(model = "Late stacked Elastic Net", split = "test", metrics_at_threshold(y[split$test], p_stack_test, thr_stack)),
  cbind(model = "Block PLS + Elastic Net", split = "validation", metrics_at_threshold(y[split$validation], p_block_val, thr_block)),
  cbind(model = "Block PLS + Elastic Net", split = "test", metrics_at_threshold(y[split$test], p_block_test, thr_block))
)
write.csv(perf, file.path(table_dir, "model_performance.csv"), row.names = FALSE)
write.csv(data.frame(
  model = c("Diet-only Elastic Net", "Metabolomics-only Elastic Net", "Late stacked Elastic Net", "Block PLS + Elastic Net"),
  threshold = c(thr_diet, thr_met, thr_stack, thr_block),
  threshold_source = "validation_youden"
), file.path(table_dir, "model_thresholds.csv"), row.names = FALSE)
write.csv(rbind(late_diet$hpo, late_met$hpo), file.path(table_dir, "elastic_net_hpo.csv"), row.names = FALSE)
write.csv(do.call(rbind, block_hpo), file.path(table_dir, "block_pls_hpo.csv"), row.names = FALSE)
write.csv(block_summary, file.path(table_dir, "block_pls_component_search.csv"), row.names = FALSE)

pred <- rbind(
  data.frame(split = "validation", Sample = df$Sample[split$validation], y_crc = y[split$validation],
             diet_prob = p_diet_val, metabolomics_prob = p_met_val, late_stack_prob = p_stack_val, block_pls_en_prob = p_block_val),
  data.frame(split = "test", Sample = df$Sample[split$test], y_crc = y[split$test],
             diet_prob = p_diet_test, metabolomics_prob = p_met_test, late_stack_prob = p_stack_test, block_pls_en_prob = p_block_test)
)
write.csv(pred, file.path(table_dir, "model_predictions.csv"), row.names = FALSE)

delong_row <- function(reference, comparator, reference_prob, comparator_prob) {
  rt <- pROC::roc(y[split$test], reference_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
  ct <- pROC::roc(y[split$test], comparator_prob, levels = c(0, 1), direction = "<", quiet = TRUE)
  test <- pROC::roc.test(rt, ct, method = "delong", paired = TRUE)
  data.frame(
    reference_model = reference, comparator_model = comparator,
    reference_auc = as.numeric(pROC::auc(rt)), comparator_auc = as.numeric(pROC::auc(ct)),
    auc_difference = as.numeric(pROC::auc(rt)) - as.numeric(pROC::auc(ct)),
    p_value = as.numeric(test$p.value)
  )
}
delong <- rbind(
  delong_row("Late stacked Elastic Net", "Diet-only Elastic Net", p_stack_test, p_diet_test),
  delong_row("Late stacked Elastic Net", "Metabolomics-only Elastic Net", p_stack_test, p_met_test),
  delong_row("Late stacked Elastic Net", "Block PLS + Elastic Net", p_stack_test, p_block_test)
)
write.csv(delong, file.path(table_dir, "delong_auc_tests.csv"), row.names = FALSE)

# XAI: direct standardized coefficients for late block models.
coef_table <- function(model_obj, prefix) {
  b <- as.matrix(coef(model_obj$best$fit, s = model_obj$best$lambda))
  values <- as.numeric(b[-1, , drop = FALSE])
  features <- rownames(b)[-1]
  data.frame(feature = features, coefficient = values, block = prefix)
}
late_coef <- rbind(coef_table(late_diet, "diet"), coef_table(late_met, "metabolomics"))
late_coef <- late_coef[late_coef$coefficient != 0, ]
late_coef$abs_coefficient <- abs(late_coef$coefficient)
write.csv(late_coef[order(-late_coef$abs_coefficient), ], file.path(table_dir, "late_model_coefficients.csv"), row.names = FALSE)

# Back-project Block PLS + EN component coefficients to original feature importance.
block_coef <- as.matrix(coef(best_block$fit$best$fit, s = best_block$fit$best$lambda))
component_coef <- as.numeric(block_coef[-1, , drop = FALSE])
component_names <- rownames(block_coef)[-1]
component_df <- data.frame(component = component_names, coefficient = component_coef)
write.csv(component_df, file.path(table_dir, "block_pls_component_coefficients.csv"), row.names = FALSE)
dp <- best_block$fit$best$diet_components
mp <- best_block$fit$best$metabolite_components
coef_lookup <- setNames(rep(0, dp + mp), c(paste0("diet_PLS", seq_len(dp)), paste0("met_PLS", seq_len(mp))))
coef_lookup[component_names] <- component_coef
diet_component_coef <- coef_lookup[paste0("diet_PLS", seq_len(dp))]
met_component_coef <- coef_lookup[paste0("met_PLS", seq_len(mp))]
diet_importance <- rowSums(abs(best_block$fit$best$wd[, seq_len(dp), drop = FALSE] %*% diag(diet_component_coef, nrow = dp)))
met_importance <- rowSums(abs(best_block$fit$best$wm[, seq_len(mp), drop = FALSE] %*% diag(met_component_coef, nrow = mp)))
block_importance <- rbind(
  data.frame(feature = best_block$fit$best$prep_diet$features, block = "diet", importance = diet_importance),
  data.frame(feature = best_block$fit$best$prep_met$features, block = "metabolomics", importance = met_importance)
)
block_importance <- block_importance[order(-block_importance$importance), ]
write.csv(block_importance, file.path(table_dir, "block_pls_backprojected_feature_importance.csv"), row.names = FALSE)

top_late <- head(late_coef[order(-late_coef$abs_coefficient), "feature"], 25)
top_block <- head(block_importance$feature, 25)
overlap <- data.frame(feature = sort(unique(c(top_late, top_block))))
overlap$late_model <- overlap$feature %in% top_late
overlap$block_pls_en <- overlap$feature %in% top_block
overlap$category <- ifelse(overlap$late_model & overlap$block_pls_en, "Common",
                           ifelse(overlap$late_model, "Late only", "Block PLS + EN only"))
write.csv(overlap, file.path(table_dir, "selected_feature_overlap.csv"), row.names = FALSE)

summary <- list(
  leakage_fix = "Modeling imputation, feature filtering, scaling, PLS scores, and inner CV preprocessing are fit within the appropriate training split only.",
  selected_block_pls_components = list(diet = best_block$fit$best$diet_components, metabolomics = best_block$fit$best$metabolite_components),
  test_auc = setNames(perf$auc[perf$split == "test"], perf$model[perf$split == "test"]),
  thresholds = c(diet = thr_diet, metabolomics = thr_met, late_stack = thr_stack, block_pls_en = thr_block)
)
writeLines(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), file.path(table_dir, "modeling_summary.json"))
print(perf[perf$split == "test", c("model", "auc", "accuracy", "sensitivity", "specificity", "precision", "f1")])

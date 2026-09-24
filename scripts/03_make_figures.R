suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grep("^--file=", args_all)]
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/03_make_figures.R"
project <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
table_dir <- file.path(project, "results", "tables")
figure_dir <- file.path(project, "results", "figures")
processed_dir <- file.path(project, "data", "processed")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

COL <- c(navy = "#243B53", teal = "#1B998B", coral = "#E76F51", gold = "#F4A261",
         gray = "#7A869A", light = "#F5F7FA", ink = "#111827")
C <- function(x) unname(COL[x])

theme_portfolio <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = C("navy"), margin = margin(b = 8)),
      plot.subtitle = element_text(hjust = 0.5, color = C("gray"), margin = margin(b = 8)),
      panel.grid = element_blank(),
      axis.line = element_line(color = C("ink"), linewidth = 0.35),
      axis.ticks = element_line(color = C("ink"), linewidth = 0.35),
      axis.text = element_text(color = C("ink")),
      axis.title = element_text(color = C("ink")),
      legend.position = "bottom",
      legend.background = element_blank(),
      legend.box.background = element_rect(fill = "white", color = C("navy"), linewidth = 0.35),
      plot.margin = margin(18, 28, 18, 22)
    )
}

save_plot <- function(name, plot, width, height) {
  ggsave(file.path(figure_dir, name), plot, width = width, height = height, dpi = 320, bg = "white")
}

clean_label <- function(x) {
  x <- sub("^diet__", "", x)
  x <- sub("^metabolite__", "", x)
  gsub("_", " ", x)
}

wrap_label <- function(x, width = 28) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n", FUN.VALUE = character(1))
}

pred <- read.csv(file.path(table_dir, "model_predictions.csv"), check.names = FALSE)
perf <- read.csv(file.path(table_dir, "model_performance.csv"), check.names = FALSE)
test_perf <- perf[perf$split == "test", ]
test_pred <- pred[pred$split == "test", ]

roc_df <- function(model, y, p) {
  r <- pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)
  data.frame(model = model, fpr = 1 - r$specificities, sensitivity = r$sensitivities, auc = as.numeric(pROC::auc(r)))
}
roc_plot <- rbind(
  roc_df("Late stacked EN", test_pred$y_crc, test_pred$late_stack_prob),
  roc_df("Block PLS + EN", test_pred$y_crc, test_pred$block_pls_en_prob)
)
late_roc_label <- sprintf("Late stacked EN (AUC = %.3f)", test_perf$auc[test_perf$model == "Late stacked Elastic Net"])
block_roc_label <- sprintf("Block PLS + EN (AUC = %.3f)", test_perf$auc[test_perf$model == "Block PLS + Elastic Net"])
roc_plot$label <- ifelse(roc_plot$model == "Late stacked EN",
                         late_roc_label,
                         block_roc_label)
roc_colors <- c(C("navy"), C("coral"))
names(roc_colors) <- c(late_roc_label, block_roc_label)
p <- ggplot(roc_plot, aes(fpr, sensitivity, color = label)) +
  geom_line(linewidth = 1.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = C("gray")) +
  coord_equal() +
  scale_color_manual(values = roc_colors, name = "Model") +
  labs(title = "Held-Out ROC Curves", x = "False positive rate", y = "Sensitivity") +
  theme_portfolio(12)
save_plot("held_out_roc_curves.png", p, 7.6, 6.2)

auc_df <- test_perf[test_perf$model %in% c("Diet-only Elastic Net", "Metabolomics-only Elastic Net", "Late stacked Elastic Net", "Block PLS + Elastic Net"), ]
auc_df$model <- factor(auc_df$model, levels = c("Diet-only Elastic Net", "Metabolomics-only Elastic Net", "Block PLS + Elastic Net", "Late stacked Elastic Net"))
auc_df$model_label <- factor(
  c("Diet-only\nEN", "Metabolomics-only\nEN", "Block PLS +\nEN", "Late stacked\nEN")[match(as.character(auc_df$model), levels(auc_df$model))],
  levels = c("Diet-only\nEN", "Metabolomics-only\nEN", "Block PLS +\nEN", "Late stacked\nEN")
)
auc_df$role <- ifelse(grepl("only", auc_df$model), "Standalone", "Integrated")
p <- ggplot(auc_df, aes(model_label, auc, fill = role)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", auc)), vjust = -0.35, fontface = "bold", color = C("navy"), size = 4.2) +
  scale_fill_manual(values = c("Standalone" = C("teal"), "Integrated" = C("coral"))) +
  coord_cartesian(ylim = c(0.50, max(auc_df$auc) + 0.04), clip = "off") +
  labs(title = "Standalone vs Integrated Model AUC", x = "Model", y = "Held-out test AUC", fill = NULL) +
  theme_portfolio(12) +
  theme(axis.text.x = element_text(size = 10.5, lineheight = 0.95), legend.position = "none")
save_plot("standalone_vs_integrated_auc.png", p, 8.4, 5.2)

make_cm <- function(model_name, prob_col, threshold, file, title) {
  y <- test_pred$y_crc
  pred_class <- ifelse(test_pred[[prob_col]] >= threshold, 1, 0)
  cm <- data.frame(
    actual = factor(c("CRC", "CRC", "Control", "Control"), levels = c("CRC", "Control")),
    predicted = factor(c("CRC", "Control", "CRC", "Control"), levels = c("CRC", "Control")),
    n = c(sum(y == 1 & pred_class == 1), sum(y == 1 & pred_class == 0),
          sum(y == 0 & pred_class == 1), sum(y == 0 & pred_class == 0)),
    cell = c("True positive", "False negative", "False positive", "True negative")
  )
  cm$label <- paste0(cm$n, "\n", cm$cell)
  cm$text_color <- ifelse(cm$n >= median(cm$n), "white", C("ink"))
  p <- ggplot(cm, aes(predicted, actual, fill = n)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = label, color = text_color), size = 4.8, fontface = "bold", lineheight = 0.95) +
    scale_color_identity() +
    scale_fill_gradient(low = "#D9F0EC", high = C("navy")) +
    coord_equal() +
    labs(title = title, subtitle = sprintf("%s; validation-selected threshold = %.3f", model_name, threshold),
         x = "Predicted class", y = "Actual class") +
    theme_portfolio(12) +
    theme(legend.position = "none", axis.ticks = element_blank(), panel.grid = element_blank())
  save_plot(file, p, 6.4, 5.6)
}
thresholds <- read.csv(file.path(table_dir, "model_thresholds.csv"), check.names = FALSE)
make_cm("Late stacked EN", "late_stack_prob", thresholds$threshold[thresholds$model == "Late stacked Elastic Net"],
        "final_late_stacked_en_confusion_matrix.png", "Final Model Confusion Matrix")
make_cm("Block PLS + EN", "block_pls_en_prob", thresholds$threshold[thresholds$model == "Block PLS + Elastic Net"],
        "intermediate_block_pls_en_confusion_matrix.png", "Intermediate Model Confusion Matrix")

overlap <- read.csv(file.path(table_dir, "selected_feature_overlap.csv"), check.names = FALSE)
n_common <- sum(overlap$category == "Common")
n_late <- sum(overlap$category == "Late only")
n_block <- sum(overlap$category == "Block PLS + EN only")
circle_points <- function(cx, cy, r, n = 360) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}
p <- ggplot() +
  geom_polygon(data = circle_points(-0.55, 0, 1.05), aes(x, y), fill = C("navy"), alpha = 0.30, color = C("navy"), linewidth = 1.1) +
  geom_polygon(data = circle_points(0.55, 0, 1.05), aes(x, y), fill = C("teal"), alpha = 0.34, color = C("teal"), linewidth = 1.1) +
  annotate("text", x = -1.18, y = 0.05, label = n_late, size = 10.5, fontface = "bold", color = C("navy")) +
  annotate("text", x = 0, y = 0.05, label = n_common, size = 11.5, fontface = "bold", color = C("ink")) +
  annotate("text", x = 1.18, y = 0.05, label = n_block, size = 10.5, fontface = "bold", color = C("teal")) +
  annotate("text", x = -0.92, y = 1.18, label = "Late stacked EN", size = 5.0, fontface = "bold", color = C("navy")) +
  annotate("text", x = 0.92, y = 1.18, label = "Block PLS + EN", size = 5.0, fontface = "bold", color = C("teal")) +
  labs(title = "Cross-Model Feature Agreement") +
  coord_equal(xlim = c(-2, 2), ylim = c(-1.25, 1.42), clip = "off") +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = C("navy"), margin = margin(b = 8)))
save_plot("cross_model_feature_agreement.png", p, 8.8, 5.4)

# XAI bars.
late_coef <- read.csv(file.path(table_dir, "late_model_coefficients.csv"), check.names = FALSE)
block_imp <- read.csv(file.path(table_dir, "block_pls_backprojected_feature_importance.csv"), check.names = FALSE)
xai <- rbind(
  data.frame(feature = late_coef$feature, model = "Late stacked EN", importance = late_coef$abs_coefficient),
  data.frame(feature = block_imp$feature, model = "Block PLS + EN", importance = block_imp$importance)
)
feature_rank <- aggregate(importance ~ feature, xai, max)
feature_rank <- feature_rank[order(-feature_rank$importance), ]
top_features <- head(feature_rank$feature, 12)
xai <- xai[xai$feature %in% top_features, ]
ord <- aggregate(importance ~ feature, xai, max)
ord <- ord[order(ord$importance), ]
xai$feature_label <- factor(wrap_label(clean_label(xai$feature), 30), levels = wrap_label(clean_label(ord$feature), 30))
p <- ggplot(xai, aes(importance, feature_label, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.62) +
  scale_fill_manual(values = c("Late stacked EN" = C("navy"), "Block PLS + EN" = C("teal"))) +
  labs(title = "XAI Importance: Top 12 Lipid-Core Variables", x = "Relative absolute importance", y = "Variable", fill = "Model") +
  theme_portfolio(12) +
  theme(axis.text.y = element_text(size = 11.2, lineheight = 0.92))
save_plot("xai_top_feature_importance.png", p, 10.4, 5.9)

# Partial Spearman network and heatmap after adjusting for CRC.
raw <- read.csv(file.path(processed_dir, "model_matrix_raw_missing.csv"), check.names = FALSE)
partial_edges <- list()
diet_cols <- grep("^diet__", names(raw), value = TRUE)
met_cols <- grep("^metabolite__", names(raw), value = TRUE)
for (d in diet_cols) {
  for (m in met_cols) {
    rd <- rank(raw[[d]], na.last = "keep", ties.method = "average")
    rm <- rank(raw[[m]], na.last = "keep", ties.method = "average")
    ok <- complete.cases(rd, rm, raw$y_crc)
    if (sum(ok) < 20) next
    d_res <- residuals(lm(rd[ok] ~ raw$y_crc[ok]))
    m_res <- residuals(lm(rm[ok] ~ raw$y_crc[ok]))
    ct <- suppressWarnings(cor.test(d_res, m_res, method = "pearson"))
    partial_edges[[length(partial_edges) + 1]] <- data.frame(
      diet_feature = d, metabolite_feature = m,
      partial_spearman_rho = unname(ct$estimate), p_value = ct$p.value
    )
  }
}
partial_edges <- do.call(rbind, partial_edges)
partial_edges$fdr_bh <- p.adjust(partial_edges$p_value, method = "BH")
write.csv(partial_edges, file.path(table_dir, "partial_spearman_diet_metabolite_edges.csv"), row.names = FALSE)
edge_plot <- partial_edges[partial_edges$fdr_bh < 0.10 & abs(partial_edges$partial_spearman_rho) >= 0.045, ]
edge_plot <- head(edge_plot[order(edge_plot$fdr_bh), ], 22)
if (nrow(edge_plot) > 0) {
  diet_nodes <- sort(unique(edge_plot$diet_feature))
  met_nodes <- sort(unique(edge_plot$metabolite_feature))
  nodes <- rbind(
    data.frame(name = diet_nodes, type = "Diet", x = 0, y = seq(0.15, 0.85, length.out = length(diet_nodes))),
    data.frame(name = met_nodes, type = "Metabolite", x = 1, y = seq(0.02, 0.98, length.out = length(met_nodes)))
  )
  e <- merge(edge_plot, nodes[, c("name", "x", "y")], by.x = "diet_feature", by.y = "name")
  e <- merge(e, nodes[, c("name", "x", "y")], by.x = "metabolite_feature", by.y = "name", suffixes = c("_diet", "_met"))
  p <- ggplot() +
    geom_curve(data = e, aes(x = x_diet + 0.04, y = y_diet, xend = x_met - 0.04, yend = y_met,
                             color = partial_spearman_rho, linewidth = abs(partial_spearman_rho)),
               curvature = 0.08, alpha = 0.85) +
    geom_point(data = nodes, aes(x, y, fill = type), shape = 21, size = 5.2, color = "white") +
    geom_text(data = nodes, aes(x = ifelse(type == "Diet", x - 0.06, x + 0.06), y = y, label = wrap_label(clean_label(name), 23),
                                hjust = ifelse(type == "Diet", 1, 0)), size = 3.0) +
    annotate("text", x = 0, y = 1.08, label = "Diet", color = C("navy"), fontface = "bold", size = 4.2) +
    annotate("text", x = 1, y = 1.08, label = "Metabolomics", color = C("navy"), fontface = "bold", size = 4.2) +
    scale_fill_manual(values = c("Diet" = C("gold"), "Metabolite" = C("teal")), name = "Node type") +
    scale_color_gradient2(low = C("navy"), mid = "#D1D5DB", high = C("coral"), midpoint = 0, name = "Partial Spearman rho") +
    scale_linewidth(range = c(0.35, 1.9), name = "Correlation strength") +
    xlim(-0.48, 1.58) + ylim(-0.08, 1.15) +
    labs(title = "Partial Spearman Correlation Network") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", color = C("navy")), legend.position = "right")
  save_plot("partial_spearman_correlation_network.png", p, 11.2, 7.0)
  heat <- edge_plot
  heat$diet_label <- factor(clean_label(heat$diet_feature), levels = rev(clean_label(diet_nodes)))
  heat$met_label <- factor(wrap_label(clean_label(heat$metabolite_feature), 24), levels = wrap_label(clean_label(met_nodes), 24))
  heat$text_color <- ifelse(heat$partial_spearman_rho <= -0.052, "white", C("ink"))
  p <- ggplot(heat, aes(met_label, diet_label, fill = partial_spearman_rho)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.3f", partial_spearman_rho), color = text_color), fontface = "bold", size = 3.0) +
    scale_color_identity() +
    scale_fill_gradient2(low = C("navy"), mid = "white", high = C("coral"), midpoint = 0, name = "Partial Spearman rho") +
    labs(title = "Partial Spearman Heatmap", x = "Metabolite", y = "Diet feature") +
    theme_portfolio(9.2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", legend.box.background = element_blank())
  save_plot("partial_spearman_heatmap.png", p, 11.8, 5.0)
}

theme_for <- function(x) {
  z <- clean_label(x)
  if (grepl("LDL|IDL", z, ignore.case = TRUE)) return("LDL/IDL cholesterol transport")
  if (grepl("HDL", z, ignore.case = TRUE)) return("HDL metabolism")
  if (grepl("Cholesterol", z, ignore.case = TRUE)) return("Cholesterol metabolism")
  if (grepl("Phospholipid", z, ignore.case = TRUE)) return("Phospholipid metabolism")
  "Other lipid/metabolite signal"
}
met_features <- grep("^metabolite__", names(raw), value = TRUE)
selected_met <- unique(c(late_coef$feature[late_coef$block == "metabolomics"], head(block_imp$feature[block_imp$block == "metabolomics"], 25)))
theme_table <- data.frame(theme = vapply(met_features, theme_for, character(1)), feature = met_features)
selected_table <- data.frame(theme = vapply(selected_met, theme_for, character(1)), feature = selected_met)
theme_counts <- aggregate(feature ~ theme, theme_table, length)
names(theme_counts)[2] <- "background_count"
selected_counts <- aggregate(feature ~ theme, selected_table, length)
names(selected_counts)[2] <- "selected_count"
enrich <- merge(theme_counts, selected_counts, by = "theme", all.x = TRUE)
enrich$selected_count[is.na(enrich$selected_count)] <- 0
total_selected <- length(unique(selected_met))
total_background <- length(unique(met_features))
enrich$p_value <- vapply(seq_len(nrow(enrich)), function(i) {
  a <- enrich$selected_count[i]
  b <- total_selected - a
  c <- enrich$background_count[i] - a
  d <- total_background - enrich$background_count[i] - b
  fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
}, numeric(1))
enrich$fdr_bh <- p.adjust(enrich$p_value, method = "BH")
write.csv(enrich, file.path(table_dir, "exploratory_metabolite_theme_enrichment.csv"), row.names = FALSE)
plot_enrich <- enrich[enrich$selected_count > 0, ]
plot_enrich$theme <- factor(plot_enrich$theme, levels = plot_enrich$theme[order(plot_enrich$selected_count)])
plot_enrich$label <- paste0(plot_enrich$selected_count, "/", plot_enrich$background_count, " selected")
p <- ggplot(plot_enrich, aes(selected_count, theme)) +
  geom_col(width = 0.52, fill = C("navy")) +
  geom_text(aes(x = selected_count + 0.25, label = selected_count), hjust = 0, size = 4.5, fontface = "bold", color = C("navy")) +
  geom_text(aes(x = max(plot_enrich$selected_count) + 2.0, label = label), hjust = 0, color = C("gray"), size = 4.0) +
  annotate("text", x = max(plot_enrich$selected_count) + 2.0, y = length(levels(plot_enrich$theme)) + 0.42,
           label = "Selected / measured", hjust = 0, size = 3.7, fontface = "bold", color = C("navy")) +
  scale_x_continuous(limits = c(0, max(plot_enrich$selected_count) + 17), breaks = seq(0, max(plot_enrich$selected_count), by = 5), expand = expansion(mult = c(0, 0))) +
  labs(title = "Exploratory Metabolite Theme Summary",
       subtitle = ifelse(any(plot_enrich$fdr_bh < 0.05), "Name-based themes; FDR shown in table", "Name-based themes; no FDR-significant enrichment"),
       x = "Selected metabolite features", y = "Functional theme",
       caption = "FDR calculated with Benjamini-Hochberg correction across name-based metabolite themes.") +
  coord_cartesian(clip = "off") +
  theme_portfolio(12.5) +
  theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, color = C("gray"), size = 10.5, margin = margin(t = 8)))
save_plot("exploratory_metabolite_theme_summary.png", p, 12.2, 5.2)

cat("Figures written to:", figure_dir, "\n")

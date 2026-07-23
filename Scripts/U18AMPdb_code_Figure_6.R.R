# ============================================================
#ROC curves for Fig.6A-B；Fig.S14A-D
# ============================================================

library(tidyverse)
library(pROC)
pred_files <- tibble(
  Dataset = c(
    "pre-NEC-Species",
    "pre-NEC-AMPs"
  ),
  Best_topn = c(
    "all",
    "4096"
  ),
  File = c(
    "C:/Users/taoch/Desktop/632751_clusters/Fig_7/new_figure_for_NEC/taxa/Pre-NEC/taxa_RF_predictions_gut.csv",
    "C:/Users/taoch/Desktop/632751_clusters/Fig_7/new_figure_for_NEC/amp/pre-NEC/amp_RF_predictions_gut_simple.csv"
  )
)

out_dir <- "C:/Users/taoch/Desktop/632751_clusters/Fig_7/new_figure_for_NEC/pre-NEC_taxa_and_amp"
dataset_cols <- c(
  "pre-NEC-Species"  = "#5F5A8E",
  "pre-NEC-AMPs" = "#CB9927"
)

calculate_mean_roc <- function(pred_raw_file, dataset_name, best_topn) {
  
  cat("\n============================================================\n")
  cat("Process the dataset: ", dataset_name, "\n", sep = "")
  cat("best_topn: ", best_topn, "\n", sep = "")
  cat("file: ", pred_raw_file, "\n", sep = "")
  cat("============================================================\n")
  
  pred_raw <- read_csv(
    pred_raw_file,
    show_col_types = FALSE
  )
  
  if (!"Repeat_ID" %in% colnames(pred_raw) && "repeat" %in% colnames(pred_raw)) {
    pred_raw <- pred_raw %>%
      rename(Repeat_ID = `repeat`)
  }
  
  required_cols <- c("Repeat_ID", "top_n", "Run", "y_true", "y_pred")
  missing_cols <- setdiff(required_cols, colnames(pred_raw))
  
  if (length(missing_cols) > 0) {
    stop(
      dataset_name,
      " prediction 文件缺少以下列: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  pred_raw <- pred_raw %>%
    mutate(
      Dataset = dataset_name,
      Repeat_ID = as.integer(Repeat_ID),
      top_n = as.character(top_n),
      y_true = as.character(y_true),
      y_pred = as.numeric(y_pred)
    )
  
  repeat_ids <- sort(unique(pred_raw$Repeat_ID))
  n_repeats <- length(repeat_ids)

  best_pred <- pred_raw %>%
    filter(top_n == best_topn)
  
  if (nrow(best_pred) == 0) {
    stop(dataset_name, ": best_topn = ", best_topn, " pred_raw not exist")
  }
  
  fpr_grid <- seq(0, 1, length.out = 101)
  
  roc_grid_list <- list()
  roc_auc_list <- list()
  
  for (r in repeat_ids) {
    
    df_r <- best_pred %>%
      filter(Repeat_ID == r)
    
    if (nrow(df_r) == 0) {
      warning(dataset_name, ": Repeat_ID = ", r, " jump")
      next
    }
    
    if (length(unique(df_r$y_true)) < 2) {
      warning(dataset_name, ": Repeat_ID = ", r, " jump")
      next
    }
    
    y_r <- factor(
      df_r$y_true,
      levels = c("Health", "Disease")
    )
    
    pred_r <- df_r$y_pred
    
    roc_r <- roc(
      response = y_r,
      predictor = pred_r,
      levels = c("Health", "Disease"),
      direction = "<",
      quiet = TRUE
    )
    
    auc_r <- as.numeric(auc(roc_r))
    
    roc_auc_list[[as.character(r)]] <- tibble(
      Dataset = dataset_name,
      Repeat_ID = r,
      top_n = best_topn,
      AUC = auc_r
    )
    
    tpr_r <- coords(
      roc_r,
      x = 1 - fpr_grid,
      input = "specificity",
      ret = "sensitivity",
      transpose = FALSE
    )
    
    if (is.data.frame(tpr_r)) {
      tpr_r <- tpr_r$sensitivity
    } else {
      tpr_r <- as.numeric(tpr_r)
    }
    
    roc_grid_list[[as.character(r)]] <- tibble(
      Dataset = dataset_name,
      Repeat_ID = r,
      top_n = best_topn,
      FPR = fpr_grid,
      TPR = tpr_r
    )
  }
  
  roc_grid_raw <- bind_rows(roc_grid_list)
  roc_auc_raw <- bind_rows(roc_auc_list)
  
  if (nrow(roc_grid_raw) == 0) {
    stop(dataset_name, "no ROC result")
  }
  
  roc_grid_summary <- roc_grid_raw %>%
    group_by(Dataset, top_n, FPR) %>%
    summarise(
      mean_TPR = mean(TPR, na.rm = TRUE),
      sd_TPR = sd(TPR, na.rm = TRUE),
      lower_TPR = pmax(mean_TPR - sd_TPR, 0),
      upper_TPR = pmin(mean_TPR + sd_TPR, 1),
      .groups = "drop"
    )
  
  auc_summary <- roc_auc_raw %>%
    summarise(
      Dataset = dataset_name,
      top_n = best_topn,
      n_repeats = n(),
      mean_AUC = mean(AUC, na.rm = TRUE),
      sd_AUC = sd(AUC, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Mean AUROC =", round(auc_summary$mean_AUC, 4), "\n")
  cat("SD AUROC   =", round(auc_summary$sd_AUC, 4), "\n")
  
  return(
    list(
      roc_grid_raw = roc_grid_raw,
      roc_grid_summary = roc_grid_summary,
      roc_auc_raw = roc_auc_raw,
      auc_summary = auc_summary
    )
  )
}

all_results <- pmap(
  list(
    pred_raw_file = pred_files$File,
    dataset_name = pred_files$Dataset,
    best_topn = pred_files$Best_topn
  ),
  calculate_mean_roc
)

roc_grid_raw_all <- map_dfr(all_results, "roc_grid_raw")
roc_grid_summary_all <- map_dfr(all_results, "roc_grid_summary")
roc_auc_raw_all <- map_dfr(all_results, "roc_auc_raw")
auc_summary_all <- map_dfr(all_results, "auc_summary")

auc_label_df <- auc_summary_all %>%
  mutate(
    Dataset_label = paste0(
      Dataset,
      "; AUROC = ",
      sprintf("%.3f", mean_AUC),
      " ± ",
      sprintf("%.3f", sd_AUC)
    )
  ) %>%
  select(Dataset, top_n, Dataset_label)

roc_grid_summary_all <- roc_grid_summary_all %>%
  left_join(auc_label_df, by = c("Dataset", "top_n"))

roc_grid_raw_all <- roc_grid_raw_all %>%
  left_join(auc_label_df, by = c("Dataset", "top_n"))

plot_cols <- dataset_cols[auc_label_df$Dataset]
names(plot_cols) <- auc_label_df$Dataset_label

plot_cols <- plot_cols[!is.na(plot_cols)]

auc_compare_df <- roc_auc_raw_all %>%
  filter(
    Dataset %in% c(
      "pre-NEC-Species",
      "pre-NEC-AMPs"
    )
  ) %>%
  select(
    Dataset,
    Repeat_ID,
    AUC
  ) %>%
  distinct(
    Dataset,
    Repeat_ID,
    .keep_all = TRUE
  ) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = AUC
  ) %>%
  filter(
    !is.na(`pre-NEC-Species`),
    !is.na(`pre-NEC-AMPs`)
  ) %>%
  arrange(Repeat_ID)

cat("successful Repeat_ID numbers:", nrow(auc_compare_df), "\n")

if (nrow(auc_compare_df) < 2) {
  stop(
    "common Repeat_ID of the two models is less than 2"
  )
}

print(auc_compare_df)

auc_compare_summary <- tibble(
  Dataset = c(
    "pre-NEC-Species",
    "pre-NEC-AMPs"
  ),
  n = c(
    sum(!is.na(auc_compare_df$`pre-NEC-Species`)),
    sum(!is.na(auc_compare_df$`pre-NEC-AMPs`))
  ),
  mean_AUC = c(
    mean(auc_compare_df$`pre-NEC-Species`, na.rm = TRUE),
    mean(auc_compare_df$`pre-NEC-AMPs`, na.rm = TRUE)
  ),
  sd_AUC = c(
    sd(auc_compare_df$`pre-NEC-Species`, na.rm = TRUE),
    sd(auc_compare_df$`pre-NEC-AMPs`, na.rm = TRUE)
  ),
  median_AUC = c(
    median(auc_compare_df$`pre-NEC-Species`, na.rm = TRUE),
    median(auc_compare_df$`pre-NEC-AMPs`, na.rm = TRUE)
  ),
  IQR_AUC = c(
    IQR(auc_compare_df$`pre-NEC-Species`, na.rm = TRUE),
    IQR(auc_compare_df$`pre-NEC-AMPs`, na.rm = TRUE)
  )
)

print(auc_compare_summary)

auc_wilcox_result <- wilcox.test(
  x = auc_compare_df$`pre-NEC-Species`,
  y = auc_compare_df$`pre-NEC-AMPs`,
  paired = TRUE,
  exact = FALSE,
  alternative = "two.sided",
  conf.int = TRUE
)

auc_p_value <- auc_wilcox_result$p.value

auc_p_adjusted <- p.adjust(
  auc_p_value,
  method = "BH"
)

auc_compare_df <- auc_compare_df %>%
  mutate(
    AUC_difference =
      `pre-NEC-Species` - `pre-NEC-AMPs`
  )

mean_auc_difference <- mean(
  auc_compare_df$AUC_difference,
  na.rm = TRUE
)

median_auc_difference <- median(
  auc_compare_df$AUC_difference,
  na.rm = TRUE
)

auc_significance_result <- tibble(
  Comparison = "pre-NEC-Species vs pre-NEC-AMPs",
  Test = "Paired Wilcoxon signed-rank test",
  n_pairs = nrow(auc_compare_df),
  Species_mean_AUC = mean(
    auc_compare_df$`pre-NEC-Species`,
    na.rm = TRUE
  ),
  Species_sd_AUC = sd(
    auc_compare_df$`pre-NEC-Species`,
    na.rm = TRUE
  ),
  AMPs_mean_AUC = mean(
    auc_compare_df$`pre-NEC-AMPs`,
    na.rm = TRUE
  ),
  AMPs_sd_AUC = sd(
    auc_compare_df$`pre-NEC-AMPs`,
    na.rm = TRUE
  ),
  Mean_AUC_difference = mean_auc_difference,
  Median_AUC_difference = median_auc_difference,
  Wilcoxon_statistic = unname(auc_wilcox_result$statistic),
  P_value = auc_p_value,
  P_adjusted_BH = auc_p_adjusted,
  Significance = case_when(
    auc_p_adjusted < 0.0001 ~ "****",
    auc_p_adjusted < 0.001  ~ "***",
    auc_p_adjusted < 0.01   ~ "**",
    auc_p_adjusted < 0.05   ~ "*",
    TRUE                    ~ "ns"
  )
)

print(auc_significance_result)

cat("\n配对 Wilcoxon 检验结果:\n")
cat(
  "Species AUROC =",
  sprintf(
    "%.3f ± %.3f",
    auc_significance_result$Species_mean_AUC,
    auc_significance_result$Species_sd_AUC
  ),
  "\n"
)

cat(
  "AMPs AUROC =",
  sprintf(
    "%.3f ± %.3f",
    auc_significance_result$AMPs_mean_AUC,
    auc_significance_result$AMPs_sd_AUC
  ),
  "\n"
)

cat(
  "Average pairing difference Species - AMPs =",
  sprintf("%.4f", mean_auc_difference),
  "\n"
)

cat(
  "P value =",
  format.pval(
    auc_p_value,
    digits = 4,
    eps = 1e-4
  ),
  "\n"
)

#P value =0.05906
cat(
  "BH-adjusted P value =",
  format.pval(
    auc_p_adjusted,
    digits = 4,
    eps = 1e-4
  ),
  "\n"
)
#P value =0.05906

cat(
  "Significance =",
  auc_significance_result$Significance,
  "\n"
)

auc_p_label <- case_when(
  auc_p_adjusted < 0.0001 ~ "P < 0.0001",
  TRUE ~ paste0(
    "P = ",
    format.pval(
      auc_p_adjusted,
      digits = 3,
      eps = 1e-4
    )
  )
)

auc_test_label <- paste0(
  "Paired Wilcoxon test, ",
  auc_p_label
)

write_csv(
  auc_significance_result,
  file.path(
    out_dir,
    "pre_NEC_Species_vs_AMPs_AUC_significance.csv"
  )
)

p_roc_all <- ggplot(
  roc_grid_summary_all,
  aes(
    x = FPR,
    y = mean_TPR,
    color = Dataset_label,
    fill = Dataset_label,
    group = Dataset_label
  )
) +
  geom_ribbon(
    aes(
      ymin = lower_TPR,
      ymax = upper_TPR
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line(
    linewidth = 1.5
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 0.97,
    y = 0.05,
    label = auc_test_label,
    hjust = 1,
    vjust = 0,
    size = 3.8,
    color = "black"
  ) +
  scale_color_manual(values = plot_cols) +
  scale_fill_manual(values = plot_cols) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = "pre-NEC",
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.line = element_line(
      linewidth = 0.7
    ),
    axis.ticks = element_line(
      linewidth = 0.7
    ),
    legend.position = c(0.62, 0.20),
    legend.background = element_rect(
      fill = "white",
      color = "grey80",
      linewidth = 0.4
    ),
    legend.key = element_blank(),
    legend.text = element_text(size = 9.5)
  )






# ============================================================
#Top 20 AMP clustersfor and species Gini for Fig.6C-D
# ============================================================
library(tidyverse)

amp_file <-read_csv(
  "C:/Users/taoch/Desktop/632751_clusters/Fig_7/new_figure_for_NEC/preNEC-GINI/AMP_RF_best_topn_4096_feature_importance_Gini_top30.csv",
  show_col_types = FALSE
)

U18AMPdb_name <-read_csv(
  "C:/Users/taoch/Desktop/632751_clusters/Supplementary_data/Supplementary Data 6_632751_cluster_summary_revise.csv",
  show_col_types = FALSE
)
U18AMPdb_name <-U18AMPdb_name%>%select(Cluster_name,U18AMPdb_name)

amp_file_U18AMPdb <-amp_file %>%inner_join(U18AMPdb_name,by=c("amp_id"="Cluster_name"))

amp_file_U18AMPdb <-amp_file_U18AMPdb%>%
arrange(desc(mean_Gini)) %>%
  slice_head(n = 20) %>%
  mutate(
    U18AMPdb_name = factor(
      U18AMPdb_name,
      levels = U18AMPdb_name
    )
  )

amp_file_U18AMPdb <- amp_file_U18AMPdb %>%
  arrange(desc(mean_Gini)) %>%
  mutate(
    U18AMPdb_name = factor(
      U18AMPdb_name,
      levels = rev(U18AMPdb_name)
    )
  )

p_gini <- ggplot(
  amp_file_U18AMPdb,
  aes(
    x = U18AMPdb_name,
    y = mean_Gini
  )
) +
  geom_col(
    width = 0.72,
    fill = "#C5952B"
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "MeanDecreaseGini"
  ) +
  
  scale_y_continuous(
    limits = c(0, 2),
    breaks = seq(0, 2, by = 0.5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  theme_classic(base_size = 12) +
  theme
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title.y = element_blank(),
    
    axis.line = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5,
      color = "black"
    ),
    
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )



# ============================================================
#AMP clusters–bacteria association network for Fig.6E-F
# ============================================================
meta <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_7/gut_amp_with_taxa/postNEC-765/new-postNEC-765-metadata.csv", stringsAsFactors = FALSE)
amp_data <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_7/gut_amp_with_taxa/postNEC-765/postNEC_amp_filt_gut.csv", row.names = 1, check.names = FALSE)
bac_data <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_7/gut_amp_with_taxa/postNEC-765/postNEC_bac_species_filt_gut.csv", row.names = 1, check.names = FALSE)

library(WGCNA)
library(multtest)
library(dplyr)
library(igraph) 

colnames(meta) <- trimws(colnames(meta))
common_samples <- Reduce(intersect, list(colnames(amp_data), colnames(bac_data), meta$Run))

amp_mat_full <- as.matrix(amp_data[, common_samples, drop = FALSE])
bac_mat_full <- as.matrix(bac_data[, common_samples, drop = FALSE])
meta_full <- meta[match(common_samples, meta$Run), ]

cat(sprintf("Original Dimensions - AMP: %d, Bacteria: %d, Samples: %d/n", 
            nrow(amp_mat_full), nrow(bac_mat_full), ncol(amp_mat_full)))

keep_amp_base <- rowSums(amp_mat_full > 0) > (0.05 * ncol(amp_mat_full))
keep_bac_base <- rowSums(bac_mat_full > 0) > (0.05 * ncol(bac_mat_full))

amp_mat_filtered <- amp_mat_full[keep_amp_base, , drop = FALSE]
bac_mat_filtered <- bac_mat_full[keep_bac_base, , drop = FALSE]

cat(sprintf("After 5%% Prevalence Filter - AMP: %d, Bacteria: %d/n", 
            nrow(amp_mat_filtered), nrow(bac_mat_filtered)))

analyze_single_cohort <- function(cohort_name, meta_df, amp_mat, bac_mat) {
  
  samples_in_cohort <- meta_df %>% filter(Study == cohort_name) %>% pull(Run)
  samples_in_cohort <- intersect(samples_in_cohort, colnames(amp_mat))
  
  if (length(samples_in_cohort) < 10) {
    cat(sprintf(" Skip (%d)/n", cohort_name, length(samples_in_cohort)))
    return(NULL)
  }
  
  amp_sub <- amp_mat[, samples_in_cohort, drop = FALSE]
  bac_sub <- bac_mat[, samples_in_cohort, drop = FALSE]
  
  cor_res <- corAndPvalue(t(amp_sub), t(bac_sub), method = "spearman")
  cor_mat <- cor_res$cor
  p_mat <- cor_res$p
  
  df_flat <- data.frame(
    AMP = rep(rownames(cor_mat), times = ncol(cor_mat)),
    Bacteria = rep(colnames(cor_mat), each = nrow(cor_mat)),
    Correlation = as.vector(cor_mat),
    P.value = as.vector(p_mat)
  )
  
  df_flat <- na.omit(df_flat)
  if (nrow(df_flat) == 0) return(NULL)
  
  fdr_res <- mt.rawp2adjp(df_flat$P.value, proc = "BH")
  adj_p_sorted <- fdr_res$adjp[order(fdr_res$index), 2]
  df_flat$FDR <- adj_p_sorted
  
  sig_edges <- df_flat %>%
    filter(FDR < 0.05 & Correlation < -0.4) %>%
    mutate(Cohort = cohort_name) 
  
  cat(sprintf("  %s:%d significant/n", cohort_name, nrow(sig_edges)))
  
  return(sig_edges)
}

cohorts <- unique(meta_full$Study)
all_results_list <- list()

cat(sprintf("/n===Queue ===/n", length(cohorts)))

for (cq in cohorts) {
  res <- analyze_single_cohort(cq, meta_full, amp_mat_filtered, bac_mat_filtered)
  if (!is.null(res)) {
    all_results_list[[cq]] <- res
  }
}


if (length(all_results_list) == 0) {
  cat("no /n")
} else {
  all_edges <- bind_rows(all_results_list)

  min_cohorts_threshold <- 2
  
  final_edges <- all_edges %>%
    group_by(AMP, Bacteria) %>%
    summarise(
      Num_Cohorts = n(), # 
      Mean_Correlation = mean(Correlation),
      Min_FDR = min(FDR)
    ) %>%
    filter(Num_Cohorts >= min_cohorts_threshold) %>%
    arrange(desc(Num_Cohorts), desc(abs(Mean_Correlation)))
  


U18AMPdb <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Supplementary_data/Supplementary Data 4_632751_cluster_summary.csv")
U18AMPdb <-U18AMPdb %>%select(Cluster_name,U18AMPdb_name)

final_edges <-final_edges %>%
  left_join(U18AMPdb,by=c("AMP"="Cluster_name"))

write.csv(final_edges,"C:/Users/taoch/Desktop/632751_clusters/Fig_7/gut_amp_with_taxa/postNEC-765/Consensus_Network_Edges_gut_postNEC_study_with_2_new.csv",row.names = FALSE) 

final_edges <-read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_7/gut_amp_with_taxa/postNEC-765/Consensus_Network_Edges_gut_postNEC_study_with_2_new.csv")

if (nrow(final_edges) > 0) {
  g <- graph_from_data_frame(final_edges %>% select(AMP, Bacteria), directed = FALSE)
  
  E(g)$weight <- abs(final_edges$Mean_Correlation)
  
  V(g)$color <- ifelse(V(g)$name %in% final_edges$AMP,"#D29D21","#625C93")
  V(g)$size <- ifelse(V(g)$name %in% final_edges$AMP, 6, 10) 
  
  # 绘图
  pdf("C:/Users/taoch/Desktop/632751_clusters/Fig_7/gut_amp_with_taxa/postNEC-765/postNEC-765—Consensus_Network_study_with_2.pdf", 
      width = 12, height = 8)
  
  plot(g, 
       main = paste("Consensus AMP-Bacteria Network/n(Reproducible in >=", min_cohorts_threshold, "Cohorts)"),
       vertex.label = V(g)$name,        
       vertex.label.cex = 0.3,        
       vertex.label.color = "black",  
       vertex.label.font = 2,           
       
       vertex.color = V(g)$color,      
       vertex.size = V(g)$size,        
       edge.color = "grey50",
       edge.width = E(g)$weight * 3,
       layout = layout_with_fr)         
  
  dev.off()
}




# ============================================================
#ROC curves for Fig.6H-I；Fig.S14G-H
# ============================================================
library(tidyverse)
library(pROC)

pred_files <- tibble(
  Dataset = c(
    "Nasal S. aureus",
    "Skin AD",
    "Skin RAG"
  ),
  Best_topn = c(
    "2048",
    "2048",
    "2048",
    "2048"
  ),
  File = c(
    "C:/Users/taoch/Desktop/632751_clusters/Fig_7/nasal/new_topn_repeat10_nasal/amp/amp_RF_repeated5fold_5repeats_predictions_nasal.csv",
    "C:/Users/taoch/Desktop/632751_clusters/Fig_7/skin_AD/new_topn_repeat10_skin_AD/amp/amp_RF_repeated5fold_5repeats_predictions_skin_AD.csv",
    "C:/Users/taoch/Desktop/632751_clusters/Fig_7/skin_RAG/new_topn_repeat10_skin_RAG/amp/amp_RF_repeated5fold_5repeats_predictions_skin_RAG.csv"
  )
)

out_dir <- "C:/Users/taoch/Desktop/632751_clusters/Fig_7/amp_cluster_4datasets_mean_ROC"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


dataset_cols <- c(
  "Nasal S. aureus"      = "#816320",
  "Skin AD"             = "#8B575A",
  "Skin RAG"            = "#4F7F52"
)

calculate_mean_roc <- function(pred_raw_file, dataset_name, best_topn) {
  
  cat("\n============================================================\n")
  cat(" ", dataset_name, "\n", sep = "")
  cat("best_topn: ", best_topn, "\n", sep = "")
  cat("file: ", pred_raw_file, "\n", sep = "")
  cat("============================================================\n")
  
  pred_raw <- read_csv(
    pred_raw_file,
    show_col_types = FALSE
  )
  if (!"Repeat_ID" %in% colnames(pred_raw) && "repeat" %in% colnames(pred_raw)) {
    pred_raw <- pred_raw %>%
      rename(Repeat_ID = `repeat`)
  }
  
  required_cols <- c("Repeat_ID", "top_n", "Run", "y_true", "y_pred")
  missing_cols <- setdiff(required_cols, colnames(pred_raw))
  
  if (length(missing_cols) > 0) {
    stop(
      dataset_name,
      "colomn not exited: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  pred_raw <- pred_raw %>%
    mutate(
      Dataset = dataset_name,
      Repeat_ID = as.integer(Repeat_ID),
      top_n = as.character(top_n),
      y_true = as.character(y_true),
      y_pred = as.numeric(y_pred)
    )
  
  repeat_ids <- sort(unique(pred_raw$Repeat_ID))
  n_repeats <- length(repeat_ids)

  best_pred <- pred_raw %>%
    filter(top_n == best_topn)
  
  if (nrow(best_pred) == 0) {
    stop(dataset_name, ": best_topn = ", best_topn, " no exit")
  }
  
  fpr_grid <- seq(0, 1, length.out = 101)
  
  roc_grid_list <- list()
  roc_auc_list <- list()
  
  for (r in repeat_ids) {
    
    df_r <- best_pred %>%
      filter(Repeat_ID == r)
    
    if (nrow(df_r) == 0) {
      warning(dataset_name, ": Repeat_ID = ", r, " jump")
      next
    }
    
    if (length(unique(df_r$y_true)) < 2) {
      warning(dataset_name, ": Repeat_ID = ", r, "jump")
      next
    }
    
    y_r <- factor(
      df_r$y_true,
      levels = c("Health", "Disease")
    )
    
    pred_r <- df_r$y_pred
    
    roc_r <- roc(
      response = y_r,
      predictor = pred_r,
      levels = c("Health", "Disease"),
      direction = "<",
      quiet = TRUE
    )
    
    auc_r <- as.numeric(auc(roc_r))
    
    roc_auc_list[[as.character(r)]] <- tibble(
      Dataset = dataset_name,
      Repeat_ID = r,
      top_n = best_topn,
      AUC = auc_r
    )
    
    tpr_r <- coords(
      roc_r,
      x = 1 - fpr_grid,
      input = "specificity",
      ret = "sensitivity",
      transpose = FALSE
    )
    
    if (is.data.frame(tpr_r)) {
      tpr_r <- tpr_r$sensitivity
    } else {
      tpr_r <- as.numeric(tpr_r)
    }
    
    roc_grid_list[[as.character(r)]] <- tibble(
      Dataset = dataset_name,
      Repeat_ID = r,
      top_n = best_topn,
      FPR = fpr_grid,
      TPR = tpr_r
    )
  }
  
  roc_grid_raw <- bind_rows(roc_grid_list)
  roc_auc_raw <- bind_rows(roc_auc_list)
  
  if (nrow(roc_grid_raw) == 0) {
    stop(dataset_name, ": no ROC")
  }
  
  
  roc_grid_summary <- roc_grid_raw %>%
    group_by(Dataset, top_n, FPR) %>%
    summarise(
      mean_TPR = mean(TPR, na.rm = TRUE),
      sd_TPR = sd(TPR, na.rm = TRUE),
      lower_TPR = pmax(mean_TPR - sd_TPR, 0),
      upper_TPR = pmin(mean_TPR + sd_TPR, 1),
      .groups = "drop"
    )
  
  auc_summary <- roc_auc_raw %>%
    summarise(
      Dataset = dataset_name,
      top_n = best_topn,
      n_repeats = n(),
      mean_AUC = mean(AUC, na.rm = TRUE),
      sd_AUC = sd(AUC, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Mean AUROC =", round(auc_summary$mean_AUC, 4), "\n")
  cat("SD AUROC   =", round(auc_summary$sd_AUC, 4), "\n")
  
  return(
    list(
      roc_grid_raw = roc_grid_raw,
      roc_grid_summary = roc_grid_summary,
      roc_auc_raw = roc_auc_raw,
      auc_summary = auc_summary
    )
  )
}


all_results <- pmap(
  list(
    pred_raw_file = pred_files$File,
    dataset_name = pred_files$Dataset,
    best_topn = pred_files$Best_topn
  ),
  calculate_mean_roc
)

roc_grid_raw_all <- map_dfr(all_results, "roc_grid_raw")
roc_grid_summary_all <- map_dfr(all_results, "roc_grid_summary")
roc_auc_raw_all <- map_dfr(all_results, "roc_auc_raw")
auc_summary_all <- map_dfr(all_results, "auc_summary")


auc_label_df <- auc_summary_all %>%
  mutate(
    Dataset_label = paste0(
      Dataset,
      "; AUROC = ",
      sprintf("%.3f", mean_AUC),
      " ± ",
      sprintf("%.3f", sd_AUC)
    )
  ) %>%
  select(Dataset, top_n, Dataset_label)

roc_grid_summary_all <- roc_grid_summary_all %>%
  left_join(auc_label_df, by = c("Dataset", "top_n"))

roc_grid_raw_all <- roc_grid_raw_all %>%
  left_join(auc_label_df, by = c("Dataset", "top_n"))

plot_cols <- dataset_cols[auc_label_df$Dataset]
names(plot_cols) <- auc_label_df$Dataset_label


repeat_tag <- paste0(
  max(auc_summary_all$n_repeats, na.rm = TRUE),
  "repeats"
)

topn_tag <- "best_topn_2048"

roc_grid_raw_file <- file.path(
  out_dir,
  paste0("amp_cluster_RF_repeated5fold_", repeat_tag, "_", topn_tag, "_ROC_grid_raw_4datasets.csv")
)

roc_grid_summary_file <- file.path(
  out_dir,
  paste0("amp_cluster_RF_repeated5fold_", repeat_tag, "_", topn_tag, "_ROC_grid_summary_4datasets.csv")
)

roc_auc_raw_file <- file.path(
  out_dir,
  paste0("amp_cluster_RF_repeated5fold_", repeat_tag, "_", topn_tag, "_AUC_raw_4datasets.csv")
)

auc_summary_file <- file.path(
  out_dir,
  paste0("amp_cluster_RF_repeated5fold_", repeat_tag, "_", topn_tag, "_AUC_summary_4datasets.csv")
)

write_csv(roc_grid_raw_all, roc_grid_raw_file)
write_csv(roc_grid_summary_all, roc_grid_summary_file)
write_csv(roc_auc_raw_all, roc_auc_raw_file)
write_csv(auc_summary_all, auc_summary_file)


p_roc_all <- ggplot(
  roc_grid_summary_all,
  aes(
    x = FPR,
    y = mean_TPR,
    color = Dataset_label,
    fill = Dataset_label
  )
) +
  geom_ribbon(
    aes(
      ymin = lower_TPR,
      ymax = upper_TPR
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line(
    linewidth = 1.5
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.8
  ) +
  scale_color_manual(values = plot_cols) +
  scale_fill_manual(values = plot_cols) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = "Mean ROC curves for AMP cluster-based random forest",
    subtitle = "(top_n = 2048)",
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    legend.position = c(0.62, 0.20),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.4),
    legend.key = element_blank(),
    legend.text = element_text(size = 9.5)
  )

print(p_roc_all)


roc_pdf <- file.path(
  out_dir,
  paste0("amp_cluster_RF_repeated5fold_", repeat_tag, "_", topn_tag, "_mean_ROC_4datasets.pdf")
)

ggsave(
  filename = roc_pdf,
  plot = p_roc_all,
  width = 8,
  height = 6
)
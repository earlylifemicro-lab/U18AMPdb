# ============================================================
#Mother-to-infant transmission for Fig.5B-C
# ============================================================

library(tidyverse)
library(vegan)   
library(coin) 
library(rstatix)
#install.packages("coin")
library(ggplot2)
library(ggpubr)   

df_infant <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/2.shanno_index/infant_shanno_index/shannon_meta_infant_six_age_groups.csv")
df_mother <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/2.shanno_index/mother_shanno_index/shannon_meta_mother_three_periods.csv")

df_all <- df_all %>%
  mutate(
    Group = factor(
      Group,
      levels = c("mother", "infant")
    ),
    Study = factor(Study)
  )
write.csv(df_all,"C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/2.shanno_index/all_shanno_index.csv")


df_all$Group <- as.factor(df_all$Group)
df_all$Study <- as.factor(df_all$Study)
try_test <- wilcox_test(amp_shannon ~ Group | Study, data = df_all)

p_value <- pvalue(try_test)
print(p_value)

my_colors <- c("mother" = "#AEB2D1", "infant" = "#7C9895") 
df_all$Group <- factor(df_all$Group, levels = c("mother", "infant"))

p <- ggplot(df_all, aes(x = Group, y = amp_shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) + 
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_manual(values = my_colors, guide = "none") + 
  scale_color_manual(values = my_colors, guide = "none") + 
  scale_y_continuous(limits = c(8, 12), breaks = seq(8, 12, 1)) +
  labs(y = "Shannon diversity index", x = "") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
  )
p_final <- p + 
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",    
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
    vjust = 1,           
    size = 5              
  )

print(p_final)


df_mother$Group <- "mother"
df_mother <- df_mother %>%
  rename(TimePoint = Period) %>%
  mutate(TimePoint = factor(TimePoint, levels = c("Pregnancy", "Delivery", "Postpartum")))
df_infant$Group <- "infant"
df_infant <- df_infant %>%
  mutate(TimePoint = factor(TimePoint, levels = c("0", "1", "3", "6", "12", "18")))
plot_data <- bind_rows(df_mother, df_infant)
write.csv(plot_data, "C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/2.shanno_index/shanno_index_monther_and_infant.csv",row.names = FALSE)
plot_data$Subject <- factor(plot_data$Group, levels = c("mother", "infant"))
plot_data$Study <- factor(plot_data$Study)
sample_sizes <- plot_data %>%
  group_by(Subject, TimePoint) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(TimePoint, " (n=", n, ")"))
plot_data <- left_join(plot_data, sample_sizes, by = c("Subject", "TimePoint"))

# A. Kruskal-Wallis
p_overall<- kruskal_test(amp_shannon ~ TimePoint | Study, data = plot_data) 
p_overall_val <- pvalue(p_overall)
print(p_overall_val)

# B.Wilcoxon test
wilcox_results <- plot_data %>%
  group_by(Subject) %>%
  wilcox_test(
    amp_shannon ~ TimePoint, 
    p.adjust.method = "BH"
  add_significance()

print(wilcox_results)

wilcox_results_sig <- wilcox_results %>%
  filter(p.adj < 0.05)%>%
  add_xy_position(x = "TimePoint", dodge = 0.8)%>%  
  as.data.frame() %>%
  select(-y.position, -groups) 

write.csv(wilcox_results_sig, 
          "C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/2.shanno_index/infant_shanno_index/wilcox_results_sig.csv", 
          row.names = FALSE)

my_colors <- c("mother" = "#AEB2D1", "infant" = "#7C9895") 
plot_data$Group <- factor(plot_data$Group, levels = c("mother", "infant"))

p <- ggplot(plot_data, aes(x = TimePoint, y = amp_shannon, fill = Subject)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) + 
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.4, size = 1) +
  facet_grid(. ~ Subject, scales = "free_x", space = "free_x", switch = "x") +
  scale_y_continuous(limits = c(5, 15), breaks = seq(5, 15, 5)) +
  scale_fill_manual(values = my_colors, guide = "none") + 
  scale_color_manual(values = my_colors, guide = "none") +
  labs(y = "Shannon diversity index", x = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    strip.background = element_blank(), 
    strip.placement = "outside",
    legend.position = "none", 
    axis.line = element_line(color = "black"), 
    panel.grid = element_blank(), 
    panel.border = element_blank()
  )


# ============================================================
#CheckM for Fig.S12B
# ============================================================
library(ggplot2)
library(patchwork)

df <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/1.checkM/merged_plot.csv")

cols <- c("HQ" = "#B83945", "MQ" = "#7C9895", "LQ" = "grey")
p_scatter <- ggplot(df, aes(x = Completeness, y = Contamination, color = Quality)) +
  geom_point(alpha = 0.6, size = 0.8) +
  scale_color_manual(values = cols) +
  theme_minimal()+
  labs(
    x = "Completeness (%)",
    y = "Contamination (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
  )+
  xlim(20, 100) + ylim(0, 20)
print(p_scatter)

p_hist_x <- ggplot(df, aes(x = Completeness, fill = Quality)) +
  geom_histogram(binwidth = 5, position = "stack", color = "white", alpha = 0.8) +
  scale_fill_manual(values = cols) +
  theme_minimal()+
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.position = "none", 
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
    axis.ticks.y =element_line() 
  ) +
  xlim(20, 100) + ylab("Count")

p_hist_y <- ggplot(df, aes(x = Contamination, fill = Quality)) +
  geom_histogram(binwidth = 2, position = "stack", color = "white", alpha = 0.8) +
  scale_fill_manual(values = cols) +
  theme_minimal()+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 8),
    axis.ticks.x =element_line(), 
    legend.position = "right",
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
  ) +
  coord_flip() + 
  xlim(0, 20) + xlab("Count")

final_plot <- p_hist_x / (p_scatter | p_hist_y) + 
  plot_layout(heights = c(1, 6), widths = c(3, 1))

df_plot <- subset(df, Quality %in% c("MQ", "HQ"))

p_Contig_count <- ggplot(df_plot, aes(x = Quality, y = Contigs, fill = Quality)) +
  geom_violin(alpha = 0.7, scale = "width", trim = FALSE) + 
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.ticks =element_line(), 
    legend.position = "right",
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
  ) +
  scale_fill_manual(values = cols) +
  labs(y = "Number of contigs")+
  scale_y_continuous(labels = scales::comma) 
print(p_Contig_count)


p_N50_contig <- ggplot(df_plot, aes(x = Quality, y = N50_contig, fill = Quality)) +
  geom_violin(alpha = 0.7, scale = "width", trim = FALSE) + 
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.ticks =element_line(), 
    legend.position = "right", 
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
  ) +
  scale_fill_manual(values = cols) +
  labs(y = "Contigs N50")+
  scale_y_log10(labels = scales::comma) 
print(p_Contig_count)



# ============================================================
#PCoA of delivery-shared AMPs for Fig.5D
# ============================================================
library(vegan)
library(dplyr)
library(tidyverse)
library(coin)
library(readxl)
library(dplyr)
library(ggplot2)
library(openxlsx)


combined_df <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/mother_242_infant_553_transmission_AMPgene_matrix_U18AMPdb_sum_merged.csv",row.names = 1)
meta_mother <-read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/Dlivery_transmission_meta_mother_242.csv")
meta_infant <-read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/Dlivery_transmission_meta_infant_553.csv")
combined_matrix <- data.matrix(combined_df)
combined_matrix  <- t(combined_matrix)


sample_sums_combined <- rowSums(combined_matrix)
valid_samples_combined <- sample_sums_combined > 0


combined_matrix_clean <- combined_matrix[valid_samples_combined, ]
dist_bc_combined <- vegdist(combined_matrix_clean, method = "bray")

pcoa_result_combined <- cmdscale(dist_bc_combined, k=2, eig=TRUE)
pcoa_points_combined <- as.data.frame(pcoa_result_combined$points)
colnames(pcoa_points_combined) <- c("PCoA1", "PCoA2")
eig_total <- sum(pcoa_result_combined$eig)
pcoa1_var <- (pcoa_result_combined$eig[1] / eig_total) * 100 
pcoa2_var <- (pcoa_result_combined$eig[2] / eig_total) * 100 
sample_ids_combined <- rownames(pcoa_points_combined)



meta_all <- bind_rows(meta_mother, meta_infant)
pcoa_plot_data <- pcoa_points_combined %>% 
  rownames_to_column(var = "metagenomes") %>% 
  left_join(meta_all, by = "metagenomes")
print(paste("total:", nrow(pcoa_plot_data)))
pcoa_plot_data <- pcoa_plot_data %>%
  mutate(PlotGroup = ifelse(Category == "mother", "Mother", as.character(Age_month_group)))
pcoa_plot_data <- pcoa_plot_data %>%
  filter(!is.na(PlotGroup))




group_order <- c("Mother", "0", "1", "3", "6", "12", "18")
pcoa_plot_data$PlotGroup <- factor(pcoa_plot_data$PlotGroup, levels = group_order)
my_colors <- c(
  "Mother" = "#8A8FB8",   
  "0"          = "#84BEC8",
  "1"          = "#7AB6C0",
  "3"          = "#65A6A8",
  "6"          = "#5A9D9A",
  "12"         = "#478074",
  "18"         = "#3D7265"    
)
eig_total <- sum(pcoa_result_combined$eig)
pcoa1_var <- round((pcoa_result_combined$eig[1] / eig_total) * 100, 1)
pcoa2_var <- round((pcoa_result_combined$eig[2] / eig_total) * 100, 1)

p_scatter <- ggplot(pcoa_plot_data, aes(x = PCoA1, y = PCoA2, color = PlotGroup)) +
  geom_point(alpha = 0.8, size = 2) + 
  scale_color_manual(values = my_colors) + 
  labs(x = paste0("PCoA1 (", pcoa1_var, "%)"),
       y = paste0("PCoA2 (", pcoa2_var, "%)"),
       color = "PlotGroup") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA),
    panel.grid = element_blank()
  )

print(p_scatter)

group_order_2 <- c("mother", "infant")
pcoa_plot_data$Category <- factor(pcoa_plot_data$Category, levels = group_order_2)
my_colors_2 <- c("mother" = "#AEB2D1", "infant" = "#7C9895") 

p_box_y <- ggplot(pcoa_plot_data, aes(x = Category, y = PCoA2)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6, aes(fill = Category)) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.5, aes(color = Category)) +
  scale_fill_manual(values = my_colors_2) +
  scale_color_manual(values = my_colors_2) +
  labs(x = NULL, y = "PCoA2") +
  theme_minimal()+
  theme(
    scale_y_continuous(position = "right"), 
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
  )
print(p_box_y)

p_box_x <- ggplot(pcoa_plot_data, aes(x = Category, y = PCoA1)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.6, aes(fill = Category)) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.5, aes(color = Category)) +
  scale_fill_manual(values = my_colors_2) +
  scale_color_manual(values = my_colors_2) +
  labs(x = NULL, y = "PCoA1") +
  coord_flip() + 
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
  )
print(p_box_x)




# ============================================================
#Cumulative curve for Fig.5E
# ============================================================
library(readxl)
library(dplyr)
library(ggplot2)
library(openxlsx)

metadata <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Supplementary_data/Supplementary Data 11_gut_metadata_Zeng_2022_month.csv")
delivery_trans<- read_excel("C:/Users/taoch/Desktop/632751_clusters/Supplementary_data/Supplementary Data 17_transmission_AMPgene_revise_final.xlsx",sheet = "Delivery_transmission")
Pregnancy_trans<- read_excel("C:/Users/taoch/Desktop/632751_clusters/Supplementary_data/Supplementary Data 17_transmission_AMPgene_revise_final.xlsx",sheet = "Pregnancy_transmission")
Postpartum_trans<- read_excel("C:/Users/taoch/Desktop/632751_clusters/Supplementary_data/Supplementary Data 17_transmission_AMPgene_revise_final.xlsx",sheet = "Postpartum_transmission")

all_transmission <- bind_rows(delivery_trans, Pregnancy_trans, Postpartum_trans)

matching_indices <- match(all_transmission$sampleID_Infant, metadata$Sample)
all_transmission$infant_DOL <- metadata$DOL[matching_indices]


df_clean <- all_transmission %>%
  filter(!is.na(infant_DOL)) %>%
  mutate(infant_DOL = as.numeric(as.character(infant_DOL)))

df_daily_count <- df_clean %>%
  group_by(infant_DOL, Mother_Status) %>%
  summarise(daily_count = n(), .groups = 'drop')
write.csv(df_daily_count,"C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/transmission_events/transmission_daily_event/df_daily_count_genes.csv",row.names = FALSE)


df_cumulative <- df_clean %>%
  group_by(infant_DOL, Mother_Status) %>%
  summarise(daily_n = n(), .groups = 'drop') %>%
  group_by(Mother_Status) %>%
  arrange(infant_DOL) %>% 
  mutate(cumulative_count = cumsum(daily_n)) %>%
  ungroup()
write.csv(df_daily_count,"C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/transmission_events/transmission_daily_event/df_cumulative_gene.csv",row.names = FALSE)


bar_width <-5
right_axis_breaks <- seq(0, 2000, by = 1000) 

p <- ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 0.4, linetype = "dashed") +
  
  geom_col(data = df_daily_count, 
           aes(x = infant_DOL, y = -daily_count, fill = Mother_Status), 
           width = bar_width, 
           alpha = 0.6,
           position = "identity") + 
  
  geom_line(data = df_cumulative, 
            aes(x = infant_DOL, y = cumulative_count, color = Mother_Status, group = Mother_Status), 
            alpha = 0.9, size = 1) +
  
  geom_point(data = df_cumulative, 
             aes(x = infant_DOL, y = cumulative_count, color = Mother_Status), 
             alpha = 0.5, size = 1) +
  
  labs(
    x = "Days of life (DOL)",
    y = "Cumulative acquisition of shared AMP genes", 
    color = "Mother Sampled Period",
    fill = "Daily transmission" 
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y.left = element_text(size = 10, color = "black"),
    axis.text.y.right = element_text(size = 10, color = "black"),
    axis.title.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    axis.ticks = element_line(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA),
    panel.grid = element_blank()
  ) +
  
  scale_color_manual(
    values = c("Pregnancy" = "#8B575A","Delivery" = "#8A8FB8", "Postpartum" = "#C3A46C")
  ) +
  scale_fill_manual(
    values = c("Pregnancy" = "#8B575A", "Delivery" = "#8A8FB8", "Postpartum" = "#C3A46C")
  ) +
  
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400)) +
  scale_y_continuous(
    name = "Cumulative acquisition of shared events",
    sec.axis = sec_axis(
      trans = ~ . * -1, 
      name = "Daily transmission events", 
      breaks = right_axis_breaks 
    )
  )


# ============================================================
#Longitudinal AMP Shannon diversity and AMP abundance for Fig.S12C-E
# ============================================================
library(tidyverse)

shannon_data_mother <- read.csv(
  paste0("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/2.shanno_index/mother_shanno_index/shannon_meta_mother_three_periods.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

subject_col <- "SubjectID"

period_levels <- c(
  "Pregnancy",
  "Delivery",
  "Postpartum"
)

paired_data <- shannon_data_mother %>%
  mutate(
    Period = trimws(
      as.character(Period)
    ),
    
    Period = factor(
      Period,
      levels = period_levels,
      ordered = TRUE
    ),
    
    Study = trimws(
      as.character(Study)
    ),
    
    mother_ID_original = trimws(
      as.character(
        .data[[subject_col]]
      )
    ),
    
    mother_ID = paste(
      Study,
      mother_ID_original,
      sep = "__"
    ),
    
    amp_shannon = as.numeric(
      amp_shannon
    )
  ) %>%
  
  filter(
    !is.na(Study),
    Study != "",
    !is.na(mother_ID_original),
    mother_ID_original != "",
    !is.na(Period),
    is.finite(amp_shannon)
  ) %>%
  
  group_by(
    Study,
    mother_ID_original,
    mother_ID,
    Period
  ) %>%
  
  summarise(
    amp_shannon = mean(
      amp_shannon,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

valid_mothers <- paired_data %>%
  group_by(
    mother_ID
  ) %>%
  summarise(
    number_of_periods = n_distinct(
      Period
    ),
    .groups = "drop"
  ) %>%
  filter(
    number_of_periods >= 2
  ) %>%
  pull(
    mother_ID
  )

paired_plot_data <- paired_data %>%
  filter(
    mother_ID %in% valid_mothers
  ) %>%
  arrange(
    Study,
    mother_ID,
    Period
  )
write.csv(paired_plot_data,"C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/2.shanno_index/paired_plot_data.csv",row.names = FALSE)

sample_number <- paired_plot_data %>%
  count(
    Period,
    name = "n"
  ) %>%
  complete(
    Period = factor(
      period_levels,
      levels = period_levels,
      ordered = TRUE
    ),
    fill = list(n = 0)
  )

print(sample_number)

period_colors <- c(
  "Pregnancy"  = "#8B575A",
  "Delivery"   = "#8A8FB8",
  "Postpartum" = "#C3A46C")


p1 <- ggplot(
  paired_plot_data,
  aes(
    x = Period,
    y = amp_shannon,
    group = mother_ID
  )
) +
  
  geom_line(
    color = "grey60",
    alpha = 0.45,
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  
  geom_point(
    aes(color = Period),
    shape = 16,
    size = 2.5,
    alpha = 0.85,
    na.rm = TRUE
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    color = "black",
    linewidth = 1.2,
    na.rm = TRUE
  ) +
  
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.5,
    fill = "white",
    color = "black",
    stroke = 0.8,
    na.rm = TRUE
  ) +
  
  scale_x_discrete(
    limits = c(
      "Pregnancy",
      "Delivery",
      "Postpartum"
    ),
    drop = FALSE
  ) +
  scale_color_manual(
    values = period_colors
  ) +
  
  labs(
    x = NULL,
    y = "Shannon index",
    title = "Within-mother change in AMP Shannon diversity"
  ) +
  
  guides(
    color = "none"
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    # 添加完整黑色边框
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.6,
      fill = NA
    ),
    
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 12
    ),
    
    axis.text.y = element_text(
      color = "black",
      size = 11
    ),
    
    axis.title.y = element_text(
      color = "black",
      size = 13
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    )
  )



# ============================================================
#Top30 shared AMP clusters for Fig.5F;Fig.S13AB
# ============================================================
library(tidyverse)


df_delivery <- read.csv(
  "C:/Users/taoch/Desktop/632751_clusters/Fig_6/linux_data/transmission/Delivery/Delivery_transmission_cluster_info.csv",
  check.names = FALSE
)

df_delivery_top30 <- df_delivery %>%
  count(U18AMPdb_name, name = "Total_shared_events") %>%
  arrange(desc(Total_shared_events)) %>%  
  head(30)                     
write.csv(df_delivery_top30,"C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/3.transmission/transmission_events/df_delivery_top30.csv",row.names = FALSE)


# 绘
df_delivery_top30 <-df_delivery_top30 %>%
  mutate(

    U18AMPdb_name = factor(U18AMPdb_name, levels = rev(U18AMPdb_name))
  ) 
p <- ggplot(
  df_delivery_top30,
  aes(
    x = Total_shared_events,
    y = U18AMPdb_name
  )
) +
  geom_col(
    width = 0.75,
    fill = "#8A8FB8",
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = Total_shared_events),
    hjust = -0.2,
    size = 3
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 6),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = "Total_shared_events",
    y = "shared_AMPs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.y = element_text(
      margin = margin(r = 10)
    ),
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    plot.margin = margin(
      t = 10,
      r = 25,
      b = 10,
      l = 10
    )
  ) 



# ============================================================
# Longitudinal prevalence and abundance for Fig.5G
# ============================================================
library(tidyverse)

df_infant <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/Dlivery_transmission_meta_infant_553.csv", header = TRUE, stringsAsFactors = FALSE)
df_matrix <- read_csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/mother_242_infant_553_transmission_AMPgene_matrix_U18AMPdb_sum_merged.csv", 
                      col_names = TRUE)%>% column_to_rownames(var = "contig") 
df_infant$metagenomes <- as.character(df_infant$metagenomes)
age_groups <- c("0", "1", "3", "6", "12", "18")
df_infant$Age_month_group <- factor(df_infant$Age_month_group, levels = age_groups)

sample_ids_matrix <- colnames(df_matrix)

matched_info <- df_infant %>%
  filter(metagenomes %in% sample_ids_matrix) %>%
  select(metagenomes, Age_month_group) %>%
  filter(Age_month_group!= "")

sample_counts <- matched_info %>%
  count(Age_month_group, name = "Sample_Count")


df_matrix_long <- df_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "metagenomes", values_to = "Abundance")
head(df_matrix_long)

stats_df_fast <- df_matrix_long %>%
  left_join(matched_info, by = "metagenomes") %>%
  filter(!is.na(Gene)) %>% 
  filter(!is.na(Age_month_group)) %>% 
  group_by(Gene, Age_month_group) %>%
  summarise(
    Prevalence = sum(Abundance > 0) / n(),
    Mean_Abundance = mean(Abundance),     
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Age_month_group,
    values_from = c(Prevalence, Mean_Abundance),
    names_sep = "_month_"
  )
write.csv(stats_df_fast, "C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/F_rpkm_preval_with_age_month/stats_df_fast_infant.csv")

stats_df_fast_mother <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/F_rpkm_preval_with_age_month/stats_df_fast_mother.csv")
stats_df_fast <-stats_df_fast%>% 
  left_join(stats_df_fast_mother,by="Gene")

top_AMP_id_30 <- stats_df_fast %>%
  rowwise() %>%
  mutate(
    avg_prev = mean(c_across(starts_with("Prevalence"))),
    avg_abun = mean(c_across(starts_with("Mean_Abundance")))
  ) %>%
  ungroup() %>%
  arrange(desc(avg_prev), desc(avg_abun)) %>%
  slice_head(n = 30) %>% # 取前30
  pull(Gene)

plot_data_Top_30 <- stats_df_fast%>%
  filter(Gene %in% top_AMP_id_30) 

write.csv(plot_data_Top_30,"C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/F_rpkm_preval_with_age_month/plot_data_Top_30.csv")


plot_data_long <- plot_data_Top_30 %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Metric_Month",
    values_to = "Value"
  ) %>%
  separate(Metric_Month, into = c("Metric", "Month"), sep = "_month_") %>%
  mutate(Month = factor(Month, levels = c("mother","0", "1", "3", "6", "12", "18")))

plot_data_long <- plot_data_long %>%
  group_by(Gene) %>%
  mutate(mean_abun_overall = mean(Value[Metric == "Mean_Abundance"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Gene = reorder(Gene, mean_abun_overall))
plot_data_long$Metric <- factor(plot_data_long$Metric, levels = c("Prevalence", "Mean_Abundance"))

max_abun <- max(plot_data_long$Value[plot_data_long$Metric == "Mean_Abundance"], na.rm = TRUE)
max_prev <- max(plot_data_long$Value[plot_data_long$Metric == "Prevalence"], na.rm = TRUE)

if(max_prev == 0) max_prev <- 1
scale_factor <- max_abun / max_prev

plot_data_long <- plot_data_long %>%
  mutate(
    Value_Scaled = ifelse(Metric == "Prevalence", Value * scale_factor, Value)
  )

p <- ggplot(plot_data_long, aes(x = Gene, y = Value_Scaled, fill = Metric)) +
  geom_col(
    data = subset(plot_data_long, Metric == "Mean_Abundance"),
    aes(y = Value), 
    fill = "#5F7A78", 
    alpha = 0.6,
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = subset(plot_data_long, Metric == "Prevalence"),
    aes(y = Value_Scaled),
    color = "#E69F00", 
    size = 2,
    show.legend = FALSE
  ) +
  geom_line(
    data = subset(plot_data_long, Metric == "Prevalence"),
    aes(y = Value_Scaled, group = 1), 
    color = "#E69F00",
    alpha = 0.5
  ) +
  facet_wrap(~ Month, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  scale_y_continuous(
    name = "Mean_Abundance",
    sec.axis = sec_axis(~ ./scale_factor, name = "Prevalence (line)")
  ) +
  labs(
    title = "Top 30 Genes: Prevalence and Abundance by Month",
    x = "AMP"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "#D3D3D3", color = "black"), 
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.y.right = element_text(color = "#56B4E9"),
    axis.title.y.right = element_text(color = "#56B4E9"),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    axis.ticks = element_line(color = "black", size = 0.4)
  )




# ============================================================
# Species-AMP clusters network for Fig.5H
# ============================================================
library(tidyverse)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)

species <- read.csv(
"C:/Users/taoch/Desktop/632751_clusters/Fig_6/linux_data/transmission/filter_transmisson_gene_rpkm_taxa_agegroup.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

top_30_delivery_cluster <- read.csv(
 "C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/transmission_events/transmission_top30/df_delivery_trans_top30.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

out_dir <- paste0(
  "C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/Species/Mother_VS_infant_carry_species/"
)

target_amp <- c(
  "U18AMPdb_008801",
  "U18AMPdb_015852"
)

layout_seed <- 20260712
set.seed(layout_seed)

required_species_columns <- c(
  "U18AMPdb_name",
  "Species",
  "RPKM"
)

missing_amp_top30 <- setdiff(
  target_amp,
  unique(top_30_delivery_cluster$U18AMPdb_name)
)

plot_data_selected <- species %>%
  filter(
    U18AMPdb_name %in% target_amp,
    !is.na(Species),
    Species != ""
  ) %>%
  mutate(
    infant_RPKM = suppressWarnings(as.numeric(RPKM))
  )

edge_summary <- plot_data_selected %>%
  group_by(U18AMPdb_name,Species) %>%
  summarise(Total_infant_RPKM = sum(infant_RPKM,na.rm = TRUE),
    
    Mean_infant_RPKM = if (
      all(is.na(infant_RPKM))
    ) {
      NA_real_
    } else {
      mean(
        infant_RPKM,
        na.rm = TRUE
      )
    },
    
    Number_of_records = n(),
    
    .groups = "drop"
  ) %>%
  arrange(
    factor(U18AMPdb_name, levels = target_amp),
    desc(Total_infant_RPKM),
    Species
  )


carrier_species_summary_selected <- edge_summary %>%
  group_by(U18AMPdb_name) %>%
  summarise(
    Number_of_carrier_species = n_distinct(Species),
    
    Carrier_species = paste(
      sort(unique(Species)),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    factor(U18AMPdb_name, levels = target_amp)
  )

print(carrier_species_summary_selected)

species_membership <- edge_summary %>%
  group_by(Species) %>%
  summarise(
    Number_of_AMP_clusters = n_distinct(U18AMPdb_name),
    
    Has_008801 = any(
      U18AMPdb_name == target_amp[1]
    ),
    
    Has_015852 = any(
      U18AMPdb_name == target_amp[2]
    ),
    
    AMP_membership = paste(
      sort(unique(U18AMPdb_name)),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Species_group = case_when(
      Has_008801 & Has_015852 ~
        "Shared by both",
      
      Has_008801 & !Has_015852 ~
        paste0(target_amp[1], " only"),
      
      !Has_008801 & Has_015852 ~
        paste0(target_amp[2], " only"),
      
      TRUE ~ "Unclassified"
    )
  ) %>%
  arrange(
    Species_group,
    Species
  )

shared_species <- species_membership %>%
  filter(Species_group == "Shared by both")

amp_nodes <- tibble(
  name = target_amp,
  Node_type = "AMP cluster",
  Node_group = target_amp,
  Number_of_AMP_clusters = NA_integer_,
  AMP_membership = target_amp,
  Node_size = 10,
  fontface = "bold"
)

species_nodes <- species_membership %>%
  transmute(
    name = Species,
    Node_type = "Species",
    Node_group = Species_group,
    Number_of_AMP_clusters,
    AMP_membership,
    Node_size = if_else(
      Species_group == "Shared by both",
      5.5,
      4.5
    ),
    fontface = "italic"
  )

nodes <- bind_rows(
  amp_nodes,
  species_nodes
)

node_group_levels <- c(
  target_amp[1],
  target_amp[2],
  "Shared by both",
  paste0(target_amp[1], " only"),
  paste0(target_amp[2], " only")
)

nodes <- nodes %>%
  mutate(
    Node_group = factor(
      Node_group,
      levels = node_group_levels
    ),
    
    Node_type = factor(
      Node_type,
      levels = c(
        "AMP cluster",
        "Species"
      )
    )
  )


make_vertical_positions <- function(
    dat,
    x_value,
    y_step = 1
) {
  
  dat <- dat %>%
    arrange(name)
  
  number_nodes <- nrow(dat)
  
  if (number_nodes == 0) {
    return(
      dat %>%
        mutate(
          x = numeric(0),
          y = numeric(0)
        )
    )
  }
  
  if (number_nodes == 1) {
    y_values <- 0
  } else {
    y_values <- seq(
      from = (number_nodes - 1) * y_step / 2,
      to = -(number_nodes - 1) * y_step / 2,
      length.out = number_nodes
    )
  }
  
  dat %>%
    mutate(
      x = x_value,
      y = y_values
    )
}

amp_positions <- nodes %>%
  filter(Node_type == "AMP cluster") %>%
  arrange(
    factor(name, levels = target_amp)
  ) %>%
  mutate(
    x = c(-2.2, 2.2),
    y = c(0, 0)
  )
shared_positions <- nodes %>%
  filter(Node_group == "Shared by both") %>%
  make_vertical_positions(
    x_value = 0,
    y_step = 1.25
  )
only_008801_positions <- nodes %>%
  filter(
    Node_group == paste0(
      target_amp[1],
      " only"
    )
  ) %>%
  make_vertical_positions(
    x_value = -5.5,
    y_step = 1
  )
only_015852_positions <- nodes %>%
  filter(
    Node_group == paste0(
      target_amp[2],
      " only"
    )
  ) %>%
  make_vertical_positions(
    x_value = 5.5,
    y_step = 1
  )

node_positions <- bind_rows(
  amp_positions,
  shared_positions,
  only_008801_positions,
  only_015852_positions
)
missing_node_positions <- setdiff(
  nodes$name,
  node_positions$name
)

if (length(missing_node_positions) > 0) {
  stop(
    paste0(
      "no axis：",
      paste(missing_node_positions, collapse = ", ")
    )
  )
}

network_edges <- edge_summary %>%
  transmute(
    from = U18AMPdb_name,
    to = Species,
    
    AMP_group = factor(
      U18AMPdb_name,
      levels = target_amp
    ),
    
    Total_infant_RPKM,
    Mean_infant_RPKM,
    Number_of_records
  )

if (
  all(
    network_edges$Total_infant_RPKM <= 0 |
    is.na(network_edges$Total_infant_RPKM)
  )
) {
  
  warning(
    "Total_infant_RPKM is 0"
  )
  
  network_edges <- network_edges %>%
    mutate(
      Edge_width_value = 1
    )
  
  edge_legend_title <- "AMP–Species association"
  
} else {
  
  network_edges <- network_edges %>%
    mutate(
      Edge_width_value = Total_infant_RPKM
    )
  
  edge_legend_title <- "Summed infant RPKM"
}

network_graph <- graph_from_data_frame(
  d = network_edges,
  directed = FALSE,
  vertices = node_positions
)

V(network_graph)$x <- node_positions$x[
  match(
    V(network_graph)$name,
    node_positions$name
  )
]

V(network_graph)$y <- node_positions$y[
  match(
    V(network_graph)$name,
    node_positions$name
  )
]

network_layout <- create_layout(
  network_graph,
  layout = "manual",
  x = V(network_graph)$x,
  y = V(network_graph)$y
)

amp_edge_colors <- setNames(
  c(
    "#B94C4C",
    "#3E79A8"
  ),
  target_amp
)

node_colors <- setNames(
  c(
    "#B94C4C",  
    "#3E79A8",  
    "#8A70B2",  
    "#D8C18C",  
    "#8DBEAD"   
  ),
  node_group_levels
)


p_network <- ggraph(network_layout) +
  
geom_edge_link(
  aes(
    edge_width = Edge_width_value,
    edge_colour = AMP_group
  ),
  alpha = 0.62,
  lineend = "round",
  show.legend = TRUE
) +
  
geom_node_point(
  aes(
    size = Node_size,
    shape = Node_type,
    fill = Node_group
  ),
  colour = "black",
  stroke = 0.45
) +
  
geom_node_text(
  aes(
    label = name,
    fontface = fontface
  ),
  repel = TRUE,
  seed = layout_seed,
  size = 3.3,
  colour = "black",
  box.padding = 0.30,
  point.padding = 0.20,
  min.segment.length = 0,
  max.overlaps = Inf,
  segment.colour = "grey65",
  segment.linewidth = 0.3
) +

scale_edge_colour_manual(
  values = amp_edge_colors,
  name = "AMP cluster"
) +
  
scale_edge_width_continuous(
  range = c(0.45, 4),
  trans = "sqrt",
  name = edge_legend_title
) +
scale_fill_manual(
  values = node_colors,
  drop = FALSE,
  name = "Node category"
) +
  scale_shape_manual(
    values = c(
      "AMP cluster" = 21,
      "Species" = 22
    ),
    guide = "none"
  ) +
  
  scale_size_identity() +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(0.16, 0.16)
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.08, 0.08)
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    title = paste0(
      target_amp[1],
      " and ",
      target_amp[2],
      " carrier-species network"
    ),
    
    subtitle = paste0(
      "Shared species are positioned between the two AMP clusters; ",
      "edge width represents summed infant RPKM"
    )
  ) +
  
  guides(
    edge_colour = guide_legend(
      order = 1,
      override.aes = list(
        edge_width = 2,
        alpha = 1
      )
    ),
    
    edge_width = guide_legend(
      order = 2
    ),
    
    fill = guide_legend(
      order = 3,
      override.aes = list(
        shape = 21,
        size = 5
      )
    )
  ) +
  
  theme_void(base_size = 12) +
  
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(
        b = 5
      )
    ),
    
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5,
      colour = "grey30",
      margin = margin(
        b = 10
      )
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    plot.margin = margin(
      t = 20,
      r = 90,
      b = 20,
      l = 90
    )
  )

ggsave(
  filename = file.path(
    out_dir,
    "AMP_species_network_008801_015852.pdf"
  ),
  plot = p_network,
  width = 13,
  height = 10,
  units = "in",
  bg = "white"
)





# ============================================================
# Time-stratified effects for Fig.5I
# ============================================================
library(vegan)
library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)

dist_df <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/G_tree/sample_amp_braycurtis_dist_delivery_share_gene_new.csv", row.names = 1)
dist_matrix <- as.dist(as.matrix(dist_df)) 
dist_labels <- attr(dist_matrix, "Labels") 
cat("Distance matrix samples:", length(dist_labels), "/n") 

meta <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/Dlivery_transmission_meta_infant_553.csv")
rownames(meta) <- meta$metagenomes
#551
meta_clean <- meta %>%
  mutate(Age_month = suppressWarnings(as.numeric(as.character(Age_month_group)))) %>%
  filter(metagenomes %in% dist_labels & !is.na(Age_month_group))
cat("Meta after cleaning:", nrow(meta_clean), "/n") #549

common_samples <- intersect(dist_labels, meta_clean$metagenomes)
dist_sub <- as.dist(
  as.matrix(dist_matrix)[common_samples, common_samples]
)

meta_sub <- meta_clean %>%
  filter(metagenomes %in% common_samples) %>%
  arrange(match(metagenomes, common_samples))

stopifnot(
  identical(rownames(as.matrix(dist_sub)), meta_sub$metagenomes)
)


covariates <- c("Delivery", "Country", "Gender", "Term", "Feed", "Study")

meta_sub$Age_month_group <- as.character(meta_sub$Age_month_group)

results <- tibble()


meta_infant <- meta_sub %>% filter(Category == "infant")

time_points <- sort(unique(meta_infant$Age_month_group))

for (time in time_points) {
  if (is.na(time)) next
  
  sub_meta <- meta_infant %>% filter(Age_month_group == time)
  n_samples <- nrow(sub_meta)
  
  cat("Processing Infant Timepoint:", time, "| n =", n_samples, "/n")
  if (n_samples < 10) next 

  sub_dist_mat <- as.matrix(dist_sub)[rownames(sub_meta), rownames(sub_meta)]
  sub_dist <- as.dist(sub_dist_mat)
  
  if (any(is.na(sub_dist))) next
  
  for (var in covariates) {
    if (!var %in% colnames(sub_meta)) next
    x <- sub_meta[[var]]
    if (all(is.na(x)) || length(unique(na.omit(x))) < 2) next
    
    sub_meta[[var]] <- as.factor(x)
    formula_str <- as.formula(paste("sub_dist ~", var))
    
    tryCatch({
      model <- adonis2(formula_str, data = sub_meta, permutations = 1000, by = "margin")
      r2 <- model$R2[1]
      pval <- model$`Pr(>F)`[1]

      time_lbl <- paste0(time, " mon.")
      
      results <- bind_rows(results, tibble(
        Group = "Infant",
        Time_Label = time_lbl,
        variable = var,
        R2 = r2,
        p = pval,
        n = n_samples
      ))
    }, error = function(e) {})
  }
}


results_final <- results %>%
  group_by(Group, Time_Label) %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    sig = case_when(
      FDR < 0.001 ~ "****",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    ),
    Time_Label = factor(Time_Label, levels = c("0 mon.", "1 mon.", "3 mon.", "6 mon.", "12 mon.", "18 mon.")),
    Group = factor(Group, levels = c("Mother", "Infant")),
    variable = factor(variable, levels = rev(c("Term", "Study", "Gender", "Feed", "Delivery", "Country")))
  )


write.csv(results_final,"C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/H_adonis2_with_age/adonis_Age_cocovarition.csv")

p <- ggplot(results_final, aes(x = R2, y = variable)) +
  geom_bar(stat = "identity", fill = "#D3D3D3", color = "black", width = 0.7) +
  geom_text(aes(label = sig), hjust = -0.1, size = 3.5, fontface = "bold") +
  facet_grid(Group ~ Time_Label, scales = "free_x", space = "free_x", switch = "y") +
  
  scale_x_continuous(
    limits = c(0, 0.20),
    breaks = seq(0, 0.2, 0.05),
    expand = expansion(mult = c(0, 0.1)) 
  ) +
  
  labs(x = "R2 (Variance Explained)", y = "") +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    strip.background = element_rect(fill = "#E0E0E0", color = "black"),
    strip.text = element_text(face = "bold", size = 10, color = "black"),

    strip.background.y = element_blank(), 
    strip.text.y = element_text(angle = 180, face = "bold", size = 14),

    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_blank(), 
    axis.line.y = element_blank(),  
    axis.ticks = element_line(),
  
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(1, 1, 1, 1, "cm") 
  )



# ============================================================
# Transmission rate for Fig.5J
# ============================================================
library(tidyr)
df <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/I-K--vaginal_vs_CS/I_transmission_protein/transmission_AMPgene.csv")
df_count <- df %>%
  group_by(sampleID_Infant, Study) %>%
  count(name = "count_100")

df_paire <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/I-K--vaginal_vs_CS/I_transmission_protein/infant_ids_paire.csv")
df_count_paire <- df_paire  %>%
  group_by(sampleID_infant, Study) %>%
  count(name = "Total_Count")  

df_rate <- df_count %>%
  inner_join(df_count_paire,by=c("sampleID_Infant"="sampleID_infant"))%>%
  select(-Study.y)%>%
  mutate(transmission_rate = (count_100/ Total_Count) *100)%>%
  mutate(transmission_rate = round(transmission_rate, 2))

infant_metadata <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/Dlivery_transmission_meta_infant_553.csv")
infant_metadata <-infant_metadata%>%
  select("Delivery","metagenomes")

df_rate_delivery <-df_rate %>%
  inner_join(infant_metadata,by=c("sampleID_Infant"="metagenomes"))

wilcox_result <- wilcox.test(transmission_rate ~ Delivery, data = df_rate_delivery, exact = FALSE)

p_value_raw <- wilcox_result$p.value

p_value_adjusted <- p.adjust(p_value_raw, method = "BH")


library(ggplot2)
library(ggpubr)
p_label <- ifelse(
  p_value_adjusted < 0.001,
  "BH-adjusted P < 0.001",
  paste0("BH-adjusted P = ", signif(p_value_adjusted, 3))
)

y_max <- max(df_rate_delivery$transmission_rate, na.rm = TRUE)
y_min <- min(df_rate_delivery$transmission_rate, na.rm = TRUE)
y_range <- y_max - y_min

y_line <- y_max + 0.08 * y_range
y_text <- y_max + 0.13 * y_range

colors_delivery <- c(
  "vaginal" = "#C25759",
  "CS" = "#377483"
)

p_simple <- ggplot(
  df_rate_delivery,
  aes(x = Delivery, y = transmission_rate, fill = Delivery, color = Delivery)
) +
  geom_violin(
    alpha = 0.5,
    trim = FALSE,
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.7
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = y_line,
    yend = y_line,
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = y_line - 0.02 * y_range,
    yend = y_line,
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = y_line - 0.02 * y_range,
    yend = y_line,
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_text,
    label = p_label,
    size = 4
  ) +
  
  scale_fill_manual(values = colors_delivery) +
  scale_color_manual(values = colors_delivery) +
  
  labs(
    x = "Delivery",
    y = "Transmission rate"
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.border = element_rect(
      colour = "black",
      linewidth = 0.6,
      fill = NA
    ),
    panel.grid = element_blank()
  )

print(p_simple)





# ============================================================
# Age-associated changes in shared AMP for Fig.5K;FigS13.E
# ============================================================
library(ggplot2)
library(tidyverse)
data <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/3.transmission/Dlivery_transmission/I-K--vaginal_vs_CS/J--gene_count_per_sample_with_age/shared_ampgene_total_rpkm_and_richness.csv")
data$log10_RPKM <- as.numeric(data$log10_RPKM)
data <- data %>% drop_na(Age_month_group)

p1 <- ggplot(data, aes(x = Age_month_group, y = log10_RPKM)) +
  geom_point(aes(color = log10_RPKM), size = 2, alpha = 0.6) +
  scale_color_gradient(low = "#299d8f", 
                       high = "#f3a361", 
                       name = "Log10(total rpkm+1)",
                       limits = c(0, 6),
                       breaks = seq(0, 6, 2) 
  )+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = TRUE, level = 0.99, linewidth = 0.9) +  
  labs(
    x = "Age (months)",
    y = "Log10(total rpkm + 1)"
  ) +
  scale_x_continuous(breaks = c(0,1, 3, 6, 12, 18))+
  scale_y_continuous(breaks = seq(0,6,2),
                     limits = c(0,6))+
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.3),
    axis.title = element_text(face = "bold", size=12),
    axis.text = element_text(color = "black", size = 11),
    axis.ticks = element_line(),
    axis.line = element_line(color="black",size=0.5),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )
print(p1)

fit <- lm(log10_RPKM ~ Age_month_group, data = data)
s <- summary(fit)

label_text <- paste0(
  "y = ",
  round(coef(fit)[2], 8), "x + ", round(coef(fit)[1], 3), "\n",
  "p = ", signif(s$coefficients[2, 4], 3)
)

p2 <- p1 + annotate(
  "text",
  x = Inf,
  y = Inf,
  label = label_text,
  hjust = 1.1,
  vjust = 1.1,
  size = 2.5
)




# ============================================================
#Species-level shared AMP for Fig.5L
# ============================================================
library(tidyverse)
library(ggplot2)
library(scales)

species_file <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_6/linux_data/transmission/filter_transmisson_gene_rpkm_taxa_agegroup.csv")
df_gene_sel <- species_file %>% select(Delivery,Species,U18AMPdb_name)%>%
  filter(!is.na(Species), Species != "")%>%
  filter(!is.na(Delivery), Delivery != "")


out_dir <- paste0(
  "C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/I-K--vaginal_vs_CS/L--gene_count_per_species/gene_count_new_final"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

df_gene_clean <- df_gene_sel %>%
  mutate(
    Species = as.character(Species),
    U18AMPdb_name = as.character(U18AMPdb_name),
    Delivery = case_when(
      str_to_lower(as.character(Delivery)) %in%
        c(
          "cs",
          "c-section",
          "c_section",
          "c section",
          "cesarean",
          "caesarean"
        ) ~ "C-section",
      
      str_to_lower(as.character(Delivery)) %in%
        c(
          "vaginal",
          "vaginal delivery"
        ) ~ "Vaginal",
      
      TRUE ~ as.character(Delivery)
    )
  ) %>%
  filter(
    !is.na(Species),
    Species != "",
    !is.na(U18AMPdb_name),
    U18AMPdb_name != "",
    Delivery %in% c("C-section", "Vaginal")
  )

loss_event_summary_species <- df_gene_clean %>%
  group_by(Species, Delivery) %>%
  summarise(
    gene_count_per_species = n_distinct(U18AMPdb_name),
    .groups = "drop"
  )

df_wide <- loss_event_summary_species %>%
  pivot_wider(
    id_cols = Species,
    names_from = Delivery,
    values_from = gene_count_per_species,
    values_fill = 0
  )

df_plot_top20 <- df_wide %>%
  mutate(
    total_count = `C-section` + Vaginal,
    diff_val = Vaginal - `C-section`
  ) %>%
  slice_max(
   order_by = diff_val,
   n = 20,
  with_ties = FALSE
  ) %>%
  arrange(
    desc(diff_val)
  ) %>%
  mutate(
    y_pos = rev(seq_len(n()))
  )

df_plot_top20 <- df_plot_top20 %>%
  mutate(
    C_section_plot = -`C-section`,
    Vaginal_plot = Vaginal
  )

write.csv(
  df_plot_top20,
  file.path(
    out_dir,
    "Top20_species_shared_AMPs_Csection_vs_Vaginal.csv"
  ),
  row.names = FALSE
)


max_value <- max(
  df_plot_top20$`C-section`,
  df_plot_top20$Vaginal,
  na.rm = TRUE
)

if (!is.finite(max_value) || max_value == 0) {
  max_value <- 1
}

x_pretty <- pretty(c(0, max_value), n = 5)
x_max <- max(x_pretty)

x_breaks <- pretty(
  c(-x_max, x_max),
  n = 8
)

x_breaks <- x_breaks[
  x_breaks >= -x_max &
    x_breaks <= x_max
]

x_breaks <- sort(unique(c(x_breaks, 0)))

n_species <- nrow(df_plot_top20)

y_breaks <- sort(df_plot_top20$y_pos)

y_labels <- df_plot_top20$Species[
  order(df_plot_top20$y_pos)
]

strip_ymin <- n_species + 0.55
strip_ymax <- n_species + 1.70
strip_ymid <- (strip_ymin + strip_ymax) / 2

strip_df <- tibble(
  xmin = c(-x_max, 0),
  xmax = c(0, x_max),
  ymin = strip_ymin,
  ymax = strip_ymax,
  x = c(-x_max / 2, x_max / 2),
  y = strip_ymid,
  label = c("C-section", "Vaginal")
)


p2 <- ggplot(
  df_plot_top20,
  aes(y = y_pos)
) +
  geom_segment(
    aes(
      x = C_section_plot,
      xend = Vaginal_plot,
      yend = y_pos
    ),
    color = "#AAA39C",
    linewidth = 0.55
  ) +
  
  # C-section
  geom_point(
    aes(x = C_section_plot),
    color = "#536F73",
    size = 2.8
  ) +
  
  # Vaginal
  geom_point(
    aes(x = Vaginal_plot),
    color = "#BC5B59",
    size = 2.8
  ) +
  
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.55
  ) +
  
  geom_rect(
    data = strip_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  geom_text(
    data = strip_df,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  
  scale_x_continuous(
    name = "No. of shared AMP genes",
    breaks = x_breaks,
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = y_labels,
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_cartesian(
    xlim = c(-x_max, x_max),
    ylim = c(0.5, n_species + 0.5),
    clip = "off"
  ) +
  
  labs(
    x = "No. of shared AMP clusters",
    y = NULL
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    panel.grid.major.x = element_line(
      color = "#E3DED8",
      linewidth = 0.4
    ),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.6
    ),
    
    axis.text.y = element_text(
      size = 8,
      color = "black",
      face = "italic",
      margin = margin(r = 4)
    ),
    
    axis.text.x = element_text(
      size = 9,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 11,
      color = "black",
      margin = margin(t = 8)
    ),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.4
    ),
    
    axis.ticks.length = unit(2, "pt"),
    
    legend.position = "none",
    
    plot.margin = margin(
      t = 35,
      r = 10,
      b = 10,
      l = 10
    ),
    plot.tag = element_text(
      size = 16,
      face = "bold"
    ),
    
    plot.tag.position = c(-0.10, 1.06)
  )





# ============================================================
#Species-AMP clusters association networks for Fig.5M
# ============================================================
library(tidyverse)
library(igraph)
library(ggraph)
library(patchwork)
library(scales)

in_file <- "C:/Users/taoch/Desktop/632751_clusters/Fig_6/linux_data/transmission/filter_transmisson_gene_rpkm_taxa_agegroup.csv"

out_dir <- "C:/Users/taoch/Desktop/632751_clusters/Fig_6/R_code/4.transmission/Dlivery_transmission_new/Species/"

top_n_species <- 10    
top_n_amp <- 50     
min_edge_RPKM <- 0  

df_raw <- read_csv(in_file, show_col_types = FALSE)

df <- df_raw %>%
  select(
    U18AMPdb_name,
    Species,
    RPKM,
    Delivery
  ) %>%
  filter(
    !is.na(U18AMPdb_name),
    !is.na(Species),
    !is.na(RPKM),
    !is.na(Delivery)
  ) %>%
  mutate(
    RPKM = as.numeric(RPKM),
    Delivery = case_when(
      str_to_lower(Delivery) %in% c("vaginal") ~ "Vaginal",
      str_to_lower(Delivery) %in% c("cs") ~ "C-section",
      TRUE ~ as.character(Delivery)
    )
  ) %>%
  filter(
    Delivery %in% c("Vaginal", "C-section"),
    RPKM > 0
  )

make_species_amp_network <- function(dat, delivery_group, plot_title) {
  
  dat_sub <- dat %>%
    filter(Delivery == delivery_group)
  
  if (nrow(dat_sub) == 0) {
    return(
      ggplot() +
        annotate(
          "text", x = 0, y = 0,
          label = paste0(plot_title, "/n no data"),
          size = 6
        ) +
        theme_void()
    )
  }
  
  
  keep_species <- dat_sub %>%
    group_by(Species) %>%
    summarise(
      Species_total_RPKM = sum(RPKM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Species_total_RPKM)) %>%
    slice_head(n = top_n_species) %>%
    pull(Species)

  keep_amp <- dat_sub %>%
    group_by(U18AMPdb_name) %>%
    summarise(
      AMP_total_RPKM = sum(RPKM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(AMP_total_RPKM)) %>%
    slice_head(n = top_n_amp) %>%
    pull(U18AMPdb_name)
  
  dat_sub <- dat_sub %>%
    filter(
      Species %in% keep_species,
      U18AMPdb_name %in% keep_amp
    )
  
  edge_df <- dat_sub %>%
    group_by(Species, U18AMPdb_name) %>%
    summarise(
      RPKM_sum = sum(RPKM, na.rm = TRUE),
      RPKM_mean = mean(RPKM, na.rm = TRUE),
      n_records = n(),
      .groups = "drop"
    ) %>%
    filter(RPKM_sum > min_edge_RPKM) %>%
    arrange(desc(RPKM_sum)) %>%
    rename(
      from = Species,
      to = U18AMPdb_name
    )
  
  if (nrow(edge_df) == 0) {
    return(
      ggplot() +
        annotate(
          "text", x = 0, y = 0,
          label = paste0(plot_title, "/nSpecies-AMP"),
          size = 6
        ) +
        theme_void()
    )
  }
  
  species_nodes <- dat_sub %>%
    group_by(Species) %>%
    summarise(
      total_RPKM = sum(RPKM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(Species %in% edge_df$from) %>%
    transmute(
      name = Species,
      type = "Species",
      total_RPKM = total_RPKM
    )
  
  amp_nodes <- dat_sub %>%
    group_by(U18AMPdb_name) %>%
    summarise(
      total_RPKM = sum(RPKM, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(U18AMPdb_name %in% edge_df$to) %>%
    transmute(
      name = U18AMPdb_name,
      type = "U18AMPdb",
      total_RPKM = total_RPKM
    )
  
  node_df <- bind_rows(species_nodes, amp_nodes) %>%
    mutate(
      log_total_RPKM = log10(total_RPKM + 1),
      
      node_size = case_when(
        type == "Species"  ~ 8.5,
        type == "U18AMPdb" ~ 5.8,
        TRUE ~ 5
      ),
      
      label_size = case_when(
        type == "Species"  ~ 3.4,
        type == "U18AMPdb" ~ 2.4,
        TRUE ~ 3
      ),
      
      label_fontface = case_when(
        type == "Species"  ~ "bold.italic",
        type == "U18AMPdb" ~ "plain",
        TRUE ~ "plain"
      )
    )
  
  g <- graph_from_data_frame(
    d = edge_df,
    vertices = node_df,
    directed = FALSE
  )
  
  set.seed(123)
  
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(edge_width = log10(RPKM_sum + 1)),
      colour = "#77AECD",
      alpha = 0.55,
      lineend = "round"
    ) +
    scale_edge_width(
      range = c(0.2, 2.8),
      name = "log10(RPKM + 1)"
    ) +
    
    geom_node_point(
      aes(
        filter = type == "Species",
        size = node_size,
        fill = type
      ),
      shape = 21,
      colour = "white",
      stroke = 0.45
    ) +
    geom_node_point(
      aes(
        filter = type == "U18AMPdb",
        size = node_size,
        fill = type
      ),
      shape = 23,
      colour = "white",
      stroke = 0.30
    ) +
    
    scale_fill_manual(
      values = c(
        "Species" = "#5E7A9C",
        "U18AMPdb" = "#C87271"
      )
    ) +
    scale_size_identity(
      guide = "none"
    ) +
    
    geom_node_text(
      aes(
        filter = type == "Species",
        label = name,
        size = label_size,
        fontface = label_fontface
      ),
      repel = TRUE,
      max.overlaps = Inf
    ) +
    
    geom_node_text(
      aes(
        filter = type == "U18AMPdb",
        label = name,
        size = label_size,
        fontface = label_fontface
      ),
      repel = TRUE,
      max.overlaps = Inf
    ) +
    
    scale_size_identity() +
    ggtitle(plot_title) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 22,
        face = "bold"
      ),
      legend.position = "right",
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

p_vaginal <- make_species_amp_network(
  dat = df,
  delivery_group = "Vaginal",
  plot_title = "Vaginal"
)

p_cs <- make_species_amp_network(
  dat = df,
  delivery_group = "C-section",
  plot_title = "C-section"
)

make_edge_node_table <- function(dat, delivery_group) {
  
  dat_sub <- dat %>%
    filter(Delivery == delivery_group)
  
  keep_species <- dat_sub %>%
    group_by(Species) %>%
    summarise(total_RPKM = sum(RPKM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_RPKM)) %>%
    slice_head(n = top_n_species) %>%
    pull(Species)
  
  keep_amp <- dat_sub %>%
    group_by(U18AMPdb_name) %>%
    summarise(total_RPKM = sum(RPKM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_RPKM)) %>%
    slice_head(n = top_n_amp) %>%
    pull(U18AMPdb_name)
  
  edge_df <- dat_sub %>%
    filter(
      Species %in% keep_species,
      U18AMPdb_name %in% keep_amp
    ) %>%
    group_by(Species, U18AMPdb_name) %>%
    summarise(
      RPKM_sum = sum(RPKM, na.rm = TRUE),
      RPKM_mean = mean(RPKM, na.rm = TRUE),
      n_records = n(),
      .groups = "drop"
    ) %>%
    filter(RPKM_sum > min_edge_RPKM) %>%
    arrange(desc(RPKM_sum))
  
  return(edge_df)
}

edge_vaginal <- make_edge_node_table(df, "Vaginal")
edge_cs <- make_edge_node_table(df, "C-section")

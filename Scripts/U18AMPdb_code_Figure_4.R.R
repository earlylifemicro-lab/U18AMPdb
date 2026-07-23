# ============================================================
#UPGMA clustering for Fig.4A
# ============================================================
#6122sample
library(readr)
library(dplyr)
library(stringr)
library(RColorBrewer)
#install.packages("ggsci")
library(ggsci)
library(itol.toolkit)
library(ape)
df <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_4CDEF_heatmap/gut_heatmap/gut_metadata_Zeng_2022.csv") 

library(ape)    
library(readr)   
library(dplyr)
library(itol.toolkit)
library(RColorBrewer)
tree_path <-"C:/Users/taoch/Desktop/632751_clusters/Fig_5/Fig_5A/6122_sample_amp_similarity_tree_gut.nwk" 
tree <- read.tree("C:/Users/taoch/Desktop/632751_clusters/Fig_5/Fig_5A/6122_sample_amp_similarity_tree_gut.nwk")
tree_tips <- ape::read.tree(tree_path)$tip.label
head(tree$tip.label)

hub <- create_hub(tree = tree_path)

gut_metadata_unique <-df %>%
  mutate(Sample = trimws(Sample)) 

Age_year_group_level_order <-gut_metadata_unique %>% dplyr::select(Sample,Age_year_group)
Delivery_level_order <-gut_metadata_unique%>% dplyr::select(Sample,Delivery)
Country_level_order <-gut_metadata_unique%>% dplyr::select(Sample,Country)
Gender_level_order <-gut_metadata_unique%>% dplyr::select(Sample,Gender)
Term_level_order <-gut_metadata_unique%>% dplyr::select(Sample,Term)
Feed_level_order <-gut_metadata_unique%>% dplyr::select(Sample,Feed)


metadata_Age_year_group <- Age_year_group_level_order %>%
  filter(Sample %in% tree_tips) %>%rename(id = Sample)%>%
  mutate(Age_year_group = ifelse(is.na(Age_year_group), "Unknown", Age_year_group))
metadata_Delivery <- Delivery_level_order %>%
  filter(Sample  %in% tree_tips) %>%rename(id = Sample)%>%
  mutate(Delivery = ifelse(is.na(Delivery), "Unknown", Delivery))
metadata_Country <- Country_level_order %>%
  filter(Sample  %in% tree_tips) %>%rename(id = Sample)%>%
  mutate(Country = ifelse(is.na(Country), "Unknown", Country))
metadata_Gender <- Gender_level_order %>%
  filter(Sample  %in% tree_tips) %>%rename(id = Sample)%>%
  mutate(Gender = ifelse(is.na(Gender), "Unknown", Gender))
metadata_Term <- Term_level_order %>%
  filter(Sample  %in% tree_tips) %>%rename(id = Sample)%>%
  mutate(Term = ifelse(is.na(Term), "Unknown",Term))
metadata_Feed <- Feed_level_order %>%
  filter(Sample  %in% tree_tips) %>%rename(id = Sample)%>%
  mutate(Feed = ifelse(is.na(Feed), "Unknown",Feed))

#Age_year_group
unit_1 <- create_unit(data = metadata_Age_year_group,
                      key = "metadata_Age_year_group", 
                      type = "DATASET_COLORSTRIP",
                      tree = tree)
groups <- unique(metadata_Age_year_group$Age_year_group)
n_groups <- length(groups)
print(groups)

my_colors <- brewer.pal(n = n_groups, name = "Set1")
names(my_colors) <- groups

unit_1@field$colors <- my_colors
print(unit_1@field$colors)

#Delivery
unit_2 <- create_unit(data = metadata_Delivery,
                      key = "metadata_Delivery", 
                      type = "DATASET_COLORSTRIP",
                      tree = tree)
groups <- unique(metadata_Delivery$Delivery)
n_groups <- length(groups)
print(groups)

my_colors <- brewer.pal(n = n_groups, name = "Set1")
names(my_colors) <- groups

unit_2@field$colors <- my_colors
print(unit_2@field$colors)


#Country
unit_3 <- create_unit(data = metadata_Country,
                      key = "metadata_Country", 
                      type = "DATASET_COLORSTRIP",
                      tree = tree)
groups <- unique(metadata_Country$Country)
n_groups <- length(groups)
print(groups)

my_colors <- pal_d3("category20", alpha = 1)(n_groups)
names(my_colors) <- groups

unit_3@field$colors <- my_colors
print(unit_3@field$colors)

#Gender
unit_4 <- create_unit(data = metadata_Gender,
                      key = "metadata_Gender", 
                      type = "DATASET_COLORSTRIP",
                      tree = tree)
groups <- unique(metadata_Gender$Gender)
n_groups <- length(groups)
print(groups)

my_colors <- brewer.pal(n = n_groups, name = "Set1")
names(my_colors) <- groups

unit_4@field$colors <- my_colors
print(unit_4@field$colors)

#Term
unit_5 <- create_unit(data = metadata_Term,
                      key = "metadata_Term", 
                      type = "DATASET_COLORSTRIP",
                      tree = tree)
groups <- unique(metadata_Term$Term)
n_groups <- length(groups)
print(groups)

my_colors <- brewer.pal(n = n_groups, name = "Set1")
names(my_colors) <- groups

unit_5@field$colors <- my_colors
print(unit_5@field$colors)

#Feed
unit_6 <- create_unit(data = metadata_Feed,
                      key = "metadata_Feed", 
                      type = "DATASET_COLORSTRIP",
                      tree = tree)

groups <- unique(metadata_Feed$Feed)
n_groups <- length(groups)
print(groups)

my_colors <- brewer.pal(n = n_groups, name = "Set1")
names(my_colors) <- groups

unit_6@field$colors <- my_colors
print(unit_6@field$colors)


## add unit into hub
hub <- hub+unit_1+unit_2+unit_3+unit_4+unit_5+unit_6

## write template file
write_hub(hub,dir="C:/Users/taoch/Desktop/632751_clusters/Fig_5/Fig_5A/")



# ============================================================
#PCoA of AMP cluster profiles for Fig.4BC; Fig.S9A-D
# ============================================================

library(ggplot2)
library(cowplot)
library(vegan)


pcoa_df <- data.frame(
  Sample = rownames(pcoa_res$points),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2]
)

plot_df <- merge(
  pcoa_df,
  metadata_sub,
  by = "Sample",
  all = FALSE
)


any(duplicated(colnames(plot_df)))
colnames(plot_df)[duplicated(colnames(plot_df))]  #character(0)

head(plot_df)


var_exp <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 2)  
xlab <- paste0("PCoA1 (", var_exp[1], "%)")  #16.45%
ylab <- paste0("PCoA2 (", var_exp[2], "%)")  #13.29%

#Age_year_group
plot_df <-plot_df[!is.na(plot_df$Age_year_group), ] 

age_cols <- c(
  "0-1 years"="#3FA7A3",  
  "1-2 years"="#E76F51", 
  "2-3 years" = "#B83945"
)

x_lim <- range(plot_df$PCoA1, na.rm = TRUE)
y_lim <- range(plot_df$PCoA2, na.rm = TRUE)


p_scatter <- ggplot(
  plot_df,
  aes(PCoA1, PCoA2, color = Age_year_group)
) +
  geom_point(size = 1, alpha = 0.5) +
  coord_fixed(
    xlim = x_lim,
    ylim = y_lim,
    expand = FALSE
  ) +
  scale_color_manual(values = age_cols, drop = FALSE) +
  theme_classic(base_size = 15) +
  theme(
    axis.line = element_line(size = 0.8),
    axis.ticks = element_line(size = 0.8),
    legend.position = "none"
  ) +
  labs(
    x = "PCoA1 (16.45%)",
    y = "PCoA2 (13.29%)"
  )


p_top <- ggplot(
  plot_df,
  aes(x = PCoA1, fill = Age_year_group, color = Age_year_group)
) +
  geom_density(alpha = 0.35, linewidth = 1) +
  scale_fill_manual(values = age_cols) +
  scale_color_manual(values = age_cols) +
  coord_cartesian(
    xlim = x_lim,
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  )


p_right <- ggplot(
  plot_df,
  aes(y = PCoA2, fill = Age_year_group, color = Age_year_group)
) +
  geom_density(alpha = 0.35, linewidth = 1) +
  scale_fill_manual(values = age_cols) +
  scale_color_manual(values = age_cols) +
  coord_cartesian(
    ylim = y_lim,
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  )


legend_plot <- ggplot(
  plot_df,
  aes(PCoA1, PCoA2, color = Age_year_group)
) +
  geom_point() +
  scale_color_manual(values = age_cols) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_blank()
  )

legend <- get_legend(legend_plot)


main_plot <- plot_grid(
  p_top, NULL,
  p_scatter, p_right,
  ncol = 2,
  rel_widths = c(1, 0.18), 
  rel_heights = c(0.18, 1) 
)

final_plot <- plot_grid(
  main_plot,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.25)
)

final_plot




# ============================================================
#Age-associated change in AMP cluster for Fig.4D; Fig.S9E-F
# ============================================================
library(ggplot2)
library(tidyverse)

data <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_5/Fig_5DEFGH/linux_data/shannon_meta_2538_sample.csv", sep=",") 

data$Age_month <- readr::parse_number(as.character(data$Age_month))

#（散点图 + 趋势线）
p1 <- ggplot(data, aes(x = Age_month_group, y = amp_shannon)) +
  geom_point(aes(color = amp_shannon), size = 2, alpha = 0.6) +
  scale_color_gradient(low = "#299d8f", 
                       high = "#f3a361", 
                       name = "shannon index",
                       limits = c(0, 12),
                       breaks = seq(0, 12, 3) 
  )+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black", se = TRUE, level = 0.99, linewidth = 0.9) +  
 
  labs(
    x = "Age (months)",
    y = "Shanno diversity"
  ) +
  scale_x_continuous(breaks = c(0,1, 3, 6, 12, 18, 24, 30, 36))+
  scale_y_continuous(breaks = seq(0,12,3),
                     limits = c(0,12))+
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
p1


fit <- lm(amp_shannon ~ Age_month_group, data = data)
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
p2




# =========================================================================
#Age-stratified effects of host and environmental covariates for Fig.4E
# ========================================================================

library(tidyverse)

df <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_5/Fig_5DEFGH/R_code/Fig_5H/adonis_Age_others5_new.csv",header = TRUE, sep = ",")

df <- df %>%
  mutate(
    Time = factor(time_bin, levels = c(0,1,3,6,12,18,24,30,36)),
    Variable = recode(variable,
                      Delivery = "Delivery mode",
                      Country  = "Country",
                      Gender   = "Gender",
                      Term     = "Term",
                      Feed     = "Feeding pattern"
    ),
    signif = FDR < 0.05
  )

df$Variable <- factor(
  df$Variable,
  levels = c(
    "Country",
    "Feeding pattern",
    "Term",
    "Delivery mode",
    "Gender"              
  )
)


var_cols <- c(
  "Country"              = "#41AB5D",
  "Feeding pattern"      = "#7B3294",
  "Term"                 = "#2EC4B6",
  "Delivery mode"       = "#2C7FB8",
  "Gender"               = "#F1B722"
)


library(ggplot2)

p <- ggplot(df, aes(x = Time, y = Variable)) +
  
  # R² 
  geom_point(
    aes(size = R2, fill = Variable),
    shape = 21,
    color = "black",
    alpha = 0.85
  ) +
  
  # FDR < 0.05
  geom_tile(
    data = df %>% filter(signif),
    fill = NA,
    color = "red",
    linewidth = 0.7
  ) +
  
  scale_fill_manual(values = var_cols) +
  
  scale_size_continuous(
    range = c(2, 10),
    breaks = c(0.05, 0.10, 0.15, 0.20, 0.25),
    name = expression(R^2)
    #limits = c(0.05, 0.25) 
  ) +
  
  labs(
    x = "Timepoints (month)",
    y = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA)
  )

p





# ========================================================================
#Divergent age-associated trajectories for Fig.4F;Fig.S10
# ========================================================================
#Fig.4F

library(tidyverse)
library(ggplot2)

df_dec <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/significant_decreasing_raw_data.csv",
  header = TRUE
)

df_inc <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/significant_increasing_raw_data.csv",
  header = TRUE
)


age_levels <- sort(unique(c(df_dec$Age_month_group, df_inc$Age_month_group)))  

df_dec$Age_month_group <- factor(df_dec$Age_month_group, levels = age_levels)
df_inc$Age_month_group <- factor(df_inc$Age_month_group, levels = age_levels)


set.seed(123) 

sample_n_by_group <- function(df, n = 400) {
  df %>%
    group_by(Age_month_group) %>%
    sample_n(size = min(n, n()), replace = FALSE) %>%   
}

df_dec_sampled <- sample_n_by_group(df_dec, n = 400)
df_inc_sampled <- sample_n_by_group(df_inc, n = 400)


df_dec_sampled$type <- "Decreasing"
df_inc_sampled$type <- "Increasing"


df_combined <- bind_rows(df_dec_sampled, df_inc_sampled)


p <- ggplot(df_combined, aes(x = Age_month_group, y = log10_RPKM, color = type)) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    alpha = 0.4
  ) +
  stat_smooth(
    aes(group = type),
    method = "loess",
    se = TRUE,
    linewidth = 1.2
  ) +
  scale_color_manual(values = c("Decreasing" = "#6FA19F","Increasing" = "#b38287"))+
  theme_classic(base_size = 15) +
  labs(
    x = "Age (month group)",
    y = expression(log[10]~"(RPKM)"),
    title = "Abundance trends with age: Decreasing vs Increasing",
    color = "Trend Type"
  )

ggsave(plot = p, filename = "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/R_code/trends_loess.png", width = 10, height = 8, dpi = 300)
ggsave(plot = p, filename = "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/R_code/trends_loess.pdf", width = 10, height = 8)



#Fig.S10  increasing_top20_pheatmap
library(tidyverse)
library(pheatmap)

df_dec <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/significant_increasing_raw_data.csv",
  stringsAsFactors = FALSE
)


df_dec_scored <- df_dec %>%
  mutate(
    score = -log10(fdr_spearman) * abs(spearman_rho)
  )

top20_amp <- df_dec_scored %>%
  distinct(U18AMPdb_name, score) %>%
  arrange(desc(score)) %>%
  slice_head(n = 20) %>%
  pull(U18AMPdb_name)

df_top20 <- df_dec_scored %>%
  filter(U18AMPdb_name %in% top20_amp)

df_top20_sum <- df_top20 %>%
  group_by(U18AMPdb_name, Age_month_group) %>%
  summarise(
    log10_RPKM = median(log10_RPKM, na.rm = TRUE),
    .groups = "drop"
  )


cols_order <- c(0, 1, 3, 6, 12, 18, 24, 30, 36)

heatmap_increasing_data <- df_top20_sum %>%
  mutate(Age_month_group = as.character(Age_month_group)) %>%
  pivot_wider(
    names_from  = Age_month_group,
    values_from = log10_RPKM,
    values_fill = NA
  ) %>%
  select(U18AMPdb_name, any_of(as.character(cols_order)))


mat <- heatmap_increasing_data %>%
  column_to_rownames("U18AMPdb_name") %>%
  as.matrix()

storage.mode(mat) <- "numeric"


my_breaks <- seq(0, 5, length.out = 51)


my_legend_breaks <- seq(0, 5, by = 1) 

# PNG
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  color = colorRampPalette(
    c("#8FB4BE", "#AFC9CF", "#D5E1E3", "#EBBFC2", "#E28187", "#D93F49")
  )(50),
  filename = "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/R_code/increasing_top20_pheatmap.png",
  width = 10,
  height = 8
)

# PDF
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  color = colorRampPalette(
    c("#8FB4BE", "#AFC9CF", "#D5E1E3", "#EBBFC2", "#E28187", "#D93F49")
  )(50),
  filename = "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/R_code/increasing_top20_pheatmap.pdf",
  width = 10,
  height = 8
)




#Fig.S10  increasing_top20_lm
library(tidyverse)

df_dec <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/significant_increasing_raw_data.csv",
  stringsAsFactors = FALSE
)

df_dec_scored <- df_dec %>%
  mutate(
    score = -log10(fdr_spearman) * abs(spearman_rho)
  )


custom_order <- c(
  "U18AMPdb_013256",
  "U18AMPdb_413900",
  "U18AMPdb_554734",
  "U18AMPdb_482549",
  "U18AMPdb_577504",
  "U18AMPdb_155485",
  "U18AMPdb_561083",
  "U18AMPdb_339095",
  "U18AMPdb_610112",
  "U18AMPdb_249068",
  "U18AMPdb_591792",
  "U18AMPdb_084115",
  "U18AMPdb_130375",
  "U18AMPdb_516827",
  "U18AMPdb_205970",
  "U18AMPdb_344098",
  "U18AMPdb_081828",
  "U18AMPdb_408415",
  "U18AMPdb_461560",
  "U18AMPdb_464989"
)


df_top20 <- df_dec_scored %>%
  filter(U18AMPdb_name %in% custom_order) %>%
  mutate(U18AMPdb_name = factor(U18AMPdb_name, levels = custom_order)) %>%
  arrange(U18AMPdb_name)


p <- ggplot(df_top20, aes(x = Age_month_group, y = log10_RPKM)) +
  geom_point(size = 0.8, alpha = 0.5, color = "#c1bebf") +
  geom_smooth(method = "lm", se = TRUE, color = "#b38287", linewidth = 0.8) +
  facet_wrap(~ U18AMPdb_name, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = c(0, 1, 3, 6, 12, 18, 24, 30, 36)) +  
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA),
    axis.ticks = element_blank()
  ) +  
  labs(
    x = "Age (month)",
    y = expression(log[10](RPKM + 1))
  )


ggsave(
  plot = p,
  filename = "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/R_code/increasing_custom_order.png",
  width = 10,
  height = 8,
  dpi = 300
)


write.csv(
  df_top20,
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGH_2538_sample/above_3_time/cluster_per_trend/R_code/increasing_custom_order_table.csv",
  row.names = FALSE
)




# ====================================================================================================
#Amino acid composition of increasing and decreasing age-associated AMP clusters for Fig.4G;Fig.S11A-B
# ====================================================================================================
#Fig.4G

library(tidyverse)
library(pheatmap)

aa_dec <- read.csv(
  "C:/Users/taoch/Desktop/632751_clusters_文章/Fig_5/Fig_5J_property/AA_frequency_decreasing.csv",
  check.names = FALSE
)

aa_inc <- read.csv(
  "C:/Users/taoch/Desktop/632751_clusters_文章/Fig_5/Fig_5J_property/AA_frequency_increasing.csv",
  check.names = FALSE
)

aa_dec_long <- aa_dec %>%
  pivot_longer(
    cols = starts_with("freq_"),
    names_to = "AA",
    values_to = "Decreasing"
  ) %>%
  mutate(
    AA = str_remove(AA, "^freq_")
  ) %>%
  select(AA, Decreasing)

aa_inc_long <- aa_inc %>%
  pivot_longer(
    cols = starts_with("freq_"),
    names_to = "AA",
    values_to = "Increasing"
  ) %>%
  mutate(
    AA = str_remove(AA, "^freq_")
  ) %>%
  select(AA, Increasing)

aa_compare <- aa_inc_long %>%
  left_join(aa_dec_long, by = "AA") %>%
  select(AA, Increasing, Decreasing)

write.csv(
  aa_compare,
  "AA_frequency_Increasing_and_Decreasing.csv",
  row.names = FALSE
)

aa_mat <- aa_compare %>%
  column_to_rownames("AA") %>%
  as.matrix()

aa_mat <- aa_mat[, c("Increasing", "Decreasing")]

pheatmap(
  aa_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  display_numbers = TRUE,
  number_format = "%.3f",
  main = "Amino acid composition: Increasing vs Decreasing",
  angle_col = 45
)



#Fig.S11A-B RPKM—weighted
library(tidyverse)
library(ggplot2)
library(ggsci)


sel <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGHIJK_2538_sample/above_3_time/cluster_per_trend/significant_decreasing_raw_data.csv",
  stringsAsFactors = FALSE
)


data <- sel %>%
  select(U18AMPdb_name, Sequence, RPKM, Age_month_group) %>%
  mutate(
    Age_month_group = as.character(Age_month_group),
    RPKM = as.numeric(RPKM)
  ) %>%
  filter(!is.na(Sequence), !is.na(RPKM), RPKM > 0)


aa_letters <- c(
  "A","C","D","E","F","G","H","I","K","L",
  "M","N","P","Q","R","S","T","V","W","Y"
)

calc_aa_freq <- function(seq){
  seq <- toupper(seq)
  n <- nchar(seq)

  sapply(aa_letters, function(aa){
    stringr::str_count(seq, aa) / n
  })
}

aa_freq <- data %>%
  distinct(U18AMPdb_name, Sequence) %>%
  mutate(
    aa_freq = map(Sequence, calc_aa_freq)
  ) %>%
  unnest_wider(aa_freq, names_sep = "_") %>%
  rename_with(
    ~ paste0("freq_", aa_letters),
    starts_with("aa_freq_")
  )


data_aa <- data %>%
  left_join(
    aa_freq,
    by = c("U18AMPdb_name", "Sequence")
  )

aa_cols <- paste0("freq_", aa_letters)


weighted_aa_by_age <- data_aa %>%
  group_by(Age_month_group) %>%
  summarise(
    total_RPKM = sum(RPKM, na.rm = TRUE),
    n_AMP = n_distinct(U18AMPdb_name),
    across(
      all_of(aa_cols),
      ~ sum(.x * RPKM, na.rm = TRUE) / sum(RPKM, na.rm = TRUE)
    ),
    .groups = "drop"
  )

write.csv(
  weighted_aa_by_age,
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGHIJK_2538_sample/above_3_time/cluster_property/RPKM_weighted_AA_frequency_by_Age_month_group_decreasing.csv",
  row.names = FALSE
)


age_order <- c("0", "1", "3", "6", "12", "18", "24", "30", "36")

weighted_aa_long <- weighted_aa_by_age %>%
  pivot_longer(
    cols = all_of(aa_cols),
    names_to = "Amino_acid",
    values_to = "Weighted_fraction"
  ) %>%
  mutate(
    Amino_acid = str_remove(Amino_acid, "^freq_"),
    Amino_acid = factor(Amino_acid, levels = aa_letters),
    Age_month_group = factor(Age_month_group, levels = age_order)
  ) %>%
  arrange(Amino_acid, Age_month_group)

colors_simpsons <- pal_simpsons(alpha = 0.8)(
  nlevels(weighted_aa_long$Age_month_group)
)


p <- ggplot(
  weighted_aa_long,
  aes(
    x = Amino_acid,
    y = Weighted_fraction,
    fill = Age_month_group
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 1.05),
    width = 0.40
  ) +
  labs(
    title = "RPKM-weighted amino acid frequency by age month group",
    x = "Amino acids",
    y = "RPKM-weighted fraction",
    fill = "Age (months)"
  ) +
  scale_fill_manual(values = colors_simpsons) +
  scale_y_continuous(
    breaks = seq(0, 0.15, by = 0.025)
  ) +
  coord_cartesian(ylim = c(0, 0.15)) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5)
  )

print(p)

ggsave(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGHIJK_2538_sample/above_3_time/cluster_property/RPKM_weighted_AA_frequency_by_Age_month_group_decreasing.pdf",
  p,
  width = 18,
  height = 7
)

ggsave(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_5DEFGHIJK_2538_sample/above_3_time/cluster_property/RPKM_weighted_AA_frequency_by_Age_month_group_decreasing.png",
  p,
  width = 18,
  height = 7,
  dpi = 300
)





# =========================================================================
#Three-dimensional physicochemical properity for Fig.4H;Fig.S11C
# =========================================================================
#Fig.4H
#x：Charge
#y：HydrophobicMoment_uH
#z：hydrophobicity_h

sample_3d <- sample_weighted_property %>%
  select(
    Sample,
    Age_month_group,
    trend_direction,
    total_RPKM,
    w_Charge,
    w_HydrophobicMoment_uH,
    w_hydrophobicity_h
  ) %>%
  drop_na()

color_map <- c(
  "Increasing" = "#b38287",
  "Decreasing" = "#6FA19F"
)


p <- plot_ly()

for(grp in c("Increasing", "Decreasing")){
  
  tmp <- sample_3d %>%
    filter(trend_direction == grp)
  
  p <- p %>%
    add_trace(
      data = tmp,
      x = ~w_Charge,
      y = ~w_HydrophobicMoment_uH,
      z = ~w_hydrophobicity_h,
      
      type = "scatter3d",
      mode = "markers",
      name = grp,
      
      marker = list(
        size = log10(tmp$total_RPKM + 1),
        color = color_map[[grp]],
        opacity = 0.65,
        showscale = FALSE,
        line = list(
          color = color_map[[grp]],
          width = 0
        )
      ),
      
      text = ~paste(
        "Trend =", trend_direction,
        "<br>Age =", Age_month_group,
        "<br>RPKM =", round(total_RPKM, 2)
      ),
      
      hoverinfo = "text",
      showlegend = TRUE
    )
}

p <- p %>%
  layout(
    scene = list(
      xaxis = list(title = "Net charge"),
      yaxis = list(title = "Hydrophobic moment (uH)"),
      zaxis = list(title = "Global hydrophobicity (h)")
    ),
    legend = list(
      title = list(text = "Trend")
    )
  )

p



#Fig.S11C
library(tidyverse)
library(ggpubr)
library(broom)


dec_sel <- read.csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_5/Fig_5J_property/significant_decreasing_selected.csv", stringsAsFactors = FALSE)
inc_sel <- read.csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_5/Fig_5J_property/significant_increasing_selected.csv", stringsAsFactors = FALSE)

dec_clu <- read.csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_5/Fig_5J_property/significant_decreasing_cluster.csv", stringsAsFactors = FALSE)
inc_clu <- read.csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_5/Fig_5J_property/significant_increasing_cluster.csv", stringsAsFactors = FALSE)


dec_data <- dec_sel %>%
  left_join(dec_clu, by = "U18AMPdb_name") %>%
  mutate(trend_direction = "Decreasing")

inc_data <- inc_sel %>%
  left_join(inc_clu, by = "U18AMPdb_name") %>%
  mutate(trend_direction = "Increasing")

plot_data <- bind_rows(dec_data, inc_data) %>%
  mutate(
    Age_month_group = as.numeric(Age_month_group),
    trend_direction = factor(trend_direction, levels = c("Increasing", "Decreasing"))
  )


w_mean <- function(x, w) {
  if(all(is.na(x)) || sum(w, na.rm = TRUE) == 0) return(NA_real_)
  sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

property_cols <- c(
  "Length", "MW", "Charge", "ChargeDensity", "pI",
  "InstabilityInd", "Aromaticity", "AliphaticInd",
  "BomanInd", "HydrophRatio",
  "HydrophobicMoment_uH", "hydrophobicity_h"
)

sample_weighted_property <- plot_data %>%
  group_by(Sample, 
           Age_month_group,  
           trend_direction) %>%  #每个样本内部，把Increasing AMP和Decreasing AMP分开计算
  summarise(
    total_RPKM = sum(10^log10_RPKM, na.rm = TRUE), #恢复原始丰度
    n_AMP = n_distinct(U18AMPdb_name),   #该Sample该Trend下检测到多少种AMP
    across(
      all_of(property_cols),
      ~ w_mean(.x, 10^log10_RPKM),
      .names = "w_{.col}"            #计算加权性质
    ),
    .groups = "drop"
  )
write.csv(sample_weighted_property,"C:/Users/taoch/Desktop/632751_clusters_文章/Fig_5/Fig_5J_property/R_code/sample_weighted_property.csv",row.names = FALSE)


stat_spearman <- property_long %>%
  group_by(trend_direction, Property) %>%
  summarise(
    spearman_rho = {
      tmp <- na.omit(data.frame(Age_month_group, Weighted_value))
      cor(tmp$Age_month_group, tmp$Weighted_value,
          method = "spearman")
    },
    p_value = {
      tmp <- na.omit(data.frame(Age_month_group, Weighted_value))
      cor.test(
        tmp$Age_month_group,
        tmp$Weighted_value,
        method = "spearman",
        exact = FALSE 
      )$p.value
    },
    .groups = "drop"
  ) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"), 
    significant = fdr < 0.05
  )

write.csv(stat_spearman,"C:/Users/taoch/Desktop/632751_clusters_文章/Fig_5/Fig_5J_property/R_code/stat_spearman.csv",row.names = FALSE)


library(tidyverse)
library(ggpubr)


stat_label <- stat_spearman %>%
  mutate(
    label = paste0(
      "rho = ", round(spearman_rho, 3),
      "\nFDR = ", signif(fdr, 3)
    )
  )


label_pos <- property_long %>%
  group_by(Property) %>%
  summarise(
    x_pos = max(Age_month_group, na.rm = TRUE),
    y_max = max(Weighted_value, na.rm = TRUE),
    y_min = min(Weighted_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  crossing(
    trend_direction = c("Increasing", "Decreasing")
  ) %>%
  mutate(
    y_pos = ifelse(
      trend_direction == "Increasing",
      y_max,
      y_min
    ),
    vjust_value = ifelse(
      trend_direction == "Increasing",
      1,
      0
    )
  )

stat_label <- stat_label %>%
  left_join(
    label_pos,
    by = c("Property", "trend_direction")
  )


p_B <- ggplot(
  property_long,
  aes(
    x = Age_month_group,
    y = Weighted_value,
    color = trend_direction,
    fill = trend_direction
  )
) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.9
  ) +
  geom_text(
    data = stat_label,
    aes(
      x = x_pos,
      y = y_pos,
      label = label,
      color = trend_direction,
      vjust = vjust_value
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 3
  ) +
  scale_color_manual(
    values = c(
      "Increasing" = "#b38287",
      "Decreasing" = "#6FA19F"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Increasing" = "#b38287",
      "Decreasing" = "#6FA19F"
    )
  ) +
  facet_wrap(
    ~ Property,
    scales = "free_y",
    ncol = 4
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 3, 6, 12, 18, 24, 30, 36)
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "Age (months)",
    y = "RPKM-weighted physicochemical property",
    color = "Trend",
    fill = "Trend"
  ) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 10),
    legend.position = "top"
  )

print(p_B)


ggsave(
  "B_RPKM_weighted_property_lm_with_spearman.pdf",
  p_B,
  width = 12,
  height = 8
)

ggsave(
  "B_RPKM_weighted_property_lm_with_spearman.png",
  p_B,
  width = 12,
  height = 8,
  dpi = 300
)
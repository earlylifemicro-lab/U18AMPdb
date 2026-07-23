# ============================================================
#Body-site distribution of U18AMPdb for Fig.2A-B; Fig.S3
# ============================================================
##log scale for Fig.S3A

library(ggplot2)
library(dplyr)
library(tidyverse)
sizes <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_3AB/linux_data/632751_cluster_summary.csv")

# 1.log scale for Fig.S3A
ggplot(sizes, aes(x = AMP_counts)) +
  geom_histogram(bins = 50, fill = "#a30543") +
  scale_x_log10() +
  labs(x = "Number of AMPs (log10 scale)",
       y = "Count of AMP clusters") +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black",linewidth=0.5),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.text.x = element_text(size=12, colour = "black"),
        panel.grid = element_blank())



# 2. Cumulative distribution for Fig.S3B
sizes_sorted <- sizes %>%
  arrange(desc(AMP_counts)) %>%
  mutate(
    cumsum = cumsum(AMP_counts),
    cumprop = cumsum / sum(AMP_counts),
    rank = row_number()  
  )

idx75 <- which(sizes_sorted$cumprop >= 0.75)[1]
cat("Top", idx75, "clusters cover 75% of all sequences\n")
##Top 66503 clusters cover 75% of all sequences

ggplot(sizes_sorted, aes(x = rank, y = cumprop)) +  
  geom_line(color = "#4DB6AC") +
  geom_vline(xintercept = idx75, linetype = "dashed", color = "red") +
  labs(
    x = "Cluster Rank (largest first)",
    y = "Cumulative Proportion of Sequences"
  ) +
  theme_minimal()


ggplot(sizes_sorted, aes(x = rank, y = cumprop)) +
  geom_line(color = "#4DB6AC") +
  geom_vline(xintercept = idx75, linetype = "dashed", color = "red") +
  geom_point(data = sizes_sorted[idx75, ], aes(x = rank, y = cumprop), 
             color = "red", size = 2) +
  annotate("text", x = idx75, y = 0.82, label = paste("Top", idx75), 
           color = "red", hjust = 0) +
  labs(
    x = " AMP cluster rank (largest first)",
    y = "Cumulative proportion of AMPs"
  ) +
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

##6*4


# 3. Top-50 AMP clusters for Fig.S3C
top50 <- sizes_sorted[1:50, ]

ggplot(top50, aes(x = reorder(U18AMPdb_name, AMP_counts), y = AMP_counts)) +
  geom_col(fill = "#4DB6AC") +
  coord_flip() +
  labs(title = "Top 50 largest AMP clusters",
       x= "AMP cluster",
       y = "Number of AMPs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank(),
  )



#4 Log-logarithmic scatter plot (power-law distribution) for Fig.2B
library(ggplot2)

ggplot(sizes_sorted, aes(x = log(rank), y = log(AMP_counts))) +
  geom_point(alpha = 0.6, color = "#4DB6AC") +
  geom_smooth(method = "lm", se = FALSE, color = "#a30543") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("α = ", round(alpha, 2)),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  labs(
    x = "Log (Rank of AMP clusters)",
    y = "Log (AMPs counts)"
  ) +
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black",linewidth=0.5))


model <- lm(log(AMP_counts) ~ log(rank), data = sizes_sorted)
alpha <- abs(coef(model)[2])  # 斜率的绝对值
cat("Power law exponent α =", round(alpha, 2), "\n")
#α = 1.24

write.csv(sizes_sorted,"C:/Users/taoch/Desktop/632751_clusters/Fig_3AB/linux_data/632751_cluster_summary_source_data.csv")


#5 UpSet plot for Fig.2A
library(tidyverse)
library(UpSetR)
library(cowplot)
library(grid)
library(gridExtra) 


df <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_3AB/linux_data/632751_cluster_body_site_multisite.csv")
head(df)


upset_data <- df%>%
  select(cluster_id, body_site) %>%
  distinct() %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = body_site, values_from = value, values_fill = FALSE) %>%
  column_to_rownames("cluster_id")


df_binary <- upset_data %>%
  mutate(across(everything(), ~ ifelse(. == TRUE, 1, 0)))


p1 <- UpSetR::upset (
  df_binary,
  nsets = 5,
  order.by = "freq", 
  # 1. 设置X轴顺序：milk, nasal, oral, skin, gut
  sets = c("milk", "nasal", "oral", "skin", "gut"),  
  keep.order = TRUE, 
  decreasing = TRUE,
  mb.ratio = c(0.7, 0.3), 
  number.angles = 20,
  point.size = 3, 
  line.size = 1, 
  mainbar.y.label = "Number of AMP clusters", 
  sets.x.label = "Body site", 
  main.bar.color="#299d8f",
  matrix.color = "#274753",
  sets.bar.color = c("#16058b", "#8ab07c", "#e66d50", "#4DB6AC", "#a30543"),  
  text.scale = c(1.25, 1.25, 1.25, 1.25, 1.25, 1) 
)
print(p1)


p2 <- UpSetR::upset(
  df_binary,
  nsets = 5,
  order.by = "freq",
  sets = c("gut", "skin", "oral", "nasal","milk"),
  keep.order = TRUE, 
  decreasing = TRUE,
  
  # 可视化参数
  mb.ratio = c(0.7, 0.3),
  number.angles = 20,
  point.size = 3,
  line.size = 1,
  mainbar.y.label = "Number of AMP clusters",
  sets.x.label = "Site-specific AMP clusters",
  main.bar.color = "#299d8f",
  matrix.color = "#274753",
  

  sets.bar.color = c("#a30543","#4DB6AC","#e66d50","#8ab07c","#16058b"),
  queries = list(
    list(query = intersects, params = list("milk"), color = "#16058b", active = FALSE),
    list(query = intersects, params = list("nasal"), color = "#8ab07c", active = FALSE),
    list(query = intersects, params = list("oral"), color = "#e66d50", active = FALSE),
    list(query = intersects, params = list("skin"), color = "#4DB6AC", active = FALSE),
    list(query = intersects, params = list("gut"), color = "#a30543", active = FALSE)
  ),
  
  text.scale = c(1.25, 1.25, 1.25, 1.25, 1.25, 1)
)
print(p2)







# ============================================================
# prevalence and mean abundance of U18AMPdb for Fig.2C;Fig.S4
# ============================================================
#Scatter_plot 

library(tidyverse)

stats <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/coverM/all_11958sample/all_11958sample_stats_with_U18AMPdb_name.csv",
  check.names = FALSE
)

out_dir <- "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/coverM/all_11958sample"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


site_order <- c("milk", "nasal", "oral", "skin", "gut")

plot_data <- stats %>%
  mutate(
    body_site = str_extract(Cluster_name, "_(gut|skin|oral|nasal|milk)$"),
    body_site = str_remove(body_site, "^_"),
    body_site = replace_na(body_site, "unknown")
  ) %>%
  filter(body_site != "unknown") %>%
  mutate(
    body_site = factor(body_site, levels = site_order)
  )

print(table(plot_data$body_site))

#  milk  nasal   oral   skin    gut
#  1802  22373 313088 186246  96726


p2_coverm_all_11958sample_Scatter_plot <- ggplot(
  plot_data,
  aes(
    x = prevalence,
    y = log10_mean_abundance,
    color = body_site
  )
) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, by = 2)
  ) +
  labs(
    x = "Prevalence (proportion of samples with >0 abundance)",
    y = "Mean abundance (log10(rpkm + 1))",
    color = "Niche of representative"
  ) +
  theme_minimal() +
  scale_color_manual(
    breaks = site_order,
    values = c(
      milk  = "#16058b",
      nasal = "#8ab07c",
      oral  = "#e66d50",
      skin  = "#4DB6AC",
      gut   = "#a30543"
    )
  ) +
  theme(
    axis.ticks = element_line(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1),
    axis.line = element_blank()
  )

ggsave(
  file.path(out_dir, "p2_coverm_all_11958sample_Scatter_plot.png"),
  plot = p2_coverm_all_11958sample_Scatter_plot,
  width = 10,
  height = 8,
  units = "in"
)






#bar_plot
library(tidyverse)

stats <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/coverM/all_11958sample/all_11958sample_stats_with_U18AMPdb_name.csv",
  check.names = FALSE
)

out_dir <- "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/coverM/all_11958sample"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

top_n <- 50

plot_data <- stats %>%
  arrange(desc(prevalence)) %>%
  slice_head(n = top_n) %>%
  mutate(
    Contig = factor(Contig, levels = rev(Contig))
  )

scale_factor <- 6  

p1_coverm_all_11958sample_bar <- ggplot(
  plot_data,
  aes(x = Contig)
) +
  geom_col(
    aes(y = log10_mean_abundance),
    fill = "#274753",
    alpha = 0.8
  ) +
  geom_line(
    aes(y = prevalence * scale_factor),  
    group = 1,
    color = "lightblue",
    linewidth = 0.8
  ) +
  geom_point(
    aes(y = prevalence * scale_factor),
    color = "lightblue",
    size = 2
  ) +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, by = 2),
    name = "Mean abundance (log10(RPKM + 1), bar)",
    sec.axis = sec_axis(
      ~ . / scale_factor,                              
      name = "Prevalence (line)",
      breaks = seq(0, 1, by = 0.2)
    )
  ) +
  labs(
    x = "AMP cluster",
    title = paste("Top", top_n, "high-prevalence AMP clusters in all samples (N = 11958)")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  coord_flip()


ggsave(
  file.path(out_dir, "p1_coverm_all_11958sample_bar.pdf"),
  plot = p1_coverm_all_11958sample_bar,
  width = 10,
  height = 8,
  units = "in"
)






# ============================================================
# specialist and generalist AMP clusters for Fig.2D-E;Fig.S5
# ============================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

df <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/linux_data/632751_cluster_prev_meanabudance.csv")
sites <- c(
  "prev_gut",
  "prev_milk",
  "prev_nasal",
  "prev_oral",
  "prev_skin"
) 

df_sel <- df %>%filter(df$n_sites == 5)
thresholds <- seq(0, 0.5, by = 0.05)
results <- map_dfr(thresholds, function(thresh) {
  sub_df <- df_sel[, sites]
  min_vals <- apply(sub_df, 1, min, na.rm = TRUE)
  is_match <- min_vals > thresh
  is_match[is.na(is_match)] <- FALSE 
  tibble(threshold = thresh, count = sum(is_match))
})

p_gen <- ggplot(results, aes(x = threshold, y = count)) +
  geom_line(color = "#2C3E50", linewidth = 1.5) + # 绿色代表通用菌
  geom_point(color = "#7BC0CD", size = 4) +
  geom_text(aes(label = scales::comma(count)), vjust = -1.5, hjust = 0.5, size = 4, fontface = "bold") +
  scale_x_continuous(breaks = thresholds, labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title = "Sensitivity Analysis: Threshold vs.Generalist Count",
    x = "Prevalence Threshold (Minimum for each site)",
    y = "Number of Generalists (n_sites = 5)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA)
  )

print(p_gen)
write.csv(results,"C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/generalist_阈值.csv")



df <- df %>%filter(df$n_sites == 5)

fixed_thresh_gen <- 0.1

sub_df <- df[, sites]
min_vals <- apply(sub_df, 1, min, na.rm = TRUE)
cond_all_high <- min_vals >= fixed_thresh_gen
keep_rows_gen <- cond_all_high
df_filter_gen <- df[keep_rows_gen, ]
df_filter_gen <- if (nrow(df_filter_gen) > 0) {
  mean_vals <- rowMeans(df_filter_gen[, sites], na.rm = TRUE)
  df_final_gen <- df_filter_gen %>%
    mutate(
      class = "generalist",            # 定义类别
      mean_prevalence = mean_vals      # 记录平均值 (用于排序参考)
    ) %>%
    arrange(desc(mean_prevalence))     # 按平均丰度降序排列
}

# Top 40 选取40
Top_40_generalist <-df_filter_gen %>%
  slice_head(n = 40) %>%
  mutate(U18AMPdb_name = factor(U18AMPdb_name, levels = rev(U18AMPdb_name))) # 逆序排列，让最大值在最上面

write.csv(Top_40_generalist,"C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/Top_40_generalist.csv")
write.csv(df_filter_gen,"C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/df_filter_gen.csv")


##for tree
library(phangorn)  
library(ape) 
seqs <- read.FASTA("C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/generalist_top40_amp_40_aligned.fasta", type = "AA") 
cat(paste("成功读取", length(seqs), "条序列\n"))
cat(paste("序列长度示例:", nchar(seqs[[1]]), "\n")) 
seqs_phy <- phyDat(seqs, type = "AA")
class(seqs_phy)
dm  <- dist.ml(seqs_phy)
treeUPGMA  <- upgma(dm)
plot(treeUPGMA, main="UPGMA")  
#treeNJ  <- NJ(dm) ##UPGMA
#plot(treeNJ, "unrooted", main="NJ")
tree_phylo <- as.phylo(treeUPGMA)


plot(treeUPGMA)
p <- get("last_plot.phylo", envir = .PlotPhyloEnv)
amp_order <- treeUPGMA$tip.label[order(p$yy)]
amp_order
write.table(
  amp_order,
  "C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/AMP_plot_order.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)



# Niche of Representative
df_top40_niche <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/Top_40_generalist.csv")
amp_order <- read.table(
  "C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/AMP_plot_order.txt",
  stringsAsFactors = FALSE
)[,1]

cluster_summary <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Supplementary_data/Supplementary Data 4_632751_cluster_summary.csv")
cluster_summary_sel <-cluster_summary%>% 
  select(Cluster_name,Cluster_body_site,U18AMPdb_name)

df_top40_niche <- df_top40_niche%>% 
  left_join(cluster_summary_sel,
            by = "U18AMPdb_name")

df_top40_niche_order <- df_top40_niche %>%
  mutate(U18AMPdb_name = factor(U18AMPdb_name, levels = amp_order)) %>%  
  arrange(U18AMPdb_name)                                        

niche_colors <- c(
  "gut"   = "#a30543",  
  "skin"  = "#4DB6AC",  
  "oral"  = "#e66d50",  
  "nasal" = "#8ab07c",  
  "milk"  = "#16058b"   
)

p4 <- ggplot(df_top40_niche_order, aes(x = 1, y = U18AMPdb_name, fill = Cluster_body_site,alpha=0.5)) +
  geom_tile(color = "white", width = 0.4, height = 0.8) + 
  scale_fill_manual(values = niche_colors, name = "Niche") + 
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),      
    axis.ticks.x = element_blank(),     
    axis.title.x = element_blank(),     
    axis.title.y = element_blank(),     
    panel.grid = element_blank(),       
    legend.position = "right",          
    plot.margin = margin(5, 10, 5, 5) 
  ) +
  labs(title = "AMP Niche Distribution")

write.csv(df_top40_niche_order,"C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/df_top40_niche_order.csv")

print(p4)



# proportion of niche
df_proportion_top_40 <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/df_top40_niche_order.csv")
##sample_n
df_proportion_40_true <-df_proportion_top_40%>%
  select(U18AMPdb_name, starts_with("prop_")) %>%
  pivot_longer(
    cols = -U18AMPdb_name,          
    names_to = "site_code",  
    values_to = "prop_samples" 
  ) %>%
  mutate(
    site = case_when(
      site_code == "prop_gut" ~ "Gut",
      site_code == "prop_breastmilk" ~ "Breastmilk",
      site_code == "prop_nasal" ~ "Nasal",
      site_code == "prop_oral" ~ "Oral",
      site_code == "prop_skin" ~ "Skin"
    ),
    site = factor(site, levels = rev(c("Breastmilk", "Nasal", "Oral", "Skin","Gut"))),
    U18AMPdb_name = factor(U18AMPdb_name, levels = amp_order)
  )

my_colors <- c(
  "Gut"   = "#a30543",  
  "Skin"  = "#4DB6AC", 
  "Oral"  = "#e66d50", 
  "Nasal" = "#8ab07c",  
  "Breastmilk"  = "#16058b"   
)


p1 <- ggplot(df_proportion_40_true, aes(x = prop_samples, y =U18AMPdb_name , fill = site)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, alpha = 0.6) +
  scale_fill_manual(values = my_colors) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.ticks.y =element_blank(),
    axis.ticks.x =element_line(color = "black"),
    axis.title.y = element_blank(),
    axis.text.y =  element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(hjust = 0.5),
    legend.box.just = "center"
  )

print(p1)
write.csv(df_proportion_40_true,"C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/df_proportion_40_true.csv",row.names = FALSE)



# Mean_abundance
df_abund_top40<- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/df_top40_niche_order.csv")


df_abund_top40 <-df_abund_top40 %>%
  select(U18AMPdb_name, starts_with("log10_mean_abundance_")) %>%
  pivot_longer(
    cols = -U18AMPdb_name,        
    names_to = "site_code", 
    values_to = "log10_mean_abundance" 
  ) %>%
  mutate(
    site = case_when(
      site_code == "log10_mean_abundance_gut" ~ "Gut",
      site_code == "log10_mean_abundance_milk" ~ "Breastmilk",
      site_code == "log10_mean_abundance_nasal" ~ "Nasal",
      site_code == "log10_mean_abundance_oral" ~ "Oral",
      site_code == "log10_mean_abundance_skin" ~ "Skin"
    ),
    site = factor(site, levels = c("Breastmilk", "Nasal", "Oral", "Skin","Gut")),
    U18AMPdb_name = factor(U18AMPdb_name, levels = amp_order)
  )

p2 <- ggplot(df_abund_top40, aes(x = site, y = U18AMPdb_name,
                                 fill = log10_mean_abundance)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradientn(
    colors = c("#2166ACFF","#4393c3ff","#92c5deff","#D1e5f0ff","#FDDBC7FF","#F4a582ff","#B83945","#D73027"),
    name = "Mean Abundance(log10(RPKM+1))",
    limits = c(0, 5),      
    breaks = seq(0, 5, by = 1))+
  theme_void() +
  theme(axis.text.x = 
          element_text(size = 10, angle = 45, hjust =1,vjust=1),
        axis.text.y = 
          element_text(size = 10),
        legend.title.position = "top",
        legend.title = element_text(vjust=0.5,hjust=0.5,color="black"),
        legend.position = "right")
print(p2)
write.csv(df_abund_top40,"C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/df_abund_top40.csv",row.names = FALSE)




# prevalance_generalist
df_preval_top40 <- read.csv("C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/df_top40_niche_order.csv")

my_colors <- c(
  "Gut"   = "#a30543",  
  "Skin"  = "#4DB6AC", 
  "Oral"  = "#e66d50", 
  "Nasal" = "#8ab07c",  
  "Breastmilk"  = "#16058b"  
)


##prevalance
df_preval_top40 <- df_preval_top40 %>%
  select(U18AMPdb_name, starts_with("prev_")) %>%
  pivot_longer(
    cols = -U18AMPdb_name,          
    names_to = "site_code",  
    values_to = "prevalence" 
  ) %>%
  mutate(
    site = case_when(
      site_code == "prev_gut" ~ "Gut",
      site_code == "prev_milk" ~ "Breastmilk",
      site_code == "prev_nasal" ~ "Nasal",
      site_code == "prev_oral" ~ "Oral",
      site_code == "prev_skin" ~ "Skin"
    ),
    site = factor(site, levels = c("Breastmilk", "Nasal", "Oral", "Skin","Gut")),
    U18AMPdb_name = factor(U18AMPdb_name, levels = amp_order)
  )
# 绘图
p3 <- ggplot(df_preval_top40, aes(x = prevalence*100, y = U18AMPdb_name, group = site)) +
  geom_line(aes(color = site), linewidth = 0.6, alpha = 0.4, orientation = "y") +
  geom_point(aes(color = site), size = 1,alpha = 0.4) +
  facet_wrap(~ site, ncol = 5, strip.position = "top") +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(), 
    axis.line = element_line(color = "black"), 
    axis.ticks.x = element_line(color = "black"),
    axis.text.y =  element_text(size = 10), 
    axis.title.y = element_blank(),
    strip.background = element_rect(color = "black", fill = NA),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none",
  ) +
  labs(x = "Prevalence (100%)", y = "")
print(p3)
write.csv(df_abund_top40,"C:/Users/taoch/Desktop/632751_clusters/Fig_3DE/R_code/generalist/df_preval_top40.csv",row.names = FALSE)

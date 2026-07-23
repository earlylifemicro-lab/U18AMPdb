# ============================================================
#biological functions: for Fig.1B 
# ============================================================

library(tidyverse)
library(scales)
library(grid)

in_file <- "C:/Users/taoch/Desktop/632751_clusters/Fig_2_biological/632751_cluster_summary_with_iampcn_function_merged.csv"

out_dir <- "C:/Users/taoch/Desktop/632751_clusters/Fig_2_biological/"


df <- read_csv(
  in_file,
  show_col_types = FALSE,
  na = c("", "NA", "NaN", "null", "NULL")
)

cat("the number of lines of the file：", nrow(df), "\n")
cat("U18AMPdb：", n_distinct(df$U18AMPdb_name), "\n")


required_cols <- c(
  "U18AMPdb_name",
  "antibacterial",
  "anti-Gram-positive",
  "anti-Gram-negative",
  "antifungal",
  "antiviral",
  "anti-HIV",
  "antibiofilm",
  "anticancer",
  "anti-MRSA",
  "antiparasitic",
  "chemotactic",
  "anti-TB",
  "endotoxin",
  "insecticidal",
  "antimalarial",
  "anticandida",
  "antiplasmodial",
  "antiprotozoal"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "The input file is missing as follows：\n",
      paste(missing_cols, collapse = ", ")
    )
  )
}


to_logical_yes <- function(x) {
  
  x <- str_to_lower(
    str_trim(
      as.character(x)
    )
  )
  
  x %in% c(
    "yes",
    "true",
    "1",
    "positive",
    "y"
  )
}

function_source_cols <- setdiff(
  required_cols,
  "U18AMPdb_name"
)

df_binary <- df %>%
  mutate(
    across(
      all_of(function_source_cols),
      to_logical_yes
    )
  )


apd_function_wide <- tibble(
  
  U18AMPdb_name = df_binary$U18AMPdb_name,
  
  Antibacterial =
    df_binary$antibacterial |
    df_binary$`anti-Gram-positive` |
    df_binary$`anti-Gram-negative`,
  
  Antibiofilm =
    df_binary$antibiofilm,
  
  `Anti-MRSA` =
    df_binary$`anti-MRSA`,
  
  `Anti-TB` =
    df_binary$`anti-TB`,
  
  `Anti-endotoxin` =
    df_binary$endotoxin,
  
  Antiviral =
    df_binary$antiviral,
  
  `Anti-HIV` =
    df_binary$`anti-HIV`,
  
  Antifungal =
    df_binary$antifungal,
  
  `Anti-Candida` =
    df_binary$anticandida,
  
  Antiparasitic =
    df_binary$antiparasitic |
    df_binary$antiprotozoal,
  
  # antiplasmodial并入Antimalarial
  Antimalarial =
    df_binary$antimalarial |
    df_binary$antiplasmodial,
  
  Anticancer =
    df_binary$anticancer,
  
  Chemotactic =
    df_binary$chemotactic,
  
  Insecticidal =
    df_binary$insecticidal
)


apd_function_long <- apd_function_wide %>%
  pivot_longer(
    cols = -U18AMPdb_name,
    names_to = "Function",
    values_to = "Positive"
  )


total_amp <- 632751L


observed_amp <- apd_function_wide %>%
  filter(
    !is.na(U18AMPdb_name),
    U18AMPdb_name != ""
  ) %>%
  summarise(
    n = n_distinct(U18AMPdb_name)
  ) %>%
  pull(n)

missing_from_function_file <- total_amp - observed_amp



function_stats <- apd_function_long %>%
  filter(
    Positive,
    !is.na(U18AMPdb_name),
    U18AMPdb_name != ""
  ) %>%
  distinct(
    U18AMPdb_name,
    Function
  ) %>%
  count(
    Function,
    name = "AMP_count"
  ) %>%
  arrange(
    desc(AMP_count)
  ) %>%
  mutate(
    Annotation_percentage =
      AMP_count / sum(AMP_count) * 100,
    
    Cluster_percentage =
      AMP_count / total_amp * 100
  )

if (nrow(function_stats) == 0) {
  stop("no function annote")
}

write_csv(
  function_stats,
  file.path(
    out_dir,
    "APD_biological_function_statistics.csv"
  )
)

print(function_stats)

cat(
  "阳性功能注释总数：",
  sum(function_stats$AMP_count),
  "\n"
)


annotated_amp <- apd_function_long %>%
  filter(
    !is.na(U18AMPdb_name),
    U18AMPdb_name != ""
  ) %>%
  group_by(U18AMPdb_name) %>%
  summarise(
    Any_APD_function = any(
      Positive,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

n_annotated <- sum(
  annotated_amp$Any_APD_function,
  na.rm = TRUE
)


n_unannotated <- total_amp - n_annotated

summary_df <- tibble(
  Item = c(
    "Total U18AMPdb clusters",
    "Clusters present in function annotation file",
    "Clusters absent from function annotation file",
    "Clusters with at least one APD function",
    "Clusters without APD-listed functions",
    "Total positive function annotations"
  ),
  Number = c(
    total_amp,
    observed_amp,
    missing_from_function_file,
    n_annotated,
    n_unannotated,
    sum(function_stats$AMP_count)
  )
)

write_csv(
  summary_df,
  file.path(
    out_dir,
    "APD_biological_function_summary.csv"
  )
)

print(summary_df)


function_order <- c(
  "Antibacterial",
  "Antibiofilm",
  "Anti-MRSA",
  "Anti-TB",
  "Anti-endotoxin",
  "Antiviral",
  "Anti-HIV",
  "Antifungal",
  "Anti-Candida",
  "Antiparasitic",
  "Antimalarial",
  "Anticancer",
  "Chemotactic",
  "Insecticidal"
)

function_stats_plot <- function_stats %>%
  mutate(
    Function = factor(
      Function,
      levels = function_order
    )
  ) %>%
  arrange(Function) %>%
  mutate(
    Legend_label = paste0(
      as.character(Function),
      "  ",
      comma(AMP_count),
      " (",
      sprintf("%.1f", Annotation_percentage),
      "%)"
    ),
    
    Internal_label = if_else(
      Annotation_percentage >= 3,
      paste0(
        comma(AMP_count),
        "\n",
        sprintf("%.1f%%", Annotation_percentage)
      ),
      ""
    )
  )

function_stats_plot <- function_stats_plot %>%
  filter(!is.na(Function))



function_colors <- c(
  
  "Antibacterial"   = "#B54A34",
  "Antibiofilm"     = "#4F9275",
  "Anti-MRSA"       = "#70A487",
  "Anti-TB"         = "#9DBBA3",
  "Anti-endotoxin"  = "#C6D4C2",
  
  "Antiviral"       = "#82A9CF",
  "Anti-HIV"        = "#B4CCE3",
  
  "Antifungal"      = "#78999A",
  "Anti-Candida"    = "#AFC5C2",
  

  "Antiparasitic"   = "#899DBC",
  "Antimalarial"    = "#AAB6D2",
  

  "Anticancer"      = "#C47764",
  "Chemotactic"     = "#8DA4B3",
  "Insecticidal"    = "#C9CECC"
)


function_colors_used <- function_colors[
  as.character(function_stats_plot$Function)
]

legend_labels <- function_stats_plot$Legend_label

names(legend_labels) <- as.character(
  function_stats_plot$Function
)


p_apd_donut <- ggplot(
  function_stats_plot,
  aes(
    x = 2,
    y = AMP_count,
    fill = Function
  )
) +
  
  geom_col(
    width = 0.88,
    color = "white",
    linewidth = 0.65
  ) +
  
  geom_text(
    aes(label = Internal_label),
    position = position_stack(vjust = 0.5),
    size = 3.1,
    lineheight = 0.88,
    color = "black"
  ) +
  
  coord_polar(
    theta = "y",
    start = 0,
    direction = 1
  ) +
  
  xlim(0.45, 2.65) +
  
  scale_fill_manual(
    values = function_colors_used,
    breaks = names(legend_labels),
    labels = legend_labels,
    drop = FALSE,
    name = "Biological function"
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 0,
    label = paste0(
      "APD functions\n",
      comma(n_annotated),
      " AMPs"
    ),
    size = 4.2,
    fontface = "bold",
    lineheight = 0.95,
    color = "#303030"
  ) +
  
  labs(
    title = "Biological functions of U18AMPdb peptides",
    subtitle = paste0(
      "Percentage of all ",
      comma(total_amp),
      " U18AMPdb clusters assigned to each function"
    ),
    x = "Percentage of all U18AMPdb clusters",
    y = NULL,
    caption = paste0(
      "Each function was calculated independently. ",
      "Functional categories are not mutually exclusive; ",
      "therefore, percentages do not sum to 100%."
    )
  )+
  
  theme_void(base_size = 12) +
  
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      color = "#202020",
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      size = 10.5,
      hjust = 0.5,
      color = "#555555",
      margin = margin(b = 10)
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 10.5,
      face = "bold",
      color = "#202020"
    ),
    
    legend.text = element_text(
      size = 8.8,
      color = "#303030"
    ),
    
    legend.key.size = unit(
      0.52,
      "cm"
    ),
    
    legend.spacing.y = unit(
      0.08,
      "cm"
    ),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )

print(p_apd_donut)



##for Fig.S2A
# ============================================================
# Cluster_percentage bar plot
# ============================================================

library(tidyverse)
library(scales)


out_dir <- "C:/Users/taoch/Desktop/Fig_2_biological"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


function_colors <- c(
  "Antibacterial"  = "#B54A34",
  "Antibiofilm"    = "#4F9275",
  "Anti-MRSA"      = "#70A487",
  "Anti-TB"        = "#9DBBA3",
  "Anti-endotoxin" = "#C6D4C2",
  "Antiviral"      = "#82A9CF",
  "Anti-HIV"       = "#B4CCE3",
  "Antifungal"     = "#78999A",
  "Anti-Candida"   = "#AFC5C2",
  "Antiparasitic"  = "#899DBC",
  "Antimalarial"   = "#AAB6D2",
  "Anticancer"     = "#C47764",
  "Chemotactic"    = "#8DA4B3",
  "Insecticidal"   = "#C9CECC"
)


function_stats_plot <- function_stats %>%
  filter(
    !is.na(Cluster_percentage),
    AMP_count > 0
  ) %>%
  arrange(Cluster_percentage) %>%
  mutate(
    Function = factor(
      Function,
      levels = Function
    ),

    Bar_label = paste0(
      comma(AMP_count),
      "  (",
      sprintf("%.2f", Cluster_percentage),
      "%)"
    )
  )


print(function_stats_plot)

x_max <- max(
  function_stats_plot$Cluster_percentage,
  na.rm = TRUE
)


x_limit <- x_max * 1.28


if (x_limit < 1) {
  x_limit <- 1
}


p_cluster_percentage <- ggplot(
  function_stats_plot,
  aes(
    x = Cluster_percentage,
    y = Function,
    fill = Function
  )
) +
  
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.4
  ) +
  

  geom_text(
    aes(label = Bar_label),
    hjust = -0.08,
    size = 3.4,
    color = "#303030"
  ) +
  
  scale_fill_manual(
    values = function_colors,
    guide = "none"
  ) +
  
  scale_x_continuous(
    labels = label_number(
      accuracy = 0.1,
      suffix = "%"
    ),
    limits = c(0, x_limit),
    expand = expansion(mult = c(0, 0))
  ) +
  
  labs(
    title = "Biological functions of U18AMPdb peptides",
    subtitle = paste0(
      "Percentage of all ",
      comma(total_amp),
      " AMP clusters assigned to each function"
    ),
    x = "Percentage of all AMP clusters",
    y = NULL,
    caption = paste0(
      "Each function was calculated independently. ",
      "Because one AMP can have multiple functions, ",
      "the percentages do not sum to 100%."
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 10.5,
      color = "#555555",
      margin = margin(
        t = 3,
        b = 12
      )
    ),
    
    axis.title.x = element_text(
      size = 11,
      face = "bold",
      margin = margin(t = 8)
    ),
    
    axis.text.x = element_text(
      size = 9.5,
      color = "#333333"
    ),
    
    axis.text.y = element_text(
      size = 10,
      color = "#222222"
    ),
    
    axis.line.y = element_blank(),
    
    axis.ticks.y = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "#E4E4E4",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.caption = element_text(
      size = 8.5,
      color = "#666666",
      hjust = 0,
      margin = margin(t = 12)
    ),
    
    plot.margin = margin(
      t = 12,
      r = 35,
      b = 10,
      l = 10
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  )

print(p_cluster_percentage)







# ============================================================
#Amino acid frequency for Fig.1D,1E
# ============================================================

AMP_public_database_AA <- read.csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_2_homology/property/AA/new_AMP_public_database_AA_frequency.csv")

AMPSphere_AA <- read.csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_2_homology/property/AA/AMPSphere_AA_frequency.csv")

AEPs_AA <- read.csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_2_homology/property/AA/AEPs_AA_frequency.csv")

U18AMPdb_AA <- read.csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_2_homology/property/AA/U18AMPdb_AA_frequency.csv")


combined_data_AA <- bind_rows(AMP_public_database_AA, AMPSphere_AA, AEPs_AA, U18AMPdb_AA)
write.csv(combined_data_AA,"C:/Users/taoch/Desktop/632751_clusters_文章/Fig_2_homology/R_code/for_AMPdatabase/AA_for_AMPdatabase.csv",row.names = FALSE)


Freq <- c("freq_A", "freq_C", "freq_D", "freq_E", "freq_F", "freq_G", "freq_H", "freq_I", "freq_K", "freq_L", "freq_M", "freq_N", "freq_P", "freq_Q", "freq_R", "freq_S", "freq_T", "freq_V", "freq_W", "freq_Y")

df_long_AA <- combined_data_AA %>%
  pivot_longer(
    cols = all_of(Freq),
    names_to = "Frequence",
    values_to = "Value"
  )

colors_AA <- c(
  "U18AMPdb" = "#AE1E24",
  
  "AMPSphere" = "#0F976F",  
  "AMP public datasets"= "#BDCCB5",
  "AEPs" = "#5f8a8b"
)

desired_order <- names(colors_AA)
df_long_AA <- df_long_AA %>%
  mutate(Group = factor(Group, levels = desired_order))

p3_freq <- ggplot(df_long_AA, aes(x = Frequence, y = Value, fill = Group)) +
  geom_col(position = "dodge",alpha = 0.8) +
  labs(title = "Amino Acid Frequency by Database", 
       x = "Amino Acid", 
       y = "Frequency") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,hjust = 1,size = 11),
    axis.text.y = element_text(size = 11),
    axis.ticks = element_line(),
    legend.position = "top",
    legend.title = element_text(size = 8),
    panel.border = element_rect(colour = "black",linewidth = 0.6,fill=NA),
    panel.grid = element_blank())+
  scale_fill_manual(values = colors_AA )

print(p3_freq)








# ============================================================
# physicochemical properties for Fig.1F,1G
# ============================================================

library(tidyverse)
library(ggplot2)
library(ggridges)
library(patchwork)
library(scales)


prop_dir <- "C:/Users/taoch/Desktop/632751_clusters_文章/Fig_2_homology/property/property"

read_prop <- function(file, db_name) {
  dat <- read.csv(
    file.path(prop_dir, file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  dat$Database <- db_name
  return(dat)
}

AMP_public_database <- read_prop(
  "new_AMP_public_database_final_all_properties.csv",
  "AMP public datasets"
)

AMPSphere <- read_prop(
  "AMPSphere_final_all_properties.csv",
  "AMPSphere"
)

AEPs <- read_prop(
  "AEPs_final_all_properties.csv",
  "AEPs"
)

U18AMPdb <- read_prop(
  "U18AMPdb_final_all_properties.csv",
  "U18AMPdb"
)


combined_data <- bind_rows(
  U18AMPdb,
  AMPSphere,
  AMP_public_database,
  AEPs
)

write.csv(combined_data,"C:/Users/taoch/Desktop/632751_clusters/Fig_2_homology/R_code/for_AMPdatabase/combined_data_properity_for_AMPdatabase.csv",row.names = FALSE)


plot_properties <- c(
  "Length",
  "Charge",
  "pI",
  "BomanInd",
  "Hydrophobicity_H",
  "HydrophobicMoment_uH"
)


missing_cols <- setdiff(c("Database", plot_properties), colnames(combined_data))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "The following does not exist in combined_data. Please check the column names：",
      paste(missing_cols, collapse = ", ")
    )
  )
}



df_long <- combined_data %>%
  select(all_of(c("Database", plot_properties))) %>%
  pivot_longer(
    cols = all_of(plot_properties),
    names_to = "Property",
    values_to = "Value"
  ) %>%
  mutate(
    Value = suppressWarnings(as.numeric(Value)),
    Property = factor(Property, levels = plot_properties)
  ) %>%
  filter(!is.na(Value), is.finite(Value))


property_labels <- c(
  "Length" = "Length",
  "MW" = "MW",
  "Charge" = "Charge",
  "ChargeDensity" = "Charge density",
  "pI" = "pI",
  "InstabilityInd" = "Instability index",
  "Aromaticity" = "Aromaticity",
  "AliphaticInd" = "Aliphatic index",
  "BomanInd" = "Boman index",
  "HydrophRatio" = "Hydrophobic ratio",
  "Hydrophobicity_H" = "Hydrophobicity",
  "HydrophobicMoment_uH" = "Hydrophobic moment"
)


database_colors <- c(
  "U18AMPdb" = "#AE1E24",
  
  # AMP 库
  "AMPSphere" = "#0F976F",
  "AMP public datasets" = "#BDCCB5",
  "AEPs" = "#5f8a8b"
)


make_plot_df <- function(dat, group_cols = "Database", max_n = 50000,
                         q_low = 0.01, q_high = 0.99) {
  
  required_cols <- c(group_cols, "Property", "Value")
  missing_cols <- setdiff(required_cols, colnames(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following columns are missing: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  dat2 <- dat %>%
    filter(!is.na(Value), is.finite(Value)) %>%
    group_by(Property) %>%
    mutate(
      xmin_plot = quantile(Value, q_low, na.rm = TRUE),
      xmax_plot = quantile(Value, q_high, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(Value >= xmin_plot, Value <= xmax_plot) %>%
    select(-xmin_plot, -xmax_plot)
  
  dat2 <- dat2 %>%
    group_by(across(all_of(c(group_cols, "Property")))) %>%
    group_modify(~ {
      if (nrow(.x) > max_n) {
        dplyr::slice_sample(.x, n = max_n)
      } else {
        .x
      }
    }) %>%
    ungroup()
  
  return(dat2)
}

compare_dbs <- setdiff(unique(df_long$Database), "U18AMPdb")

wilcox_results <- map_dfr(compare_dbs, function(db) {
  
  map_dfr(plot_properties, function(prop) {
    
    x <- df_long %>%
      filter(Database == "U18AMPdb", Property == prop) %>%
      pull(Value)
    
    y <- df_long %>%
      filter(Database == db, Property == prop) %>%
      pull(Value)
    
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    
    if (length(x) < 2 | length(y) < 2) {
      return(tibble(
        Property = prop,
        Comparison = paste0("U18AMPdb vs ", db),
        Database_ref = "U18AMPdb",
        Database_compared = db,
        n_U18AMPdb = length(x),
        n_compared = length(y),
        median_U18AMPdb = median(x, na.rm = TRUE),
        median_compared = median(y, na.rm = TRUE),
        delta_median = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
        p_value = NA_real_
      ))
    }
    
    test <- wilcox.test(
      x,
      y,
      alternative = "two.sided",
      exact = FALSE
    )
    
    tibble(
      Property = prop,
      Comparison = paste0("U18AMPdb vs ", db),
      Database_ref = "U18AMPdb",
      Database_compared = db,
      n_U18AMPdb = length(x),
      n_compared = length(y),
      median_U18AMPdb = median(x, na.rm = TRUE),
      median_compared = median(y, na.rm = TRUE),
      delta_median = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
      p_value = test$p.value
    )
  })
}) %>%
  group_by(Property) %>%  ##```r

# Wilcoxon test: U18AMPdb vs each external AMP databasewith BH correction
compare_dbs <- setdiff(unique(df_long$Database), "U18AMPdb")


wilcox_results <- map_dfr(compare_dbs, function(db) {
  
  map_dfr(plot_properties, function(prop) {
    
    x <- df_long %>%
      filter(Database == "U18AMPdb", Property == prop) %>%
      pull(Value)
    
    y <- df_long %>%
      filter(Database == db, Property == prop) %>%
      pull(Value)
    
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    
    if (length(x) < 2 | length(y) < 2) {
      return(tibble(
        Property = prop,
        Comparison = paste0("U18AMPdb vs ", db),
        Database_ref = "U18AMPdb",
        Database_compared = db,
        n_U18AMPdb = length(x),
        n_compared = length(y),
        median_U18AMPdb = median(x, na.rm = TRUE),
        median_compared = median(y, na.rm = TRUE),
        delta_median = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
        p_value = NA_real_
      ))
    }
    
    test <- wilcox.test(
      x,
      y,
      alternative = "two.sided",
      exact = FALSE
    )
    
    tibble(
      Property = prop,
      Comparison = paste0("U18AMPdb vs ", db),
      Database_ref = "U18AMPdb",
      Database_compared = db,
      n_U18AMPdb = length(x),
      n_compared = length(y),
      median_U18AMPdb = median(x, na.rm = TRUE),
      median_compared = median(y, na.rm = TRUE),
      delta_median = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
      p_value = test$p.value
    )
  })
}) %>%
  group_by(Property) %>%  ##每个Property内部对 3 个比较做 BH 校正
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    direction = case_when(
      delta_median > 0 ~ "Higher in U18AMPdb",
      delta_median < 0 ~ "Lower in U18AMPdb",
      TRUE ~ "No median difference"
    )
  ) %>%
  ungroup()


print(wilcox_results)

  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance = case_when(
      is.na(p_adj_BH) ~ "NA",
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    direction = case_when(
      delta_median > 0 ~ "Higher in U18AMPdb",
      delta_median < 0 ~ "Lower in U18AMPdb",
      TRUE ~ "No median difference"
    )
  ) %>%
  ungroup()


print(wilcox_results)


write.csv(
  wilcox_results,
  file = "C:/Users/taoch/Desktop/632751_clusters_文章/Fig_2_homology/R_code/for_AMPdatabase/U18AMPdb_vs_AMP_databases_wilcox_BH.csv",
  row.names = FALSE
)



u18_median <- df_long %>%
  filter(Database == "U18AMPdb") %>%
  group_by(Property) %>%
  summarise(
    U18_median = median(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Property = factor(Property, levels = plot_properties)
  )

if (nrow(u18_median) == 0) {
  stop("没有找到 U18AMPdb 数据，请检查 Database 名称是否为 U18AMPdb。")
}


db_order_top_to_bottom <- c(
  "U18AMPdb",
  "AMPSphere",
  "AMP public datasets",
  "AEPs"
)

db_order_present <- db_order_top_to_bottom[
  db_order_top_to_bottom %in% unique(as.character(df_long$Database))
]

df_detail <- df_long %>%
  filter(Database %in% db_order_present)

df_detail_plot <- make_plot_df(
  df_detail,
  group_cols = "Database",
  max_n = 50000
)

df_detail_plot$Database <- factor(
  df_detail_plot$Database,
  levels = rev(db_order_present)
)

df_detail_plot <- df_detail_plot %>%
  mutate(
    Outline = if_else(as.character(Database) == "U18AMPdb", "U18AMPdb", "Others"),
    Ridge_alpha = if_else(as.character(Database) == "U18AMPdb", "U18AMPdb", "Others")
  )

database_label <- c(
  "U18AMPdb" = "18AMPdb",
  "AMPSphere" = "AMPSphere",
  "AEPs" = "AEPs",
  "AMP public datasets" = "AMP public datasets"
)

P_detail <- ggplot(
  df_detail_plot,
  aes(x = Value, y = Database, fill = Database)
) +
  geom_density_ridges(
    aes(color = Outline, alpha = Ridge_alpha),
    scale = 1.05,
    linewidth = 0.28,
    quantile_lines = TRUE,
    quantiles = 2,
    rel_min_height = 0.01
  ) +
  geom_vline(
    data = u18_median,
    aes(xintercept = U18_median),
    inherit.aes = FALSE,
    color = "#AE1E24",
    linetype = "dashed",
    linewidth = 0.45,
    alpha = 0.9
  ) +
  scale_fill_manual(
    values = database_colors,
    breaks = db_order_present,
    labels = database_label[db_order_present]
  ) +
  scale_color_manual(
    values = c(
      "U18AMPdb" = "black",
      "Others" = "white"
    ),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c(
      "U18AMPdb" = 0.95,
      "Others" = 0.55
    ),
    guide = "none"
  ) +
  scale_y_discrete(labels = function(x) database_label[x]) +
  facet_wrap(
    ~Property,
    scales = "free_x",
    ncol = 3,
    labeller = labeller(Property = property_labels)
  ) +
  labs(
    x = "Value",
    y = NULL,
    fill = "Database"
  ) +
  theme_ridges() +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey85", color = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11)
  )

print(P_detail)


#plot_sub_properties
plot_sub_properties <- c(
  "MW",
  "ChargeDensity",
  "InstabilityInd",
  "Aromaticity",
  "AliphaticInd",
  "HydrophRatio"
)


missing_cols <- setdiff(c("Database", plot_sub_properties ), colnames(combined_data))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "The following does not exist in combined_data：",
      paste(missing_cols, collapse = ", ")
    )
  )
}


df_long_sub <- combined_data %>%
  select(all_of(c("Database", plot_sub_properties ))) %>%
  pivot_longer(
    cols = all_of(plot_sub_properties ),
    names_to = "Property",
    values_to = "Value"
  ) %>%
  mutate(
    Value = suppressWarnings(as.numeric(Value)),
    Property = factor(Property, levels = plot_sub_properties )
  ) %>%
  filter(!is.na(Value), is.finite(Value))



property_labels <- c(
  "Length" = "Length",
  "MW" = "MW",
  "Charge" = "Charge",
  "ChargeDensity" = "Charge density",
  "pI" = "pI",
  "InstabilityInd" = "Instability index",
  "Aromaticity" = "Aromaticity",
  "AliphaticInd" = "Aliphatic index",
  "BomanInd" = "Boman index",
  "HydrophRatio" = "Hydrophobic ratio",
  "Hydrophobicity_H" = "Hydrophobicity",
  "HydrophobicMoment_uH" = "Hydrophobic moment"
)


database_colors <- c(
  "U18AMPdb" = "#AE1E24",
  
  # AMP datasets
  "AMPSphere" = "#0F976F",
  "AMP public datasets" = "#BDCCB5",
  "AEPs" = "#5f8a8b"
)



make_plot_df <- function(dat, group_cols = "Database", max_n = 50000,
                         q_low = 0.01, q_high = 0.99) {
  
  required_cols <- c(group_cols, "Property", "Value")
  missing_cols <- setdiff(required_cols, colnames(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "The following columns are missing: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  

  dat2 <- dat %>%
    filter(!is.na(Value), is.finite(Value)) %>%
    group_by(Property) %>%
    mutate(
      xmin_plot = quantile(Value, q_low, na.rm = TRUE),
      xmax_plot = quantile(Value, q_high, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(Value >= xmin_plot, Value <= xmax_plot) %>%
    select(-xmin_plot, -xmax_plot)
  
  dat2 <- dat2 %>%
    group_by(across(all_of(c(group_cols, "Property")))) %>%
    group_modify(~ {
      if (nrow(.x) > max_n) {
        dplyr::slice_sample(.x, n = max_n)
      } else {
        .x
      }
    }) %>%
    ungroup()
  
  return(dat2)
}


u18_median <- df_long_sub %>%
  filter(Database == "U18AMPdb") %>%
  group_by(Property) %>%
  summarise(
    U18_median = median(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Property = factor(Property, levels = plot_sub_properties )
  )

if (nrow(u18_median) == 0) {
  stop("No U18AMPdb data was found")
}


db_order_top_to_bottom <- c(
  "U18AMPdb",
  "AMPSphere",
  "AMP public datasets",
  "AEPs"
)

db_order_present <- db_order_top_to_bottom[
  db_order_top_to_bottom %in% unique(as.character(df_long$Database))
]

df_detail <- df_long_sub %>%
  filter(Database %in% db_order_present)

df_detail_plot <- make_plot_df(
  df_detail,
  group_cols = "Database",
  max_n = 50000
)

df_detail_plot$Database <- factor(
  df_detail_plot$Database,
  levels = rev(db_order_present)
)

df_detail_plot <- df_detail_plot %>%
  mutate(
    Outline = if_else(as.character(Database) == "U18AMPdb", "U18AMPdb", "Others"),
    Ridge_alpha = if_else(as.character(Database) == "U18AMPdb", "U18AMPdb", "Others")
  )

database_label <- c(
  "U18AMPdb" = "18AMPdb",
  "AMPSphere" = "AMPSphere",
  "AEPs" = "AEPs",
  "AMP public datasets" = "AMP public datasets"
)

P_detail_sub <- ggplot(
  df_detail_plot,
  aes(x = Value, y = Database, fill = Database)
) +
  geom_density_ridges(
    aes(color = Outline, alpha = Ridge_alpha),
    scale = 1.05,
    linewidth = 0.28,
    quantile_lines = TRUE,
    quantiles = 2,
    rel_min_height = 0.01
  ) +
  geom_vline(
    data = u18_median,
    aes(xintercept = U18_median),
    inherit.aes = FALSE,
    color = "#AE1E24",
    linetype = "dashed",
    linewidth = 0.45,
    alpha = 0.9
  ) +
  scale_fill_manual(
    values = database_colors,
    breaks = db_order_present,
    labels = database_label[db_order_present]
  ) +
  scale_color_manual(
    values = c(
      "U18AMPdb" = "black",
      "Others" = "white"
    ),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c(
      "U18AMPdb" = 0.95,
      "Others" = 0.55
    ),
    guide = "none"
  ) +
  scale_y_discrete(labels = function(x) database_label[x]) +
  facet_wrap(
    ~Property,
    scales = "free_x",
    ncol = 3,
    labeller = labeller(Property = property_labels)
  ) +
  labs(
    x = "Value",
    y = NULL,
    fill = "Database"
  ) +
  theme_ridges() +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey85", color = NA),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_line(),
    panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11)
  )

print(P_detail_sub)

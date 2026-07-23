# ============================================================
#Hierarchical taxonomic flow of U18AMPdb for Fig.3A; Fig.S6
# ============================================================

lineage_wide_change_gut <- read.table(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4A_taxa_sankey/gut_taxa_sankey/figure_6_sankey/domain_to_species/lineage_wide_gut_gtdbtk.bac120_76133.tsv",
  header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE
)

library(dplyr)
library(tidyverse)
library(networkD3)
library(htmlwidgets)
data_gut <- lineage_wide_change_gut %>% select(-any_of("user_genome"))

data_counted <- data_gut %>%
  group_by(Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(count = n(), .groups = "drop")

threshold <- 8000   


Domain_total  <- data_counted %>% group_by(Domain)  %>% summarise(total=sum(count),.groups="drop")
phylum_total  <- data_counted %>% group_by(Phylum)  %>% summarise(total=sum(count),.groups="drop")
class_total   <- data_counted %>% group_by(Class)   %>% summarise(total=sum(count),.groups="drop")
order_total   <- data_counted %>% group_by(Order)   %>% summarise(total=sum(count),.groups="drop")
family_total  <- data_counted %>% group_by(Family)  %>% summarise(total=sum(count),.groups="drop")
genus_total   <- data_counted %>% group_by(Genus)   %>% summarise(total=sum(count),.groups="drop")
species_total <- data_counted %>% group_by(Species) %>% summarise(total = sum(count), .groups = "drop")

# Top & low 
top_domain    <- head(Domain_total$Domain,5)
low_domain    <- Domain_total$Domain[Domain_total$total < threshold]

top_phyla    <- head(phylum_total$Phylum,5)
low_phyla    <- phylum_total$Phylum[phylum_total$total < threshold]

top_classes  <- head(class_total$Class, 6)
low_classes  <- class_total$Class[class_total$total < threshold]

top_orders   <- head(order_total$Order, 8)
low_orders   <- order_total$Order[order_total$total < threshold]

top_families <- head(family_total$Family, 10)
low_families <- family_total$Family[family_total$total < threshold]

top_genus  <- head(genus_total$Genus, 10)
low_genus  <- genus_total$Genus[genus_total$total < threshold]

top_species  <- head(species_total$Species, 20)
low_species  <- species_total$Species[species_total$total < threshold]


data_counted_simplified <- data_counted %>%
  mutate(
    Phylum  = ifelse(!Phylum  %in% top_phyla    & Phylum  %in% low_phyla,    "Other_Phylum",  Phylum),
    Class   = ifelse(!Class   %in% top_classes  & Class   %in% low_classes,  "Other_Class",   Class),
    Order   = ifelse(!Order   %in% top_orders   & Order   %in% low_orders,   "Other_Order",   Order),
    Family  = ifelse(!Family  %in% top_families & Family  %in% low_families, "Other_Family",  Family),
    Genus   = ifelse(!Genus   %in% top_genus    & Genus   %in% low_genus,   "Other_Genus",   Genus),
    Species = ifelse(!Species %in% top_species  & Species %in% low_species, "Other_Species", Species)
  ) %>%
  group_by(Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    Class  = ifelse(Phylum == "Other_Phylum",  "Other_Class",  Class),
    Order  = ifelse(Class  == "Other_Class",   "Other_Order",  Order),
    Family = ifelse(Order  == "Other_Order",   "Other_Family", Family),
    Genus  = ifelse(Family == "Other_Family",  "Other_Genus",  Genus),
    Species= ifelse(Genus  == "Other_Genus",   "Other_Species",Species)
  ) %>%
  group_by(Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(count = sum(count), .groups = "drop")


links_char <- data_counted_simplified %>%
  transmute(source = Domain, target = Phylum, Value = count) %>%
  bind_rows(
    data_counted_simplified %>% transmute(source = Phylum,  target = Class,  Value = count),
    data_counted_simplified %>% transmute(source = Class,   target = Order,  Value = count),
    data_counted_simplified %>% transmute(source = Order,   target = Family, Value = count),
    data_counted_simplified %>% transmute(source = Family,  target = Genus,  Value = count),
    data_counted_simplified %>% transmute(source = Genus,   target = Species,Value = count)
  ) %>%
  group_by(source, target) %>%
  summarise(Value = sum(Value), .groups = "drop")

nodes <- data.frame(name = unique(c(links_char$source, links_char$target)),
                    stringsAsFactors = FALSE)

domain_order  <- c("Bacteria")

phylum_order <- c("Fusobacteriota","Pseudomonadota","Bacillota",
                  "Bacteroidota","Actinomycetota","Campylobacterota",
                  "Other_Phylum") 

class_order   <- c("Fusobacteriia","Gammaproteobacteria","Bacilli",
                   "Clostridia","Negativicutes","Campylobacteria",
                   "Bacteroidia","Actinomycetes","Other_Class")

order_order <- c("Fusobacteriales","Cardiobacteriales","Burkholderiales",
                 "Staphylococcales","Enterobacterales","Lactobacillales",
                 "Veillonellales","Flavobacteriales","Bacteroidales",
                 "Actinomycetales","Other_Order")

family_order <- c("Mycobacteriaceae","Actinomycetaceae","Staphylococcaceae",
                  "Lactobacillaceae","Pseudomonadaceae","Enterobacteriaceae",
                  "Other_Family")

genus_order <- c("Bifidobacterium", "Escherichia", "Staphylococcus", "Lactobacillus",
                 "Pseudomonas", "Bacteroides", "Fusobacterium", "Other_Genus")

species_order <- c("Bifidobacterium breve", "Escherichia coli", "Staphylococcus aureus",
                   "Lactobacillus rhamnosus", "Pseudomonas aeruginosa", "Bacteroides fragilis",
                   "Fusobacterium nucleatum", "Other_Species")


custom_order <- c(domain_order, phylum_order, class_order, order_order, family_order, genus_order, species_order)


nodes <- nodes %>%
  mutate(order_index = match(name, custom_order)) %>%
  arrange(order_index) %>%
  mutate(
    order_index = ifelse(is.na(order_index),
                         max(order_index, na.rm = TRUE) + seq_len(sum(is.na(order_index))),
                         order_index)
  ) %>%
  arrange(order_index) %>%
  mutate(id = row_number() - 1)


links <- links_char %>%
  mutate(
    source = match(source, nodes$name) - 1,
    target = match(target, nodes$name) - 1
  )

stopifnot(!any(is.na(links$source)), !any(is.na(links$target)))


links$group <- nodes$name[links$source + 1]


fixed_colors <- list(
"Bacteria" = "#457B9D",
"Other_Phylum"="#c1bebf",
"Other_Class"="#c1bebf",
"Other_Order"="#c1bebf",
"Other_Family"="#c1bebf",
"Other_Genus"="#c1bebf",
"Other_Species"="#c1bebf",

"Actinomycetota"="#06D6A0",
"Actinomycetes"="#43BEB5",
"Actinomycetales"="#5EADCA",
"Bifidobacteriaceae"="#7C9FBF",
"Bifidobacterium"="#9991B4",

"Pseudomonadota"="#687BA9",
"Gammaproteobacteria"="#B38287",
"Enterobacterales"="#927B8E",
"Enterobacteriaceae"="#7A7594",

"Bacillota"="#DAB977",
"Bacilli"="#BFC58A",
"Lactobacillales"="#AABD9B",

"Clostridia"="#96B4A8",
"Oscillospirales"="#82ABA4",
"Acutalibacteraceae"="#6FA19F",

"Lachnospirales"="#5B979A",
"Lachnospiraceae"="#4A8F95",

"Bacteroidota"="#46A48A",
"Bacteroidia"="#6BAED6",
"Bacteroidales"="#98BAF4",
"Bacteroidaceae"="#6794FC"
)

colourScale <- JS(sprintf(
  'd3.scaleOrdinal().domain(%s).range(%s.concat(d3.schemeCategory10))',
  jsonlite::toJSON(names(fixed_colors), auto_unbox = TRUE),
  jsonlite::toJSON(unname(fixed_colors), auto_unbox = TRUE)
))

# onRender JS
customJS <- JS("
function(el) {
  var nodeColors = {};
  d3.select(el).selectAll('.node').each(function(d){
    nodeColors[d.name] = d3.select(this).style('fill');
  });
  d3.select(el).selectAll('.link').each(function(d){
    var c = nodeColors[d.source.name];
    d3.select(this).style('stroke', c).style('stroke-opacity', 1.0);
  });
}
")

nodes <- nodes %>%
  mutate(level = case_when(
    name %in% domain_order  ~ 0,
    name %in% phylum_order  ~ 1,
    name %in% class_order   ~ 2,
    name %in% order_order   ~ 3,
    name %in% family_order  ~ 4,
    name %in% genus_order   ~ 5,
    TRUE                    ~ 6
  ))

# sankeyNetwork
p <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value  = "Value",
  NodeID = "name",
  NodeGroup = "name",
  LinkGroup = "group",
  colourScale = colourScale,
  fontFamily = "Arial",
  fontSize = 15,
  nodeWidth = 6,
  nodePadding = 12,
  iterations = 64,
  sinksRight = FALSE
)

# onRender
p <- p %>% onRender("
function(el, x) {
  var svg = d3.select(el).select('svg');
  var bb = el.getBoundingClientRect();
  var width = bb.width || 800;
  var margin = 20;
  var nodeWidth = x.options.nodeWidth || 15;
  var nodePadding = x.options.nodePadding || 10;

  var levels = Array.from(new Set(x.nodes.map(function(n){ return +n.level; }))).sort(function(a,b){return a-b;});
  if (levels.length === 0) return;

  var xScale = d3.scalePoint().domain(levels).range([margin, width - margin - nodeWidth]);

  x.nodes.forEach(function(n){
    n.x0 = xScale(+n.level);
    n.x1 = n.x0 + nodeWidth;
  });

  var levelMap = {};
  x.nodes.forEach(function(n){
    levelMap[n.level] = levelMap[n.level] || [];
    levelMap[n.level].push(n);
  });

  Object.keys(levelMap).forEach(function(l){
    var arr = levelMap[l].sort(function(a,b){
      return (b.value || 0) - (a.value || 0);
    });
    var y = 10;
    arr.forEach(function(n){
      var h = Math.max(1, n.y1 - n.y0);
      n.y0 = y;
      n.y1 = y + h;
      y += h + nodePadding;
    });
  });

  svg.selectAll('.node')
     .attr('transform', function(d){ return 'translate(' + d.x0 + ',' + d.y0 + ')'; });

  svg.selectAll('.node rect')
     .attr('width', nodeWidth)
     .attr('height', function(d){ return Math.max(1, d.y1 - d.y0); });

  svg.selectAll('.link')
     .attr('d', function(d){
        var x0 = d.source.x1;
        var x1 = d.target.x0;
        var y0 = (d.source.y0 + d.source.y1) / 2;
        var y1 = (d.target.y0 + d.target.y1) / 2;
        var xi = d3.interpolateNumber(x0, x1);
        var x2 = xi(0.5);
        return 'M' + x0 + ',' + y0 + 'C' + x2 + ',' + y0 + ' ' + x2 + ',' + y1 + ' ' + x1 + ',' + y1;
     })
     .style('stroke-opacity', 0.7);

  svg.selectAll('.node text')
     .attr('x', function(d){ return d.x1 + 6; })
     .attr('y', function(d){ return (d.y0 + d.y1) / 2; })
     .attr('text-anchor', 'start')
     .attr('dominant-baseline', 'middle');
}
")

p
saveNetwork(p, file = "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4A_taxa_sankey/gut_taxa_sankey/figure_6_sankey/domain_to_species/sankey_gut_76133.html")





# ============================================================
#Taxonomic resolution  for Fig.3B; Fig.S7
# ============================================================
library(dplyr)

INPUT_CSV <- "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4B_amp_with_taxonomy/gut/amp_with_taxonomy_gut_1973944.csv"

OUT_DIR <- "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4B_amp_with_taxonomy/gut/R_code"

OUTPUT_RANK_COUNT <- file.path(OUT_DIR, "rank_gut_amp_taxa_count.csv")
OUTPUT_GENUS_COUNT <- file.path(OUT_DIR, "genus_gut_count.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Step 1: 读取 CSV 文件
# ----------------------------
data <- read.csv(
  INPUT_CSV,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A", "na", "n/a")
)


rank_cols <- c(
  "Domain",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

for (col in rank_cols) {
  if (!col %in% colnames(data)) {
    data[[col]] <- NA_character_
  }
}

clean_tax <- function(x) {
  x <- trimws(as.character(x))

  x[x %in% c("", "NA", "N/A", "na", "n/a", "Unclassified", "unclassified")] <- NA_character_

  x[grepl("^[dpcfgos]__$", x)] <- NA_character_

  return(x)
}

data <- data %>%
  mutate(
    across(all_of(rank_cols), clean_tax),
    matched_contig = clean_tax(matched_contig),
    match_status = trimws(as.character(match_status))
  )


matched_data <- data %>%
  filter(
    match_status == "matched",
    !is.na(matched_contig),
    matched_contig != ""
  )

matched_data <- matched_data %>%
  mutate(
    rank = case_when(
      !is.na(Species) ~ "species",
      !is.na(Genus)   ~ "genus",
      !is.na(Family)  ~ "family",
      !is.na(Order)   ~ "order",
      !is.na(Class)   ~ "class",
      !is.na(Phylum)  ~ "phylum",
      !is.na(Domain)  ~ "domain",
      TRUE            ~ NA_character_
    ),
    scientific_name = coalesce(
      Species,
      Genus,
      Family,
      Order,
      Class,
      Phylum,
      Domain
    )
  )


rank_levels <- c(
  "domain",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)

rank_result <- matched_data %>%
  filter(!is.na(rank), rank != "") %>%
  mutate(rank = factor(rank, levels = rank_levels)) %>%
  group_by(rank) %>%
  summarise(
    amp_count = n(),
    .groups = "drop"
  ) %>%
  arrange(rank)

write.csv(
  rank_result,
  OUTPUT_RANK_COUNT,
  row.names = FALSE
)

print(rank_result)
cat("Rank count saved to:", OUTPUT_RANK_COUNT, "\n")


genus_result <- matched_data %>%
  filter(
    !is.na(Genus),
    Genus != ""
  ) %>%
  group_by(scientific_name = Genus) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

write.csv(
  genus_result,
  OUTPUT_GENUS_COUNT,
  row.names = FALSE
)

print(head(genus_result, 20))
cat("Genus count saved to:", OUTPUT_GENUS_COUNT, "\n")





# =========================================================================================
#Heatmap of age-associated changes in skin AMP-family abundance for Fig.3C,E; Fig.S8A,C,E
# =========================================================================================

library(tidyverse)
library(pheatmap)

decreasing_clusters_df <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/skin_heatmap/skin_top500_clusters_across_timepoints.csv",
  stringsAsFactors = FALSE
)

target_clusters <- decreasing_clusters_df %>%
  pull(family_id) %>%      
  as.character() %>%             
  head(500)                     

cat("Number of target clusters loaded:", length(target_clusters), "\n")


target_Species_df <- read.csv(
  "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/skin_heatmap/skin_top40_species_across_timepoints.csv",
  stringsAsFactors = FALSE
)


target_Species  <- target_Species_df %>%
  pull(Species) %>%        
  as.character() %>%             
  head(30)                     


meta <- read.csv("/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/skin_heatmap/metadata_skin.csv", stringsAsFactors = FALSE) 
final_data_filter <- read.csv("/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/skin_heatmap/merge_amp_with_cluster_reduced8aa_skin.csv", stringsAsFactors = FALSE)

tp <- "0-0.5years"

meta_clean <- meta %>% mutate(Age_year_group = trimws(Age_year_group))
samples_tp <- meta_clean %>% filter(Age_year_group %in% tp) %>% pull(Run)
cat("Number of samples_tp loaded:", length(samples_tp), "\n")

data_use <- final_data_filter %>%
  filter(
    sample_id %in% samples_tp,
    !is.na(family_id),
    Species %in% target_Species
  )

taxon_cluster_expr <- data_use %>%
  group_by(Species, family_id) %>%
  summarise(total_RPKM = sum(RPKM, na.rm = TRUE), .groups = "drop")

mat_full <- taxon_cluster_expr %>%
  pivot_wider(
    names_from = family_id,
    values_from = total_RPKM,
    values_fill = 0
  ) %>%
  column_to_rownames("Species") %>%
  as.matrix()

mat_complete <- matrix(0, nrow = length(target_Species), ncol = length(target_clusters))
rownames(mat_complete) <- target_Species
colnames(mat_complete) <- as.character(target_clusters)

# 提取实际存在的行列交集
common_rows <- intersect(target_Species, rownames(mat_full))
common_cols <- intersect(as.character(target_clusters), colnames(mat_full))

if (length(common_rows) > 0 && length(common_cols) > 0) {
  mat_complete[common_rows, common_cols] <- mat_full[common_rows, common_cols]
}

mat_log <- log10(mat_complete + 1)

cat("Matrix dimensions:", dim(mat_log), "\n")
cat("Value range in log10(RPKM + 1):", round(range(mat_log), 3), "\n")

out_dir <- "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/skin_heatmap/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_val <- max(mat_log, na.rm = TRUE)
breaks_vec <- seq(0, ceiling(max_val + 0.5), length.out = 101)

draw_heatmap <- function() {
  pheatmap(
    mat_log,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("black", "blue", "lightblue", "red", "orange"))(100),
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 18,
    fontsize_col = 12,
    main = "Top 30 Species × Top 500 clusters (0-0.5years)",
    border_color = NA,
    legend = TRUE,
    legend_title = expression(log[10](RPKM + 1)),
    breaks = seq(0, 9, length.out = 101)
  )
}

pdf(file.path(out_dir, "0-0.5years_decreasing_skin.pdf"), width = 12, height = 10)
draw_heatmap()
dev.off()

png(file.path(out_dir, "0-0.5years_decreasing_skin.png"), width = 3600, height = 3000, res = 300)
draw_heatmap()
dev.off()

write.csv(mat_log, file.path(out_dir, "0-0.5years_decreasing_skin.csv"))

cat("✅ Done! Output saved to:", out_dir, "\n")

#Matrix dimensions: 30 500
#Value range in log10(RPKM + 1): 0 6.135



# =========================================================================================
#taxonomic_richness in gut for Fig.3D,F; Fig.S8B,D,F
# =========================================================================================

#1)Step1:
nano gut_species_rpkm_all_timepoints.R

library(tidyverse)

meta <- read.csv("/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/heatmap/gut_heatmap_true/gut_metadata_Zeng_2022.csv", 
                 stringsAsFactors = FALSE) 

final_data_filter <- read_csv("/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/heatmap/gut_heatmap_true/merge_amp_with_cluster_reduced8aa_gut.csv") %>%
  filter(!is.na(family_id))

meta_clean <- meta %>%
  mutate(Age_year_group = trimws(Age_year_group))

time_points <- c("0-1 years", "1-2 years", "2-3 years")

get_all_species <- function(tp, meta_df, expr_df) {
  cat("Processing time point:", tp, "\n")
  
  samples_tp <- meta_df %>%
    filter(Age_year_group == tp) %>%
    pull(Sample) 
  
  if (length(samples_tp) == 0) {
    warning("No samples for time point: ", tp)
    return(tibble()) 
  }

  data_tp <- expr_df %>%
    filter(sample_id %in% samples_tp)
  
  if (nrow(data_tp) == 0) {
    warning("No expression data for time point: ", tp)
    return(tibble())
  }

  taxon_cluster_expr <- data_tp %>%
    group_by(Species, family_id) %>%         
    summarise(total_RPKM = sum(RPKM, na.rm = TRUE), .groups = "drop")

  mat_wide <- taxon_cluster_expr %>%
    pivot_wider(
      names_from = family_id,
      values_from = total_RPKM,
      values_fill = 0
    )     
  return(mat_wide)
}

results_list <- map(time_points, ~ get_all_species(.x, meta_clean, final_data_filter)) %>%
  set_names(time_points)

output_dir <- "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/taxonomic_richness/gut_richness/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

walk2(results_list, time_points, function(df, tp) {
  safe_tp_name <- str_replace_all(tp, "[^a-zA-Z0-9]", "_")  # 将文件名中的特殊符号替换为下划线
  file_path <- paste0(output_dir, "species_", safe_tp_name, ".csv")
  write.csv(df, file_path, row.names = FALSE)
  cat("已单独导出文件:", file_path, "\n")
})

all_species_data <- bind_rows(results_list, .id = "TimePoint")

output_path_all <- paste0(output_dir, "gut_species_rpkm_all_timepoints.csv")
write.csv(
  all_species_data,
  output_path_all,
  row.names = FALSE
)

cat("\n全部分析完成！总表已保存至:", output_path_all, "\n")


#1)Step2
library(data.table)
library(ggplot2)
library(rcartocolor)

infile <- "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/taxonomic_richness/gut_richness/gut_species_rpkm_all_timepoints.csv"

out_dir <- "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4CDFEF_heatmap/taxonomic_richness/gut_richness/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

time_levels <- c("0-1 years", "1-2 years", "2-3 years")
meta_cols <- c("TimePoint", "Species")

dt <- fread(
  infile,
  na.strings = c("", "NA", "NaN"),
  showProgress = TRUE,
  nThread = max(1L, parallel::detectCores() - 1L)
)

missing_meta <- setdiff(meta_cols, names(dt))
if (length(missing_meta) > 0) {
  stop("The input file is missing these columns: ", paste(missing_meta, collapse = ", "))
}

value_cols <- setdiff(names(dt), meta_cols)

if (length(value_cols) == 0) {
  stop("The RPKM value column for calculating richness was not found")
}


is_num <- vapply(dt[, ..value_cols], is.numeric, logical(1))

if (any(!is_num)) {
  non_num_cols <- value_cols[!is_num]
  message(
    "numeric: ",
    paste(head(non_num_cols, 10), collapse = ", "),
    ifelse(length(non_num_cols) > 10, " ...", "")
  )
  
  dt[, (non_num_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))),
     .SDcols = non_num_cols]
}


chunk_size <- 500L
richness_count <- integer(nrow(dt))

for (i in seq(1L, length(value_cols), by = chunk_size)) {
  cols_i <- value_cols[i:min(i + chunk_size - 1L, length(value_cols))]
  
  mat_i <- as.matrix(dt[, ..cols_i])
  
  richness_count <- richness_count + rowSums(!is.na(mat_i) & mat_i > 0)
  
  rm(mat_i)
  gc(verbose = FALSE)
}

summary_table <- dt[, .(TimePoint, Species)]
summary_table[, count_gt_above_0 := richness_count]
summary_table[, TimePoint := factor(TimePoint, levels = time_levels)]

summary_table <- summary_table[!is.na(TimePoint)]

fwrite(
  summary_table,
  file.path(out_dir, "gut_taxonomy_richness_summary.csv")
)

rm(dt)
gc(verbose = FALSE)

present_levels <- time_levels[time_levels %in% as.character(unique(summary_table$TimePoint))]

if (length(present_levels) < 2) {
  stop("there are less than two TimePoints")
}

safe_wilcox <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }
  
  if (length(unique(c(x, y))) < 2) {
    return(1)
  }
  
  suppressWarnings(
    wilcox.test(x, y, exact = FALSE)$p.value
  )
}

comparisons <- combn(present_levels, 2, simplify = FALSE)

pval_results <- rbindlist(lapply(comparisons, function(comp) {
  d1 <- summary_table[TimePoint == comp[1], count_gt_above_0]
  d2 <- summary_table[TimePoint == comp[2], count_gt_above_0]
  
  data.table(
    group1 = comp[1],
    group2 = comp[2],
    p.value = safe_wilcox(d1, d2)
  )
}))

pval_results[, p.adj := p.adjust(p.value, method = "BH")]

pval_results[, sig_label := fcase(
  is.na(p.adj), NA_character_,
  p.adj < 0.001, "***",
  p.adj < 0.01,  "**",
  p.adj < 0.05,  "*",
  default = "ns"
)]

fwrite(
  pval_results,
  file.path(out_dir, "gut_taxonomy_richness_pairwise_wilcox.csv")
)


y_min <- min(summary_table$count_gt_above_0, na.rm = TRUE)
y_max <- max(summary_table$count_gt_above_0, na.rm = TRUE)

y_step <- max(1, (y_max - y_min) * 0.08)
bracket_tip <- y_step * 0.25

stat_annotations <- pval_results[!is.na(p.adj) & sig_label != "ns"]

if (nrow(stat_annotations) > 0) {
  stat_annotations[, `:=`(
    x = match(group1, time_levels),
    xend = match(group2, time_levels),
    y.position = y_max + seq_len(.N) * y_step
  )]
} else {
  stat_annotations <- data.table(
    x = numeric(),
    xend = numeric(),
    y.position = numeric(),
    sig_label = character()
  )
}

p <- ggplot(
  summary_table,
  aes(x = TimePoint, y = count_gt_above_0, fill = TimePoint)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = TimePoint),
    width = 0.15,
    size = 0.35,
    alpha = 0.20,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 2.5,
    fill = "white"
  ) +
  
  geom_segment(
    data = stat_annotations,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  
  geom_segment(
    data = stat_annotations,
    aes(x = x, xend = x, y = y.position - bracket_tip, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  
  geom_segment(
    data = stat_annotations,
    aes(x = xend, xend = xend, y = y.position - bracket_tip, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  
  geom_text(
    data = stat_annotations,
    aes(x = (x + xend) / 2, y = y.position, label = sig_label),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 5
  ) +
  
  scale_fill_carto_d(palette = "NCS") +
  scale_color_carto_d(palette = "NCS") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.18))) +
  labs(
    title = "Gut Taxonomic Richness Across Time Points",
    x = "Time Point",
    y = "Taxonomic richness, count of RPKM columns > 0",
    fill = "Time Point"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 12),
    axis.ticks = element_line(),
    panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(
  file.path(out_dir, "gut_taxonomy_richness.pdf"),
  plot = p,
  width = 8,
  height = 6
)

ggsave(
  file.path(out_dir, "gut_taxonomy_richness.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\n taxonomy richness ", out_dir, "\n")



#3)step3:for top30

library(ggsci)
library(tidyverse)
library(rcartocolor)
library(rstatix)
library(ggpubr)

summary_table <- read_csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_4CDEF_heatmap/gut_richness/gut_taxonomy_richness_summary.csv")

taxanomy_top40 <- read_csv("C:/Users/taoch/Desktop/632751_clusters_文章/Fig_4CDEF_heatmap/gut_richness/gut_top40_species_across_timepoints.csv") 

taxanomy_top30 <- taxanomy_top40 %>% slice_head(n = 30)
top30_species_names <- taxanomy_top30$Species


plot_data <- summary_table %>%
  filter(Species %in% top30_species_names) %>%        
  select(Species, TimePoint, count_gt_above_0) %>%   
  mutate(
    Species = factor(Species, levels = top30_species_names),
    TimePoint = factor(TimePoint, levels = c("0-1 years", "1-2 years", "2-3 years")) 

stat_test <- plot_data %>%
  wilcox_test(count_gt_above_0 ~ TimePoint, p.adjust.method = "BH") %>% 
  add_significance("p.adj") %>%                      
  add_xy_position(x = "TimePoint", dodge = 0)        


p <- ggplot(plot_data, aes(x = TimePoint, y = count_gt_above_0, fill = Species)) +
  geom_col(position = "stack") +                     
  stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  labs(
    title = "Top 30 Species richness",
    x = "TimePoint",
    y = "Richness (the count of family AMP rpkm >0)",
    fill = "Species"
  ) +
  scale_fill_carto_d(palette = "NCS") +              
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    plot.title = element_text(hjust = 0.5),         
    legend.position = "right",
    legend.title = element_text(size = 12),
    axis.ticks = element_line(),
    panel.border = element_rect(colour = "black", linewidth = 0.6, fill = NA),
    panel.grid = element_blank()                    


# =========================================================================================
#AMP production changes with host age for Fig.3G,H
# =========================================================================================
library(dplyr)
library(ggplot2)

early_metadat <- read.csv("/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4GHIJ_per_species/early_top_3/early_top_3_spcies_sum_rpkm_per_AMP_and_dol_select.csv")

target_species <- c("Escherichia coli", "Enterococcus faecalis", "Bifidobacterium longum")

lapply(target_species, function(spc) {
  
  spc_data <- early_metadat %>%
    filter(Species == spc)
  
  safe_name <- gsub(" ", "_", spc) 
  write.csv(spc_data, 
            file.path("/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4GHIJ_per_species/early_top_3", 
                      paste0("early_metadat_", safe_name, ".csv")), 
            row.names = FALSE)
  
  p <- ggplot(spc_data, aes(x = DOL, y = log_total_rpkm)) +
    geom_point(alpha = 0.6, size = 0.2, color = "#299d8f") +
    geom_smooth(method = "lm", se = TRUE, color = "#f3a361") +
    scale_x_continuous(  
    breaks = c(0, 400, 800, 1200), 
    limits = c(0, 1200)
     ) +
    scale_y_continuous(  
    breaks = c(0, 5, 10, 15, 20), 
    limits = c(0, 20)
     ) +
    labs(
      title = paste("Early Colonizing Bacteria:", spc),
      x = "DOL",
      y = "log10(total RPKM+1)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(color = "black", size = 12)
    )
  
  # 3.4 打印并保存图表
ggsave(
  filename = paste0("early_plot_", safe_name, ".png"),
  plot = p,
  path = "/mnt/data/taochunlin/project/amp/R_new_for_632751cluster/Fig_4GHIJ_per_species/early_top_3",
  width = 6, 
  height = 4,
  dpi = 300  
)
})

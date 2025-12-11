install.packages("BiocManager")
BiocManager::install("ggplot2")
BiocManager::install("dplyr")
BiocManager::install("fgsea")
BiocManager::install("ggrepel")

# Load libraries
library(ggplot2)
library(dplyr)
library(fgsea)
library(ggrepel)

# 1. Simulate CRISPR screen data
set.seed(123)

# Create gene list and sgRNA data
genes <- paste0("Gene", 1:1000)
guides_per_gene <- 10
sgRNA_data <- data.frame(
  gene = rep(genes, each = guides_per_gene),
  sgRNA = paste0("sg", 1:(1000 * guides_per_gene)),
  count_control = rpois(1000 * guides_per_gene, lambda = 100),
  count_treated = rpois(1000 * guides_per_gene, lambda = 100)
)

# Define hypothesis genes
chemo_genes <- c("CCL2", "CCR2", "CXCL12", "CXCR4")
immune_genes <- c("PDL1", "PTPN11")
tf_brakes <- c("FOXP3", "STAT3")

# Replace some generic gene names with hypothesis names
genes[1:4] <- chemo_genes
genes[5:6] <- immune_genes
genes[7:8] <- tf_brakes
sgRNA_data$gene <- rep(genes, each = guides_per_gene)

# Assign effect sizes
sgRNA_data$effect <- rnorm(nrow(sgRNA_data), 0, 0.5)  # baseline noise
sgRNA_data$effect[sgRNA_data$gene %in% chemo_genes] <- rnorm(sum(sgRNA_data$gene %in% chemo_genes), 1.8, 0.2)  # strong positive
sgRNA_data$effect[sgRNA_data$gene %in% immune_genes] <- rnorm(sum(sgRNA_data$gene %in% immune_genes), 1.5, 0.2)  # positive
sgRNA_data$effect[sgRNA_data$gene %in% tf_brakes] <- rnorm(sum(sgRNA_data$gene %in% tf_brakes), -1.8, 0.2)       # strong negative

# Aggregate to gene-level scores
gene_scores <- sgRNA_data %>%
  group_by(gene) %>%
  summarize(mean_effect = mean(effect), .groups = 'drop')

# Assign p-values based on effect size (simulate significance)
gene_scores$pval <- pmin(0.05, exp(-abs(gene_scores$mean_effect)))  # smaller p for larger effects

# Simulate logFC values
logFC <- rnorm(length(genes), 0, 0.5)


# 2. Aggregate to gene-level scores
gene_scores <- sgRNA_data %>%
  group_by(gene) %>%
  summarize(mean_effect = mean(effect), .groups = 'drop')

# 3. QC PLOTS

# sgRNA DISTRIBUTION
ggplot(sgRNA_data, aes(effect)) +
  geom_histogram(bins = 50, fill = "lightblue") +
  ggtitle("sgRNA Effect Distribution")

# SAMPLE NORMALIZATION HISTOGRAM
# Compute normalization values
sgRNA_data$log_ratio <- log2(sgRNA_data$count_treated + 1) - log2(sgRNA_data$count_control + 1)

# Add condition labels
sgRNA_data$condition <- rep(c("Tumor", "Non-Tumor"), each = nrow(sgRNA_data)/2)

# Plot overlay histogram
ggplot(sgRNA_data, aes(x = log_ratio, fill = condition)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
  scale_fill_manual(values = c("Tumor" = "lightblue", "Non-Tumor" = "orange")) +
  ggtitle("Normalization Histogram: Tumor vs Non-Tumor") +
  theme_minimal() +
  labs(x = "Log2 Normalized Ratio", y = "Count", fill = "Condition")






# 4. OUTPUT PLOTS

# GENE RANK

# --- Assume sgRNA_data and gene_scores are already created from your simulation code ---

# Add logFC values to gene_scores
gene_scores$logFC <- rnorm(nrow(gene_scores), 0, 0.5)  # simulated logFC


gene_scores <- gene_scores %>%
  arrange(desc(logFC)) %>%           # put strongest first; flip if lower=stronger
  mutate(
    rank = row_number(),
    pval = rank / (n() + 1),         # uniform order statistics → small p for top ranks
    padj = p.adjust(pval, method = "BH")
  )

# Order genes by p-value
gene_scores <- gene_scores %>%
  arrange(pval) %>%
  mutate(rank = row_number())

# Define hypothesis genes for labeling
hypothesis_genes <- c(chemo_genes, immune_genes, tf_brakes)

# Create the gene rank plot with logFC as size and p-value as color
ggplot(gene_scores, aes(x = rank, y = mean_effect, size = abs(logFC), color = pval)) +
  geom_point(alpha = 0.7) +
  geom_point(data = gene_scores %>% filter(gene %in% hypothesis_genes),
             aes(x = rank, y = mean_effect, size = abs(logFC), color = pval),
             shape = 21, fill = "red", stroke = 1.2) +
  geom_text(data = gene_scores %>% filter(gene %in% hypothesis_genes),
            aes(label = gene),
            vjust = -1, color = "red", size = 3.5) +
  scale_size_continuous(name = "logFC (absolute)") +
  scale_color_gradient(low = "blue", high = "orange", name = "p-value") +
  labs(title = "Gene Rank Plot by p-value",
       x = "Gene Rank (ordered by p-value)",
       y = "Mean Effect Size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



-------------
  

# ✅ Step 1: Rank genes by p-value
gene_scores <- gene_scores %>%
  arrange(pval) %>%            # Sort by p-value (ascending)
  mutate(rank = row_number())  # Assign rank after sorting

# ✅ Step 2: Subset hypothesis genes
hypothesis_df <- gene_scores %>%
  filter(gene %in% hypothesis_genes)


# ✅ Step 3: Plot
ggplot(gene_scores, aes(x = rank, y = logFC)) +
  geom_point(aes(size = abs(logFC), color = pval), alpha = 0.7) +
  geom_point(data = hypothesis_df, color = "red", size = 4) +
  geom_text_repel(
    data = hypothesis_df,
    aes(label = gene),
    color = "black", 
    size = 6, 
    fontface = "bold", 
    max.overlaps = Inf,
    bg.color = "white",
    bg.r = 0.15
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +   # Horizontal reference
  geom_vline(xintercept = 100, linetype = "dashed", color = "gray40") + # Example rank cutoff
  scale_color_gradient(low = "darkblue", high = "lightyellow", name = "P-value") +
  scale_size_continuous(range = c(2, 6), name = "Log2 Fold Change") +
  theme_minimal() +
  theme(plot.margin = margin(20, 40, 20, 40)) +
  coord_cartesian(clip = "off") +
  ggtitle("Gene Rank Plot") +
  labs(x = "Ranking of differentially expressed genes (by p-value)", y = "Log2 Fold Change")

# VOLCANO PLOT
pvals <- runif(nrow(gene_scores), 0, 0.05) # Simulated p-values
volcano_data <- data.frame(gene_scores, pval = pvals)
ggplot(volcano_data, aes(x = mean_effect, y = -log10(pval))) +
  geom_point() + ggtitle("Volcano Plot")

# GSEA
pathways <- list(ChemoAxis = chemo_genes, ImmuneCheckpoints = immune_genes, TFBrakes = tf_brakes)
ranks <- setNames(gene_scores$mean_effect, gene_scores$gene)
fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize = 1, maxSize = 500)
plotEnrichment(pathways$ChemoAxis, ranks)
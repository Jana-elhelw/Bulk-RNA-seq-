#running fgsea following bioconductor steps

# Step 1: Install and load msigdbr
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

library(msigdbr)

# Step 2: Download GO Biological Process gene sets for human
msig_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

# View the first few rows
head(msig_go)

# Step 3: Convert to fgsea-compatible list
pathways_list <- split(msig_go$gene_symbol, msig_go$gs_name)

# Check the first few gene sets
str(pathways_list[1:2], max.level = 1)

# Rebuild the gene_ranks vector (in case R was restarted)
res_clean_pou3f2_sox1 <- na.omit(res_pou3f2_sox1[, c("SYMBOL", "log2FoldChange")])
ranks_df <- res_clean_pou3f2_sox1[order(abs(res_clean_pou3f2_sox1$log2FoldChange), decreasing = TRUE), ]
ranks_df <- ranks_df[!duplicated(ranks_df$SYMBOL), ]
gene_ranks <- ranks_df$log2FoldChange
names(gene_ranks) <- ranks_df$SYMBOL
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

library(fgsea)

# Re-run Step 4 using fgseaMultilevel (no nperm)
set.seed(123)
fgsea_res <- fgsea(
  pathways = pathways_list,
  stats    = gene_ranks,
  minSize  = 15,
  maxSize  = 500
)

# View top results
head(fgsea_res[order(padj)][, .(pathway, padj, NES)])

# Save full fgsea results
# Convert fgsea result to a data frame
fgsea_res_df <- as.data.frame(fgsea_res)

# Sort by adjusted p-value
fgsea_res_df <- fgsea_res_df[order(fgsea_res_df$padj), ]

# Remove 'leadingEdge' column (it's a list and can't be written to CSV)
fgsea_res_df_no_list <- fgsea_res_df[, !(names(fgsea_res_df) %in% "leadingEdge")]

# Save to CSV
write.csv(fgsea_res_df_no_list, "results/enrichment/FGSEA_GO_BP_POU3F2_vs_SOX1.csv", row.names = FALSE)


library(ggplot2)

# Get the name of the top enriched pathway
top_pathway <- fgsea_res_df$pathway[1]

# Plot the enrichment curve
png("results/enrichment/FGSEA_enrichment_plot_top.png", width = 1200, height = 900, res = 150)
plotEnrichment(pathways_list[[top_pathway]], gene_ranks) +
  labs(title = top_pathway)
dev.off()


dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

#-----------------------------
#step 1: install/load packages
#-----------------------------
# Install packages only if needed:
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "fgsea", "org.Hs.eg.db"))

library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
library(cowplot)


#-----------------------------------------
#step 2: load POU3F2 vs SOX1 results table
#-----------------------------------------
#pou3f2 vs sox1 had the biggest DE contrast

res_pou3f2_sox1 <- read.csv("results/DEGs/POU3F2_vs_SOX1_DEGs.csv")


#----------------------
#Step 3. Run ORA (GO)
#---------------------
#Pick significant DE genes (padj < 0.05):

# 1) get significant gene symbols
sig_genes <- na.omit(res_pou3f2_sox1$SYMBOL
                     [res_pou3f2_sox1$padj < 0.05])

# 2) run GO enrichment
ego <- enrichGO(gene          = sig_genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",          # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# 3) preview results
head(ego)

# 4) clean the enrichment table
ego_df <- as.data.frame(ego)

# keep useful columns, sort by FDR (p.adjust)
ego_clean <- ego_df[, c("ID","Description","GeneRatio","BgRatio",
                        "pvalue","p.adjust","qvalue","Count","geneID")]
ego_clean <- ego_clean[order(ego_clean$p.adjust), ]

# keep top 30
ego_top30 <- head(ego_clean, 30)

# save both full and top30
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
write.csv(ego_clean,  "results/enrichment/GO_BP_POU3F2_vs_SOX1_full.csv", row.names = FALSE)
write.csv(ego_top30,  "results/enrichment/GO_BP_POU3F2_vs_SOX1_top30.csv", row.names = FALSE)

# 5) Dotplot of top pathways

library(enrichplot)
library(ggplot2)

# top 10 terms, ordered by FDR
p <- dotplot(ego, showCategory = 10, orderBy = "p.adjust") +
  ggtitle("GO BP enrichment — POU3F2 vs SOX1 (padj)") +
  theme(plot.title = element_text(hjust = 0.5))

# view
print(p)

# save
ggsave("results/enrichment/GO_BP_dotplot_top10.png", p, width = 10, height = 7, dpi = 150)

#key findings: Ribosome-related processes (ribosome biogenesis, rRNA processing, translation) 
#are strongly enriched → suggests differences in protein synthesis machinery 
#between POU3F2 and SOX1 cells.
#Cell division processes (mitotic sister chromatid segregation, 
#mitotic nuclear division) are enriched → may mean one condition has 
#higher proliferation activity.
#Forebrain development shows up → directly links to neural 
#differentiation, consistent with transcription factors (SOX1 is a 
#neural progenitor marker, POU3F2 is neurogenic).

#---------------------------
#Step 4. Create a ranked gene list (for fgsea)
#--------------------------
#use all genes ranked by log2foldchange

# keep rows with both SYMBOL and log2FC
res_clean_pou3f2_sox1 <- na.omit(res_pou3f2_sox1[, c("SYMBOL",
                                                     "log2FoldChange")])

# named vector: names = gene symbols, values = log2FC
gene_ranks <- res_clean_pou3f2_sox1$log2FoldChange
names(gene_ranks) <- res_clean_pou3f2_sox1$SYMBOL

# order decreasing
gene_ranks <- sort(gene_ranks, decreasing = TRUE)


#----------------------
# Step 5. Run fgsea
#----------------------
#running fgsea following bioconductor steps

# 1: Install and load msigdbr

library(msigdbr)

# 2: Download GO Biological Process gene sets for human
msig_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

# View the first few rows
head(msig_go)

# 3: Convert to fgsea-compatible list
pathways_list <- split(msig_go$gene_symbol, msig_go$gs_name)

# Check the first few gene sets
str(pathways_list[1:2], max.level = 1)

# 4: Rebuild the gene_ranks vector (in case R was restarted)
res_clean_pou3f2_sox1 <- na.omit(res_pou3f2_sox1[, c("SYMBOL", "log2FoldChange")])
ranks_df <- res_clean_pou3f2_sox1[order(abs(res_clean_pou3f2_sox1$log2FoldChange), decreasing = TRUE), ]
ranks_df <- ranks_df[!duplicated(ranks_df$SYMBOL), ]
gene_ranks <- ranks_df$log2FoldChange
names(gene_ranks) <- ranks_df$SYMBOL
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# 5: run fgsea and save to a csv file
library(fgsea)

# Run fgseaMultilevel (no nperm)
set.seed(123)
fgsea_res <- fgsea(
  pathways = pathways_list,
  stats    = gene_ranks,
  minSize  = 15,
  maxSize  = 500
)

# View top results
head(fgsea_res[order(padj)][, .(pathway, padj, NES)])

# key findings: [+ve NES are pathways upregulated in POU3F2, 
#-VE NES are pathways upregulated in SOX1] So, central nervous system neuron 
#differentation is upregulated in POU3F2 while regulation of cell killing in SOX1

# Save full fgsea results
# Convert fgsea result to a data frame
fgsea_res_df <- as.data.frame(fgsea_res)

# Sort by adjusted p-value
fgsea_res_df <- fgsea_res_df[order(fgsea_res_df$padj), ]

# Remove 'leadingEdge' column (it's a list and can't be written to CSV)
fgsea_res_df_no_list <- fgsea_res_df[, !(names(fgsea_res_df) %in% "leadingEdge")]

# Save to CSV
write.csv(fgsea_res_df_no_list, "results/enrichment/FGSEA_GO_BP_POU3F2_vs_SOX1.csv", row.names = FALSE)


# 6: plot results
#Enrichment curve for the top upregulated pathway in POU3F2 (to visually understand how the genes in the top pathway
#are distributed in our ranked gene list)

library(ggplot2)

# Get the name of the top enriched pathway
top_pathway <- fgsea_res_df$pathway[1]

# Plot the enrichment curve
png("results/enrichment/FGSEA_enrichment_plot_top.png", width = 1200, height = 900, res = 150)
plotEnrichment(pathways_list[[top_pathway]], gene_ranks) +
  labs(title = top_pathway)
dev.off()

# key finding: genes in this pathway are upregulated in POU3F2 and enriched near the top of the ranked list


#Enrichment curve for the top upregulated pathway in SOX1 (so top downregulated in POU3F2)

# Select top downregulated pathway (lowest padj among NES < 0)
top_down_pathway <- fgsea_res_df[fgsea_res_df$NES < 0, ]
top_down_pathway <- top_down_pathway[order(top_down_pathway$padj), ]
top_down_pathway_name <- top_down_pathway$pathway[1]
# Plot top pathway upregulated in SOX1
png("results/enrichment/FGSEA_enrichment_plot_top_SOX1_upregulated.png", width = 1200, height = 900, res = 150)
plotEnrichment(pathways_list[[top_down_pathway_name]], gene_ranks) +
  labs(title = top_down_pathway_name)
dev.off()

# key finding: Genes involved in cell killing (e.g., apoptosis regulation) are 
#more active in SOX1-expressing cells

# Pre-processing and DESEQ2 Analysis (downstream analysis)
# ============================
# RNA-seq Metadata Preparation
# ============================

#make a folder layout
dir.create("data",    showWarnings = FALSE)
dir.create("scripts", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/DEGs/plots",  showWarnings = FALSE)

# Clear environment - to prevent old variables from interfering with the analysis
rm(list = ls())

# Load libraries : for data manipulation as we will use mutate() and filter()
library(dplyr)

# ----------------------------
# Step 1: Load the count matrix (the raw data DESeq2 will use to compare gene expression.)
# ----------------------------
# Use read.delim if the file is tab-delimited (it was eventhough the file was.csv)
counts <- read.delim("C:/R_Projects/Bulk RNA seq/data/counts.csv",
                     row.names = 1, check.names = FALSE)

# Quick check
dim(counts)   # Expect ~36601 rows × 39 columns

# --------------------------------
# Step 2: Load and parse metadata
# --------------------------------
meta_raw <- readLines("C:/R_Projects/Bulk RNA seq/data/GSE303918_series_matrix.txt")

# Extract sample characteristics lines
grep("!Sample_", meta_raw, value = TRUE) #we can see that in line 10 there was 
#"sample_characteristics_ch1" which contains the sample info we want
characteristics <- grep("!Sample_characteristics_ch1", meta_raw, value = TRUE)


# Clean (remove labels) and split each line into fields (structured rows)
characteristics_clean <- lapply(characteristics, function(line) {
  fields <- strsplit(sub("!Sample_characteristics_ch1\t", "", line), "\t")[[1]]
  return(fields)
})

#check if rows represent samples, if NOT then transpose in next step
lengths <- sapply(characteristics_clean, length)
print(lengths)
#what to look for: output showed list of 8 with 39 values each, 
#that means each line represents one characteristic across 39 samples, 
#meaning 8 rows × 39 columns, so need to transpose to get 39 rows × 8 columns.


# build the metadata matrix: combine all cleaned lines into a dataframe and transpose it
meta_matrix <- as.data.frame(t(do.call(rbind, characteristics_clean)),
                             stringsAsFactors = FALSE)
#(now we have 1 row/sample and characteristics are in columns)

# Add sample IDs (extract from the meta_raw, remove labels, add as new column in meta_matrix)
sample_ids <- grep("!Sample_geo_accession", meta_raw, value = TRUE)
sample_ids <- strsplit(sub("!Sample_geo_accession\t", "", sample_ids), "\t")[[1]]
meta_matrix$sample_id <- sample_ids

# Check dimensions
dim(meta_matrix)   # Expect 39 rows × 9 columns

# --------------------------------
# Step 3: Filter for day 14 samples (to know which TF to use in the analysis)
# --------------------------------
#in the meta_matrix, v6 was the characteristic for time so we filter by it 
#and extract only the samples in day 14
day14_samples <- meta_matrix[grepl("day 14", meta_matrix$V6), ]

# Check dimensions
dim(day14_samples) # Expect 6 rows (3 KO and 3 WT)× 9 columns

# --------------------------------
# Step 4: See which TFs are present
# --------------------------------
table(day14_samples$V4)
# Output was:
# "altered gene_symbol: EMX2"           3
# "altered gene_symbol: Not applicable" 3 (these are wildtype)

#---------------
#[will do the next steps as the column names in counts are not GSM.., 
#so they don't match the sampleIDs in the meta_matrix]
#so reconstructing a new metadata frame (char_df) that produces a “match label”, 
#which can be matched to the sample names in the counts matrix.
#SO BASICALLY, WILL CREATE A MATCH_LABEL COLUMN INSTEAD OF THE SAMPLE ID ONE

# Step 4.5: Build char_df and match_label (char_df is the new meta_matrix)
# --------------------------------
# Clean and Restructure Metadata from GEO
.clean <- function(x) trimws(gsub("\"", "", x))

# Extract characteristics from meta_raw again
char_lines <- grep("^!Sample_characteristics_ch1", meta_raw, value = TRUE)
# Build a list of the values for each sample
char_list <- lapply(char_lines, function(line) {
  fields <- strsplit(sub("^!Sample_characteristics_ch1\t", "", line), "\t")[[1]]
  .clean(fields)
})
# You now have a list where each element is a characteristic (e.g., time, genotype, cell line) across all samples.

# Keys and values
#extract the column names (keys) from the meta_raw
get_key <- function(x) sub("^([^:]+):.*$", "\\1", x[1])
keys <- vapply(char_list, get_key, character(1))
#char_list is a list of character vectors, where each vector holds the values for one 
#characteristic (like "altered gene_symbol: EMX2", etc.)
#We want to extract just the key name from the first element — e.g., get "altered gene_symbol" from "altered gene_symbol: EMX2"
#so the regex is used, and then the keys vector will become the column names in the final metadata frame

#extract the row names (values)
strip_values <- function(v) sub("^[^:]+:\\s*", "", v)
vals_list <- lapply(char_list, strip_values)
#we want to remove from char_list the part before the colon, and keep just the values like "EMX2" or "Not applicable"
#vals_list becomes a list of just the values for each characteristic across all samples

#then we combine these clean values into a matrix where rows are characteristics and columns are samples
#then transpose
# Build tidy metadata
char_mat <- do.call(rbind, vals_list)
rownames(char_mat) <- make.names(keys, unique = TRUE)
char_df <- as.data.frame(t(char_mat), stringsAsFactors = FALSE)
#now char_df is the clean and complete metadata table

# Add sample titles
titles_raw <- grep("^!Sample_title", meta_raw, value = TRUE)
titles <- .clean(strsplit(sub("!Sample_title\t", "", titles_raw), "\t")[[1]])
char_df$Sample_title <- titles

# Create match_label to align with counts colnames
#create a label (called match_label) for each sample in the metadata 
#that matches the format used as column names in the counts table.
# So we want to recreate the EMX2_PTC_AB2_Day14 (colnames) part using metadata, so we can match it.

# 1-genotype label
gene_label <- ifelse(char_df$altered.gene_symbol == "Not applicable", "WT", char_df$altered.gene_symbol)
#result: gene_label = EMX2 or WT

# 2-extract the well name from the cell line
well <- sub(".*-", "", char_df$cell.line)
#result: will take what is after the dash, so H9-AB2 will be AB2

# 3-format the time
day <- gsub("day ", "Day", char_df$time)
# will capitalize day

#4-put it all together: combine all parts into 1 string
char_df$match_label <- paste0(gene_label, "_PTC_", well, "_", day)
#result (the match label): "EMX2" + "_PTC_" + "AB2" + "_Day14" = "EMX2_PTC_AB2_Day14"

#now the match label values match the beginning of the column names 
#in counts making it easy to align metadata to the counts data

# --------------------------------
# Step 5: Align metadata to counts- match the match_label to the column names
# --------------------------------
#Get the part before _GT from count column names
count_prefix <- sub("_GT.*", "", colnames(counts))

#Filter counts to just samples present in metadata:
keep_cols <- count_prefix %in% char_df$match_label
counts <- counts[, keep_cols]

#Reorder metadata to match the order of the counts:
meta_aligned <- char_df[match(sub("_GT.*", "", colnames(counts)), char_df$match_label), ]

#check if they align
stopifnot(identical(meta_aligned$match_label, sub("_GT.*", "", colnames(counts))))
#If this runs without error, then you’ve successfully matched the metadata to the counts!

#confirm that the final metadata (meta_aligned) and counts have the same sample order
head(meta_aligned$match_label)
head(sub("_GT.*", "", colnames(counts)))

# --------------------------------
# Step 6: Create coldata 
# --------------------------------
coldata <- meta_aligned

#based on the meta_aligned table, must change analysis plan as the WT samples are not in counts file
# New Analysis plan: work on the day 11 TFs (HES5, POU3F2, SOX1)

dim(counts)
dim(coldata)
head(colnames(counts))
head(rownames(coldata))
# set metadata rownames to the sample names (counts colnames)
rownames(coldata) <- colnames(counts)

# quick check (should be TRUE)
identical(rownames(coldata), colnames(counts))

# --------------------------------------
#step 7: create the condition col in coldata
# ---------------------------------------
# we now have ready the 2 ingredients DESeq2 needs: counts and coldata
#we want to compare different KO groups in day 11
#so will give DEseq a clean condition variable to split samples into groups

# 1) make condition = TF name 
coldata$condition <- factor(coldata$altered.gene_symbol)
#create a new column in the coldata with the TF names (altered.gene_symbol)
#wrap it in factor as DESeq expects categorical groups as factors

# 2) keep only Day 11 samples (case-insensitive) by filtering the time col.
is_day11 <- grepl("^\\s*day\\s*11\\s*$", tolower(coldata$time))
coldata  <- droplevels(coldata[is_day11, , drop = FALSE])
counts   <- counts[, rownames(coldata), drop = FALSE]

# 3) sanity check: how many per group?
table(coldata$condition)

#now we have a grouping variable "condition"

# ---------------------
# step 8: build and run DESeq object (dds)
# ---------------------
library (DESeq2)

#the counts table contained decimals from formatting, so are stored as numeric not integers
# convert to integer storage (safe because they are raw counts)
counts <- round(as.matrix(counts))
storage.mode(counts) <- "integer"

# quick check: first few values should be integers now
counts[1:5, 1:3]

# 1) build DESeqDataSet (design by KO group)
dds <- DESeqDataSetFromMatrix(
  countData = counts,     # raw counts (Day 11 only)
  colData   = coldata,    # has 'condition'
  design    = ~ condition
)

# 2) light filter: drop genes with near-zero counts across all samples
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 3) fit the DESeq2 model (estimates size factors, dispersions, etc.)
dds <- DESeq(dds)

# 4) make transformed matrix for EDA (variance stabilizing transform)
vsd <- vst(dds, blind = TRUE)

# quick sanity checks (should run without error)
dim(dds)
assay(vsd)[1:3, 1:3]  # peek at transformed values

# ---------------------
# step 9: Exploratory data analysis (EDA)
# ---------------------
#library size barplot
# 1) calculate library sizes (total counts per sample)
lib_sizes <- colSums(counts)

# 2)barplot
barplot(lib_sizes,
        main = "Library sizes per sample",
        ylab = "Total counts (reads)",
        las = 2,            # make sample names vertical
        cex.names = 0.7)    # shrink labels so they fit
#what to look for: all samples should have similar bar length (no. of reads)
#export as png
png("results/DEGs/plots/LibrarySizes_barplot.png", width=1200, height=800, res=150)
barplot(colSums(counts), main="Library sizes per sample", ylab="Total reads", las=2, cex.names=0.6)
dev.off()

# PCA by your grouping variable (condition = TF KO)
p <- plotPCA(vsd, intgroup = "condition")

# show it
p
# what to look for: each KO group should cluster and the 3 KO groups should be separated
# export as png
png("results/DEGs/plots/PCA_plot.png", width=1200, height=1000, res=150)
plotPCA(vsd, intgroup="condition")
dev.off()

# sample distance heatmap

library(pheatmap)
library(RColorBrewer)
#to export as png
png("results/DEGs/plots/SampleDistance_heatmap.png", width = 1400, height = 1100, res = 150)

# 1) calculate distances 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# 2) make short labels: condition + replicate index
short_labels <- paste0(coldata$condition, "_rep", seq_len(nrow(coldata)))
rownames(sampleDistMatrix) <- short_labels
colnames(sampleDistMatrix) <- short_labels

# 3) plot with cleaner labels
heat_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = heat_colors,
         main = "Sample-to-sample distances (short labels)")

dev.off()
# what to look for: replicates within the same group are similar, 
#high difference across groups, replicates from the same group 
#cluster together, and replicates from the same KO group cluster tightly

# ------------------------------------
# step 10: Differential expression analysis
# ------------------------------------

#we already created the DESseq object before,
#so will now extract the DE results table for 3 comparisons

# 1- HES5 vs POU3F2

# build the results table (Wald test, FDR at 0.05)
res_HES5_vs_POU3F2 <- results(dds, contrast = c("condition",
                            "HES5","POU3F2"), alpha = 0.05)
# make a tidy data.frame ordered by FDR (padj)
res_HES5_vs_POU3F2_tbl <- as.data.frame(res_HES5_vs_POU3F2)
res_HES5_vs_POU3F2_tbl$gene <- rownames(res_HES5_vs_POU3F2_tbl)
res_HES5_vs_POU3F2_tbl <- res_HES5_vs_POU3F2_tbl[order(res_HES5_vs_POU3F2_tbl$padj), ]

# quick peek
summary(res_HES5_vs_POU3F2)

#result: 86 genes with padj < 0.05 and positive LFC → up in HES5 relative to POU3F2.
#80 genes with padj < 0.05 and negative LFC → up in POU3F2 relative to HES5.

# 2- HES5 VS SOX1

res_HES5_vs_SOX1 <- results(dds, contrast = c("condition",
                                  "HES5","SOX1"), alpha = 0.05)
res_HES5_vs_SOX1_tbl <- as.data.frame(res_HES5_vs_SOX1)
res_HES5_vs_SOX1_tbl$gene <- rownames(res_HES5_vs_SOX1_tbl)
res_HES5_vs_SOX1_tbl <- res_HES5_vs_SOX1_tbl[order(res_HES5_vs_SOX1_tbl$padj), ]
summary(res_HES5_vs_SOX1)

#result: 1130 genes upregulated in HES5, 934 genes upregulated in SOX1

# 3- POU3F2 VS SOX1

res_POU3F2_vs_SOX1 <- results(dds, contrast = c("condition",
                                "POU3F2","SOX1"), alpha = 0.05)
res_POU3F2_vs_SOX1_tbl <- as.data.frame(res_POU3F2_vs_SOX1)
res_POU3F2_vs_SOX1_tbl$gene <- rownames(res_POU3F2_vs_SOX1_tbl)
res_POU3F2_vs_SOX1_tbl <- res_POU3F2_vs_SOX1_tbl[order(res_POU3F2_vs_SOX1_tbl$padj), ]
summary(res_POU3F2_vs_SOX1)

#result: 2041 genes upregulated in POU3F2 and 2061 genes upregulated in SOX1
# total 4132 DE genes, these 2 KOs are the most different in their expression profile

#key finding: SOX1 KO seems to drive the largest expression differences
# compared to the other TFs, HES5 and POU3F2 are relatively similar, 
# maybe SOX1 KO causes a stronger developmental shift at Day 11 compared to the others

# ------------------------------------
# step 11: Plotting results
# ------------------------------------

# MA plot
# MA plot for HES5 vs POU3F2

png("results/DEGs/plots/MAplot_HES5_vs_POU3F2.png", width = 1200, height = 1000, res = 150)
plotMA(res_HES5_vs_POU3F2, ylim = c(-5, 5), 
       main = "MA Plot: HES5 vs POU3F2")
dev.off()

# MA plot for HES5 vs SOX1

png("results/DEGs/plots/MAplot_HES5_vs_SOX1.png", width = 1200, height = 1000, res = 150)
plotMA(res_HES5_vs_SOX1, ylim = c(-5, 5), 
       main = "MA Plot: HES5 vs SOX1")
dev.off()

# MA plot for POU3F2 vs SOX1

png("results/DEGs/plots/MAplot_POU3F2_vs_SOX1.png", width = 1200, height = 1000, res = 150)
plotMA(res_POU3F2_vs_SOX1, ylim = c(-5, 5), 
       main = "MA Plot: POU3F2 vs SOX1")
dev.off()

#Key findings: MA plots show that SOX1 KO samples are most distinct (huge number of DE genes).
#HES5 vs POU3F2 is the most similar comparison (few DE genes).

# Top genes heatmap

BiocManager::install("genefilter")
library(pheatmap)
library(genefilter)
library(SummarizedExperiment)

png("results/DEGs/plots/Top20Genes_heatmap.png", width = 1400, height = 900, res = 150)

## 1) pick top 20 most-variable genes (VST data)
mat_vsd <- assay(vsd)                           # VST matrix: genes x samples
topVarGenes <- head(order(rowVars(mat_vsd), decreasing = TRUE), 20)

## 2) extract those genes and center each gene (subtract its mean) across samples
mat <- mat_vsd[topVarGenes, , drop = FALSE]
mat <- mat - rowMeans(mat)

## 3) build the column annotation 
## use 'condition' (HES5/POU3F2/SOX1).
anno <- data.frame(condition = coldata$condition)
rownames(anno) <- colnames(mat)       # IMPORTANT: rownames must match columns

## short column labels so it’s readable
cond <- coldata[colnames(mat), "condition", drop = TRUE]
rep_idx <- ave(seq_along(cond), cond, FUN = seq_along)
short <- paste0(cond, "_rep", rep_idx)

# define annotation colors for conditions
# Annotation colors: purples/blues
anno_colors <- list(
  condition = c(
    HES5   = "#6A3D9A",  # dark purple
    POU3F2 = "#1F78B4",  # royal blue
    SOX1   = "#A6CEE3"   # light blue
  )
)
## 4) plot heatmap (gene clustering + sample clustering + annotation)
pheatmap(mat,
         annotation_col = anno,
         labels_col     = short,                       # short sample labels
         show_rownames  = TRUE,                        # show the 20 gene names
         main           = "Top 20 most-variable genes (VST, centered)"
         ,color= colorRampPalette(rev(brewer.pal(n = 8, "RdBu")))(255),
         annotation_colors= anno_colors)
dev.off()

#key findings:The top 20 genes cluster samples primarily by condition, 
#with SOX1 KO forming a distinct cluster from HES5 and POU3F2 KOs. 
#This agrees with our DESeq2 results showing SOX1 KO has the largest number of DE genes. 
#HES5 and POU3F2 samples appear more similar to each other, as reflected by their 
#closer clustering and fewer distinct gene blocks.”

#---------------------------------
#step 12: Annotate results
#---------------------------------

#1- load annotation packages


library(AnnotationDbi)
library(org.Hs.eg.db)

##2- annotate HES5 vs POU3F2 results

# 0) pick the results object to annotate (can use res_HES5_vs_POU3F2 OR its *_tbl)
res_obj <- res_HES5_vs_POU3F2   # DESeqResults

# 1) get Ensembl IDs (strip version suffix like ".1" if present)
ens <- sub("\\.\\d+$", "", rownames(res_obj))

# 2) map to SYMBOL and full GENENAME
sym  <- AnnotationDbi::mapIds(org.Hs.eg.db,
                              keys = ens,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")

gname <- AnnotationDbi::mapIds(org.Hs.eg.db,
                               keys = ens,
                               column = "GENENAME",
                               keytype = "ENSEMBL",
                               multiVals = "first")

# 3) build an annotated data.frame, order by FDR
res_annot_HES5_vs_POU3F2 <- as.data.frame(res_obj)
res_annot_HES5_vs_POU3F2$ENSEMBL  <- ens
res_annot_HES5_vs_POU3F2$SYMBOL   <- unname(sym[ens])
res_annot_HES5_vs_POU3F2$GENENAME <- unname(gname[ens])

# put nice columns first and sort by padj (FDR)
res_annot_HES5_vs_POU3F2 <- res_annot_HES5_vs_POU3F2[order(res_annot_HES5_vs_POU3F2$padj), ]
res_annot_HES5_vs_POU3F2 <- res_annot_HES5_vs_POU3F2[, c("ENSEMBL","SYMBOL","GENENAME",
                                                         "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

# quick peek at the top rows
head(na.omit(res_annot_HES5_vs_POU3F2), 10)

## HES5 vs SOX1
ens <- sub("\\.\\d+$", "", rownames(res_HES5_vs_SOX1))
sym  <- mapIds(org.Hs.eg.db, keys=ens, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
gname<- mapIds(org.Hs.eg.db, keys=ens, column="GENENAME", keytype="ENSEMBL", multiVals="first")

res_annot_HES5_vs_SOX1 <- as.data.frame(res_HES5_vs_SOX1)
res_annot_HES5_vs_SOX1$ENSEMBL  <- ens
res_annot_HES5_vs_SOX1$SYMBOL   <- unname(sym[ens])
res_annot_HES5_vs_SOX1$GENENAME <- unname(gname[ens])
res_annot_HES5_vs_SOX1 <- res_annot_HES5_vs_SOX1[order(res_annot_HES5_vs_SOX1$padj), ]
res_annot_HES5_vs_SOX1 <- res_annot_HES5_vs_SOX1[, c("ENSEMBL","SYMBOL","GENENAME",
                                                     "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]

## POU3F2 vs SOX1
ens <- sub("\\.\\d+$", "", rownames(res_POU3F2_vs_SOX1))
sym  <- mapIds(org.Hs.eg.db, keys=ens, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
gname<- mapIds(org.Hs.eg.db, keys=ens, column="GENENAME", keytype="ENSEMBL", multiVals="first")

res_annot_POU3F2_vs_SOX1 <- as.data.frame(res_POU3F2_vs_SOX1)
res_annot_POU3F2_vs_SOX1$ENSEMBL  <- ens
res_annot_POU3F2_vs_SOX1$SYMBOL   <- unname(sym[ens])
res_annot_POU3F2_vs_SOX1$GENENAME <- unname(gname[ens])
res_annot_POU3F2_vs_SOX1 <- res_annot_POU3F2_vs_SOX1[order(res_annot_POU3F2_vs_SOX1$padj), ]
res_annot_POU3F2_vs_SOX1 <- res_annot_POU3F2_vs_SOX1[, c("ENSEMBL","SYMBOL","GENENAME",
                                                         "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
#save the annotated results tabels
write.csv(res_annot_HES5_vs_POU3F2, "results/DEGs/HES5_vs_POU3F2_DEGs.csv", row.names = FALSE)
write.csv(res_annot_HES5_vs_SOX1,   "results/DEGs/HES5_vs_SOX1_DEGs.csv",   row.names = FALSE)
write.csv(res_annot_POU3F2_vs_SOX1, "results/DEGs/POU3F2_vs_SOX1_DEGs.csv", row.names = FALSE)

sessionInfo()

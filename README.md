# Bulk-RNA-seq-
Bulk RNA seq for hipsc cortical brain organoids with 8 different TF KOs. Analyzed the KOs HES5, POU3F2, and SOX1 at day 11 for differential expression and enrichment analysis. To understand how these TF affect cortical brain development.
[README.md](https://github.com/user-attachments/files/22548253/README.md)
**Goal:** Basic downstream RNA-seq analysis to compare KO groups in cortical organoids (Day 11).  
**Data:** GEO GSE303918 (subset); 3 replicates each for HES5, POU3F2, SOX1.

## PART 1: Differential expression analysis
## Pipeline (DESeq2)
1. Load counts & metadata; align sample names.
2. EDA: library sizes, PCA, sample-to-sample distances.
3. DE (design = ~ condition): HES5 vs POU3F2, HES5 vs SOX1, POU3F2 vs SOX1.
4. Annotate results (ENSEMBL → SYMBOL/GENENAME).
5. Save tables and plots.

## Key findings
- **SOX1 KO** is most distinct at Day 11 (largest number of DEGs vs other KOs).
- **HES5 vs POU3F2** shows relatively few DEGs (more similar profiles).
- PCA and distance heatmap: replicates cluster by KO; groups separate cleanly.

## Outputs
- **Results (CSV):** `results/DEGs/*_DEGs.csv` 
- **Plots (PNG):**
  - `PCA_plot.png`
  - `LibrarySizes_barplot.png`
  - `SampleDistance_heatmap.png`
  - `MAplot_*` (one per contrast)
  - `Top20VarGenes_heatmap.png`


--------------------------------------------------
## PART 2: Enrichemnt analysis
Among the pairwise comparisons between TF knockout conditions, the contrast between POU3F2 vs SOX1 showed the most pronounced transcriptional differences (highest DE contrast), with the largest number of differentially expressed genes and strongest fold changes.
As a result, enrichment analysis was focused on this comparison.

## Enrichment Analysis — POU3F2 vs SOX1

This script performs both ORA and GSEA-based GO enrichment analysis on DESeq2 differential expression results for POU3F2 vs SOX1.

### Steps:
1. Load DESeq2 results
2. Run Over-Representation Analysis (ORA) on DE genes (padj < 0.05)
3. save GO Dotplot for top 10 pathways
3. Create ranked gene list (log2FoldChange) for GSEA
4. Retrieve GO:BP gene sets from MSigDB using `msigdbr`
5. Run GSEA using `fgseaMultilevel`
6. Save full enrichment results and generate enrichment plots for top POU3F2- and SOX1-upregulated pathways

### Geenral Key Findings:
- ORA revealed enrichment for ribosome biogenesis, cell cycle, and forebrain development terms.
- GSEA identified pathways upregulated in POU3F2 (e.g., central nervous system neuron differentiation) and in SOX1 (e.g., regulation of cell killing).

Output files saved to `results/enrichment/`

### GO Over-Representation Analysis (ORA)

Over-representation analysis was performed using `clusterProfiler::enrichGO()` on significantly differentially expressed genes (padj < 0.05) from the POU3F2 vs SOX1 comparison.

The top enriched GO Biological Process terms included:

- Ribosome biogenesis
- Translation and rRNA processing
- Cell division and chromatid segregation
- Forebrain and nervous system development

These results suggest that POU3F2-expressing cells show elevated protein synthesis and proliferation, along with a shift toward neurogenic differentiation. This aligns with POU3F2's role as a neurogenic transcription factor, contrasting with SOX1's role in maintaining neural progenitor identity.


### GSEA Analysis (fgsea)

Gene Set Enrichment Analysis (GSEA) of POU3F2 vs SOX1 differential expression results revealed significant enrichment of neurodevelopmental pathways.

The top upregulated pathway in POU3F2 was:

GO: Central nervous system neuron differentiation
(NES = 2.22, padj = 1.66e-09)

This indicates that genes promoting neuronal differentiation are more active in POU3F2-expressing cells, consistent with its role in promoting neurogenesis.
Conversely, pathways upregulated in SOX1 included several cell cycle and stemness-related processes, suggesting maintenance of a progenitor-like state.

The pathway “Regulation of Cell Killing” (GO:BP) was significantly enriched in genes upregulated in SOX1, with genes clustering at the bottom of the ranked list (NES = negative, padj < 0.05).
This suggests that SOX1-expressing progenitor-like cells may exhibit increased activity in genes involved in programmed cell death regulation, potentially reflecting developmental checkpoint mechanisms.

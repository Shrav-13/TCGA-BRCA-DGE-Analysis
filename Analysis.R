# ============================================
# Step 1: Install required packages
# ============================================
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install(c('TCGAbiolinks', 'DESeq2', 'EnhancedVolcano',
                       'clusterProfiler', 'org.Hs.eg.db', 'ComplexHeatmap'))

install.packages(c('ggplot2', 'dplyr', 'pheatmap'))

# Load libraries
library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# ============================================
# Step 2: Query and download TCGA-BRCA RNA-seq data
# ============================================
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)

GDCdownload(query)

data <- GDCprepare(query)

# ============================================
# Step 3: Prepare data and run DESeq2
# ============================================
counts <- assay(data, 'unstranded')
coldata <- colData(data)

# TCGA sample barcodes: 01 = tumour, 11 = normal
coldata$condition <- ifelse(
  substr(coldata$sample_type, 1, 6) == 'Primar',
  'tumour', 'normal'
)

keep <- !is.na(coldata$condition)
counts <- counts[, keep]
coldata <- coldata[keep, ]

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Remove genes with very low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2 (takes 5-20 minutes)
dds <- DESeq(dds)

res <- results(dds,
               contrast = c('condition', 'tumour', 'normal'))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

head(res_df[order(res_df$padj), ])

# ============================================
# Step 4: Visualise results
# ============================================

# Volcano plot
EnhancedVolcano(res_df,
                lab = res_df$gene,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Tumour vs Normal - Differential Gene Expression',
                subtitle = 'TCGA-BRCA',
                legendPosition = 'right'
)

# Heatmap of top 50 DEGs
sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig_genes$gene, ]

mat_scaled <- t(scale(t(mat)))

annotation <- data.frame(
  condition = coldata$condition,
  row.names = colnames(mat)
)

pheatmap(mat_scaled,
         annotation_col = annotation,
         show_rownames = TRUE, show_colnames = FALSE,
         main = 'Top 50 DEGs: Tumour vs Normal'
)

# ============================================
# Step 5: GO enrichment analysis
# ============================================
up_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>%
  pull(gene)

# Strip Ensembl version suffix (e.g. .10) before mapping
up_genes_clean <- sub("\\..*", "", up_genes)

gene_ids <- bitr(up_genes_clean,
                 fromType = 'ENSEMBL',
                 toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db)

ego <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego, showCategory = 20,
        title = 'GO Biological Process Enrichment')
# TCGA-BRCA-DGE-Analysis
Differential gene expression analysis of breast tumour vs normal tissue using TCGA RNA-seq data (DESeq2, clusterProfiler). Identifies significantly dysregulated genes and their associated biological processes

## Results

Of the genes tested, 11,464 were significantly differentially expressed (padj < 0.05, |log2FC| > 1): 7,522 upregulated and 3,942 downregulated in tumour tissue relative to normal. This is a large number relative to typical DGE studies, which reflects the statistical power of the TCGA-BRCA cohort's sample size rather than an unusually strong biological effect on its own.

### Upregulated genes

GO Biological Process enrichment of the upregulated gene set was dominated by cell division and chromatin-related processes, including:

| GO term | Description | Adjusted p-value | Gene count |
|---|---|---|---|
| GO:0006334 | Nucleosome assembly | 4.84e-13 | 51 |
| GO:0061644 | Protein localization to CENP-A containing chromatin | 4.84e-13 | 18 |
| GO:0034728 | Nucleosome organization | 3.12e-12 | 54 |
| GO:0000280 | Nuclear division | 3.59e-12 | 117 |
| GO:0071459 | Protein localization to chromosome, centromeric region | 3.59e-12 | 26 |
| GO:0048285 | Organelle fission | 1.29e-09 | 119 |
| GO:0098813 | Nuclear chromosome segregation | 1.30e-09 | 89 |
| GO:0051321 | Meiotic cell cycle | 1.30e-09 | 82 |

This is consistent with the high proliferation rate of breast tumour tissue: genes controlling chromosome segregation, nucleosome assembly and cell division are upregulated as tumour cells divide more frequently than normal breast tissue.

### Downregulated genes

GO Biological Process enrichment of the downregulated gene set was dominated by muscle, vascular and circulatory processes:

| GO term | Description | Adjusted p-value | Gene count |
|---|---|---|---|
| GO:0008015 | Blood circulation | 1.75e-23 | 123 |
| GO:0003012 | Muscle system process | 3.42e-14 | 104 |
| GO:0006936 | Muscle contraction | 6.66e-14 | 88 |
| GO:0060047 | Heart contraction | 6.84e-13 | 64 |
| GO:0003015 | Heart process | 8.41e-13 | 65 |
| GO:1903522 | Regulation of blood circulation | 2.98e-12 | 66 |
| GO:0003018 | Vascular process in circulatory system | 4.91e-12 | 65 |
| GO:0006935 | Chemotaxis | 4.91e-12 | 93 |

Rather than reflecting an active tumour suppressive process, this most likely reflects a difference in tissue composition. Normal breast tissue contains a mix of epithelial, vascular, smooth muscle and stromal cells, whereas tumour tissue is dominated by proliferating epithelial tumour cells. Genes associated with muscle and vasculature therefore appear depleted in bulk RNA-seq of tumour samples partly because that tissue type is underrepresented, not necessarily because those pathways are being suppressed within individual cells. This is a known limitation of bulk (as opposed to single-cell) RNA-seq differential expression analysis.
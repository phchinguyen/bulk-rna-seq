library(SummarizedExperiment)
# Sample metadata
samples_metadata <- data.frame(sample_id = colnames(count_matrix_rounded),
                            condition = c("normal", "normal", "normal", "tumor", "tumor", "tumor"),
                            row.names = colnames(count_matrix_rounded)
                            )
# Features metadata
gene_metadata <- data.frame(gene_id = rownames(count_matrix_rounded),
                            row.names = rownames(count_matrix_rounded)
                            )
# Summarized Experiment
se <- SummarizedExperiment(
  assays = list(counts = count_matrix_rounded),
  colData = samples_metadata,
  rowData = gene_metadata
)

# Remove low expressed genes
se <- se[rowSums(assay(se, "counts")) >= 10]

library(DESeq2)
dds <- DESeq2::DESeqDataSet(se, design = ~ condition)
dds <- estimateSizeFactors(dds)

library(vsn)
vsd <- vst(dds, blind = TRUE)
meanSdPlot(assay(vsd))

library(ggplot2)
pcaData <- DESeq2::plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = condition, shape = condition), size = 5) + 
  theme_minimal() + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  labs(title = "PCA: Tumor vs Normal Samples") 
  




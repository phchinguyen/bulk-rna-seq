library(SummarizedExperiment)
# Sample metadata
samples_metadata <- data.frame(sample_id = colnames(count_matrix),
                            condition = c("normal", "normal", "normal", "tumor", "tumor", "tumor"),
                            row.names = colnames(count_matrix)
                            )
# Features metadata
gene_metadata <- data.frame(gene_id = rownames(count_matrix),
                            row.names = rownames(count_matrix)
                            )
# Summarized Experiment
se <- SummarizedExperiment(
  assays = list(counts = count_matrix),
  colData = samples_metadata,
  rowData = gene_metadata
)

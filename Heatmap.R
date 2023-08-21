library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(latex2exp)

setwd("Z:/UDP- Research/RNAseq data/DE")

# Import all .csv files in folder
participantfiles <- list.files(pattern = ".*?\\.csv")

# Create empty df with 14150 rows == # of genes
FCs <- as.data.frame(matrix(nrow = 14150, ncol = length(participantfiles)))
genes <- c()

# Import gene reads for all participants, keep vector of significant ones
for (i in 1:length(participantfiles)) {
    participant <- strsplit(participantfiles[i], "[.]")
    participant <- participant[[1]][1]
    import <- read_csv(participantfiles[i]) %>%
        dplyr::select(c(log2FoldChange, hgnc_symbol, sig)) %>%
        arrange(hgnc_symbol)
    FCs[,i] <- import$log2FoldChange
    colnames(FCs)[i] <- participant
    genes <- append(genes,
                    import$hgnc_symbol[which(import$sig == 1)]) %>%
        unique()
}
FCs <- cbind(FCs, import$hgnc_symbol) %>%
    rename(Gene = `import$hgnc_symbol`)

# Merge significant genes with fold changes
FC_sig <- as.data.frame(genes) %>%
    rename(Gene = genes) %>%
    merge(., FCs) %>%
    na.omit() %>%
    column_to_rownames(var = "Gene")

# Melt dataframe for heatmap
forHeatmap <- FC_sig %>%
    rownames_to_column() %>%
    reshape::melt()
colnames(forHeatmap) <- c("Gene", "Patient", "Log2FC")
forHeatmap <- forHeatmap %>%
    arrange(Gene)

# Plot heatmap
ggplot(forHeatmap, aes(Patient, Gene, fill = Log2FC)) +
    geom_tile() +
    scale_fill_gradient2(low = "#075AFF",
                        mid = "#FFFFCC",
                        high = "#FF0000") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 8),
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)) +
    scale_y_discrete(limits = rev) +
    labs(fill = TeX("$Log_{2} Change$"))
ggsave("Heatmap.png", bg = "white")

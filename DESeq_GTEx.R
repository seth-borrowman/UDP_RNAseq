library(tximeta)
library(DESeq2)
library(tidyverse)
library(ashr)
library(biomaRt)
library(ggrepel)

# Set wd
setwd("Z:/UDP- Research/RNAseq data")


# Load data
coldata <- as.data.frame(matrix(nrow = 6, ncol = 2))
colnames(coldata) <- c("names", "files")
coldata[,1] <- c("714051", "715986", "715991", "731949", "731952", "731956")
for (i in seq_along(1:length(coldata[,1]))) {
    coldata[i, 2] <- file.path("quants",
                               paste(coldata[i, 1], "_quant", sep = ""),
                               "quant.sf")
}

case <- read_delim("quants/UIRDB20230008/quant.sf")
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
geneNames <- getBM(filters = "ensembl_transcript_id_version",
                   attributes = c("ensembl_gene_id_version",
                                  "ensembl_transcript_id_version",
                                  "hgnc_symbol"),
                   values = case$Name,
                   mart = mart)
case <- merge(case, geneNames,
              by.x = "Name", by.y = "ensembl_transcript_id_version")

case_sum <- case %>%
    group_by(ensembl_gene_id_version) %>%
    summarise(Reads = sum(NumReads))

# Import GTEx controls
gtex <- read_delim("Readsgtex.gct",
                   delim = '\t', col_names = T, skip = 2) %>%
    rename(gene_name = Name) %>%
    dplyr::select(-c(id, Description))

# Merge with case
gtex_case <- merge(gtex,
                   case_sum %>% dplyr::select(Reads, ensembl_gene_id_version),
                   by.x = "gene_name", by.y = "ensembl_gene_id_version")
row.names(gtex_case) <- gtex_case[,1]
gtex_case <- gtex_case[,-1]



condition <- rep("Control", ncol(gtex_case))
condition[which(colnames(gtex_case) == "Reads")] <- "Case"
ddsTxi <- DESeqDataSetFromMatrix(round(gtex_case), DataFrame(condition),
                                 ~ condition, tidy = F)

dds <- DESeq(ddsTxi)
resLFC <- lfcShrink(dds, coef = "condition_Control_vs_Case", type = "ashr")
results <- results(dds)
resultsDF <- as.data.frame(resLFC)
resultsDF$GeneID <- resLFC@rownames


resultsDF <- left_join(resultsDF,
                       geneNames %>%
                           dplyr::select(c("ensembl_gene_id_version", "hgnc_symbol")) %>%
                           unique(),
                   by = c("GeneID" = "ensembl_gene_id_version"))

resultsDF <- resultsDF %>%
    mutate(sig = case_when(
        -log(pvalue) >= -log(0.05/nrow(resultsDF)) & abs(log2FoldChange) >= 2 ~ 1,
        .default = 0
    ))
resultsDF$sig <- factor(resultsDF$sig, ordered = is.ordered(c(0, 1)))
ggplot(data = resultsDF, aes(x = log2FoldChange, y = -log(pvalue), color = sig)) +
    geom_point(show.legend = F) +
    geom_vline(xintercept = c(-2, 2), color = "pink") +
    geom_hline(yintercept = -log(0.05/nrow(resultsDF)), color = "lightblue") +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    geom_label_repel(data = resultsDF %>%
                        filter(sig == 1),
        aes(label = hgnc_symbol)) +
    theme(legend.position = "none")

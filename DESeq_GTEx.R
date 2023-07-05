library(tximeta)
library(DESeq2)
library(tidyverse)
library(apeglm)
library(biomaRt)

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

# Import GTEx controls
gtex <- read_delim("Readsgtex.gct",
                   delim = '\t', col_names = T, skip = 2) %>%
    rename(gene_name = Name) %>%
    dplyr::select(-id)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_to_transcript <- getBM(attributes = c("ensembl_gene_id_version",
                                           "ensembl_transcript_id_version"),
                            mart = mart)
colnames(gene_to_transcript) <- c("gene_name", "trans_name")

gtex <- gtex %>%
    left_join(gene_to_transcript, ., by = "gene_name")
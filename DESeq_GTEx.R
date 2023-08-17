library(tximeta)
library(DESeq2)
library(tidyverse)
library(apeglm)
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

case <- read_delim("quants/731952_quant/quant.sf")

# Import GTEx controls
gtex <- read_delim("Readsgtex.gct",
                   delim = '\t', col_names = T, skip = 2) %>%
    rename(gene_name = Name) %>%
    dplyr::select(-id)

## Get gene names
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
geneNames <- getBM(filters = "ensembl_gene_id_version",
                   attributes = c("ensembl_gene_id_version",
                                  "ensembl_transcript_id_version",
                                 "hgnc_symbol"),
                   values = gtex$gene_name,
                   mart = mart)

## Merge case and controls
case_counts <- case %>% select(Name, NumReads)
colnames(case_counts) <- c("Gene", "Case")
case_counts <- merge(case_counts, geneNames,
                     by.x = "Gene", by.y = "ensembl_transcript_id_version")
gtex_case <- merge(gtex, case_counts,
                   by.x = "gene_name", by.y = "ensembl_gene_id_version") %>%
    select(-c(hgnc_symbol, Description, gene_name))

gtex1 <- gtex_case %>%
    select(-Gene) %>%
    t() %>%
    as.data.frame
colnames(gtex1) <- make.names(gtex_case$Gene, unique = T)
gtex1$CaseControl <- c(rep("Control", nrow(gtex1)))
gtex1$CaseControl[which(row.names(gtex1) == "Case")] <- "Case"

## DE analysis
DE_results <- data.frame(GeneID = colnames(gtex1[,1:(ncol(gtex1) - 1)]),
                         Mean = rep(NA, (ncol(gtex1) - 1)),
                         SE = rep(NA, (ncol(gtex1) - 1)),
                         FoldChange = rep(NA, (ncol(gtex1) - 1)),
                         pval = rep(NA, (ncol(gtex1) - 1)))


DE_results[,c(2,4)] <- gtex1 %>%
    group_by(CaseControl) %>%
    summarise(across(everything(), c(mean))) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    janitor::row_to_names(row_number = 1) %>% 
    apply(., 2, as.numeric) %>%
    as.data.frame() %>%
    mutate(FC = case_when(
        is.na(Case) | is.na(Control) ~ NA,
        Control != 0 & Case != 0 ~ (Case - Control) / Control,
        Case == 0 | Control == 0 ~ 0)) %>%
    select(c(Control, FC))

DE_results[,3] <- gtex1 %>%
    group_by(CaseControl) %>%
    summarise(across(everything(), c(sd))) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    janitor::row_to_names(row_number = 1) %>% 
    apply(., 2, as.numeric) %>%
    as.data.frame() %>%
    select(Control)

for (i in 1:(ncol(gtex1) - 1)) {
    x <- gtex1[which(gtex1$CaseControl == "Case"), i]
    y <- gtex1[which(gtex1$CaseControl == "Control"), i]
    test <- wilcox.test(x, y)
    DE_results$pval[i] <- test$p.value
}

DE_results <- DE_results %>%
    mutate(Log2FoldChange = case_when(
        FoldChange < 0 ~ -log(abs(FoldChange), base = 2),
        FoldChange > 0 ~ log(FoldChange, base = 2),
        FoldChange == 0 ~ 0
    )) %>%
    merge(., geneNames, by.x = "GeneID", by.y = "ensembl_transcript_id_version")

ggplot(data = DE_results, aes(x = Log2FoldChange, y = -log(pval))) +
    geom_point(show.legend = F) +
    geom_vline(xintercept = c(-2, 2), color = "pink") +
    geom_hline(yintercept = -log(0.05/nrow(DE_results)), color = "lightblue") +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    geom_label_repel(data = DE_results %>% mutate(lab = case_when(
        -log(pval) >= -log(0.05/nrow(DE_results)) & abs(Log2FoldChange) > 2 ~ 1,
        .default = 0
    )) %>%
        filter(lab == 1),
    aes(label = hgnc_symbol)) +
    theme(legend.position = "none")

library(tximeta)
library(DESeq2)
library(tidyverse)
library(ashr)
library(biomaRt)
library(ggrepel)
library(tximport)

# Set wd
setwd("Z:/UDP_Research/RNAseq")
ifelse(!dir.exists('DE'), dir.create('DE'), FALSE)

# Import GTEx controls
gtex <- read_delim("Readsgtex.gct",
                   delim = '\t', col_names = T, skip = 2) %>%
    dplyr::rename(gene_name = Name) %>%
    dplyr::select(-c(id, Description))

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


# Load data
participants <- list.files("quants/")

for (i in participants) {
    
    # Load case data
    case <- tximport(paste("quants/", i, "/quant.sf", sep = ""),
                     type = "salmon",
                     txOut = T) %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(ensembl_transcript_id_version = rowname, NumReads = counts)
    
    geneNames <- getBM(filters = "ensembl_transcript_id_version",
                       attributes = c("ensembl_transcript_id_version",
                                      "ensembl_gene_id_version",
                                      "hgnc_symbol"),
                       values = case$ensembl_transcript_id_version,
                       mart = mart)
    case <- merge(case, geneNames)
    
    # IIHG gives results by transcript ID, GTEx is by gene - convert to gene
    case_sum <- case %>%
        group_by(ensembl_gene_id_version) %>%
        summarise(Reads = sum(NumReads))


    # Merge with controls
    gtex_case <- merge(gtex, case_sum,
                       by.x = "gene_name", by.y = "ensembl_gene_id_version")
    row.names(gtex_case) <- gtex_case[,1] # Make genes column name
    gtex_case <- gtex_case[,-1]

    # Set case/control condition
    condition <- rep("Control", ncol(gtex_case))
    condition[which(colnames(gtex_case) == "Reads")] <- "Case"

    # Create DSEq dataset
    ddsTxi <- DESeqDataSetFromMatrix(round(gtex_case), DataFrame(condition),
                                 ~ condition, tidy = F)

    # DESeq analysis with ASHR shrinkage
    dds <- DESeq(ddsTxi)
    resLFC <- lfcShrink(dds, coef = "condition_Control_vs_Case", type = "ashr")
    resultsDF <- as.data.frame(resLFC)
    resultsDF$GeneID <- resLFC@rownames

    # Move results to df with gene symbol
    resultsDF <- left_join(resultsDF,
                       geneNames %>%
                           dplyr::select(c("ensembl_gene_id_version",
                                           "hgnc_symbol")) %>%
                           unique(),
                   by = c("GeneID" = "ensembl_gene_id_version"))

    # Find significant results - Bonferroni p, arbitrary log2FC >= 2
    resultsDF <- resultsDF %>%
        mutate(sig = case_when(
            -log(pvalue) >= -log(0.05/nrow(resultsDF)) &
                abs(log2FoldChange) >= 2 ~ 1,
            .default = 0
        ))
    resultsDF$sig <- factor(resultsDF$sig, ordered = is.ordered(c(0, 1)))
    write_csv(resultsDF, paste("DE/", i, ".csv", sep = ""))

    # Plot
    options(ggrepel.max.overlaps = Inf)
    ggplot(data = resultsDF, aes(x = log2FoldChange, y = -log(pvalue), color = sig)) +
        geom_point(show.legend = F) +
        #geom_vline(xintercept = c(-2, 2), color = "pink") +
        geom_hline(yintercept = -log(0.05/nrow(resultsDF)), color = "lightblue") +
        scale_color_manual(values = c("grey", "red")) +
        theme_minimal() +
        geom_text_repel(data = resultsDF %>%
                            filter(sig == 1),
            aes(label = hgnc_symbol)) +
        theme(legend.position = "none") +
        ggtitle(i)

    ggsave(paste("DE/", i, "_volcano.png", sep = ""), bg = "white")
}

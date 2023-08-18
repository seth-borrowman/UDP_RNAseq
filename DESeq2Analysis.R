library(tximeta)
library(DESeq2)
library(apeglm)
library(vsn)
library(biomaRt)
library(tidyverse)
library(ggrepel)

# Set wd
setwd("Z:/UDP- Research/RNAseq data")

# Load data
coldata <- as.data.frame(matrix(nrow = 12, ncol = 2))
colnames(coldata) <- c("names", "files")
coldata[,1] <- c("UIRDB20230008",
                 "UIRDB20230003", "UIRDB20230005", "UIRDB20230006",
                 "UIRDB20230007", "UIRDB20230009", "UIRDB20230009-A",
                 "UIRDB20230009-B", 'UIRDB20230013', "UIRDB20230015",
                 "UIRDB20230016","UIRDB20230018")
for (i in seq_along(1:length(coldata[,1]))) {
    coldata[i, 2] <- file.path("quants",
                               paste(coldata[i, 1], "_quant", sep = ""),
                               "quant.sf")
}

### 731952 as outgroup
coldata$treat <- c("case", rep("control", nrow(coldata) - 1))
coldata$treat %>%
    factor() %>%
    relevel(ref = "control")
se <- tximeta(coldata)
ddsTxi <- DESeqDataSet(se, design = ~ treat)

# DE analysis
#keep <- rowSums(counts(ddsTxi)) >= 10 # Filter low count genes
#ddsTxi <- ddsTxi[keep, ]
dds <- DESeq(ddsTxi)

# vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)
# ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))
# meanSdPlot(assay(vsd))
# meanSdPlot(assay(rld))

resLFC <- lfcShrink(dds, coef = "treat_control_vs_case", type = "apeglm")

plotMA(resLFC, ylim = c(-5, 5))
results <- results(dds)
resultsDF <- as.data.frame(results)
resultsDF$TranscriptID <- results@rownames

# Get gene symbols
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_transcript_id","hgnc_symbol"),values=results@rownames,mart= mart)
resultsDF <- merge(resultsDF, G_list, by.x = "TranscriptID", by.y = "ensembl_transcript_id")





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
    geom_label_repel(data = resultsDF %>% mutate(lab = case_when(
        -log(pvalue) >= -log(0.05/nrow(resultsDF)) & log2FoldChange <= -2 ~ 1,
        -log(pvalue) >= 20 & log2FoldChange > 3 ~ 1,
        .default = 0
        )) %>%
        filter(lab == 1),
        aes(label = hgnc_symbol)) +
    theme(legend.position = "none")

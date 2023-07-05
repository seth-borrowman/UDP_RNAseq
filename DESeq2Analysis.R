library(tximeta)
library(DESeq2)
library(magrittr)
library(apeglm)
library(vsn)

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

### 714051 as outgroup
coldata51 <- coldata
coldata51$treat <- c("case", rep("control", 5))
coldata51$treat %>%
    factor() %>%
    relevel(ref = "control")
se <- tximeta(coldata51)
ddsTxi <- DESeqDataSet(se, design = ~ treat)

# DE analysis
#keep <- rowSums(counts(ddsTxi)) >= 10 # Filter low count genes
#ddsTxi <- ddsTxi[keep, ]
dds <- DESeq(ddsTxi)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

resLFC <- lfcShrink(dds, coef = "treat_control_vs_case", type = "apeglm")

plotMA(resLFC, ylim = c(-5, 5))
results <- results(dds)
results

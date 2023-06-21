library(tximeta)
library(DESeq2)

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

# controls <- as.data.frame(matrix(nrow = 4, ncol = 2))
# colnames(controls)  <- c("names", "files")
# controls [,1] <- c("ENCFF446MZM", "ENCFF545DTG", "ENCFF676GUE", "ENCFF730OTJ")
# for ( i in seq_along(1:length(controls[,1]))) {
#     controls[i, 2] <- file.path("Controls",
#                                 paste(controls[i, 1], ".tsv", sep = ""))
# }
# coldata <- rbind(coldata, controls)
# coldata$treat <- c(rep("case", 6), rep("control", 4))
# coldata$treat <- factor(coldata$treat)
#cts <- coldata[,1]

se <- tximeta(coldata)
ddsTxi <- DESeqDataSet(se, design = ~ names)

# DE analysis
keep <- rowSums(counts(ddsTxi)) >= 10 # Filter low count genes
ddsTxi <- ddsTxi[keep, ]
dds <- DESeq(ddsTxi)

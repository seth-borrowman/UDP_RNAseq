library(dplyr)
library(biomaRt)
library(tibble)
library(tximeta)
library(naivebayes)
library(DESeq2)
library(ashr)
library(latex2exp)
library(ggplot2)
library(ggfortify)

setwd("Z:/UDP_Research/RNAseq")

# Import GTEx controls
gtex <- readr::read_delim("Readsgtex.gct",
                          delim = '\t', col_names = T, skip = 2) %>%
    dplyr::rename(gene_name = Name) %>%
    dplyr::select(-c(id, Description))

# Use biomart for gene neames
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

geneid <- getBM(filters = "ensembl_gene_id_version",
                attributes = c("ensembl_gene_id_version",
                               "ensembl_gene_id"),
                values = gtex$gene_name,
                mart = mart)
gtex <- merge(gtex, geneid, by.x = "gene_name", by.y = "ensembl_gene_id_version")
rownames(gtex) <- gtex$ensembl_gene_id
gtex <- gtex %>% dplyr::select(-c(gene_name, ensembl_gene_id))

# Load data
participants <- list.files("quants/")

# Load all cases
cases <- matrix()
for (i in participants) {
    
    # Load case data
    case <- tximport::tximport(paste("quants/", i, "/quant.sf", sep = ""),
                               type = "salmon",
                               txOut = T) %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(ensembl_transcript_id_version = rowname, NumReads = counts)
    
    geneNames <- getBM(filters = "ensembl_transcript_id_version",
                       attributes = c("ensembl_transcript_id_version",
                                      "hgnc_symbol", "ensembl_gene_id"),
                       values = case$ensembl_transcript_id_version,
                       mart = mart)
    case <- merge(case, geneNames)
    
    # IIHG gives results by transcript ID, GTEx is by gene - convert to gene
    case_sum <- case %>%
        group_by(ensembl_gene_id) %>%
        summarise(Reads = sum(NumReads))
    colnames(case_sum)[2] <- i
    
    cases <- cbind(cases, case_sum)
}

# Make sure genes are the same in both sets
rownames(cases) <- cases$ensembl_gene_id
cases <- cases[,which(colnames(cases) %in% participants)]

gtex <- gtex %>%
    filter(rownames(gtex) %in% rownames(cases))
cases <- cases %>%
    filter(rownames(cases) %in% rownames(gtex))

### Naive bayes ----
# Combine cases and controls
combined <- merge(cases, gtex, by = 0)
rownames(combined) <- combined$Row.names
combined <- combined %>%
    dplyr::select(-Row.names) %>%
    t() %>%
    as.data.frame()
# Let's take a look  at PCA
combined <- combined[,which(colSums(combined) > 0)]
pca <- prcomp(combined, scale. = T)
combined1 <- combined
combined1$case <- c(rep("UIRDB", 15), rep("GTEx", 755))
combined1$case <- factor(combined1$case, levels = c("GTEx", "UIRDB"))
autoplot(pca, data = combined1, color = "case") +
    theme_minimal()

# Find parameters that are very different from both UIRDB and GTEx
# Use Naive Bayes parameter estimates
parameters <- data.frame(row.names = colnames(combined))
for (i in 1:15) {
    case <- c(rep("UIRDB", 15), rep("GTEx", 755))
    case[i] <- "Case"
    combined$case <- case
    combined$case <- factor(combined$case, levels = c("Case", "UIRDB", "GTEx"))
    mod <- multinomial_naive_bayes(x = dplyr::select(combined, -case),
                                   y = combined$case,
                                   laplace = 1)
    parameters <- cbind(parameters, sqrt((coef(mod)[,1] - coef(mod)[,2])^2 + (coef(mod)[,1] - coef(mod)[,3])^2))
    colnames(parameters)[i] <- rownames(combined)[i]
}


zadj <- (parameters - rowMeans(parameters)) / apply(parameters, 1, sd)
zadj$max <- apply(zadj, 1, function(x) max(abs(x)))
selectedz <- zadj[which(zadj$max > zadj$max[order(zadj$max, decreasing = T)[21]]),]#2.5758),]

### Get differential expression ----
# Import all .csv files in folder
participantfiles <- list.files(path = "Z:/UDP_Research/RNAseq/DE/",
                               pattern = ".*?\\.csv")

# Create empty df with 14150 rows == # of genes
FCs <- as.data.frame(matrix(nrow = 14150, ncol = length(participantfiles)))
genes <- c()

# Import gene reads for all participants, keep vector of significant ones
for (i in 1:length(participantfiles)) {
    participant <- strsplit(participantfiles[i], "[.]")
    participant <- participant[[1]][1]
    import <- readr::read_csv(paste("Z:/UDP_Research/RNAseq/DE/",
                                    participantfiles[i], sep = "")) %>%
        dplyr::select(c(log2FoldChange, GeneID, sig)) %>%
        arrange(GeneID)
    FCs[,i] <- import$log2FoldChange
    colnames(FCs)[i] <- participant
    genes <- append(genes,
                    import$GeneID[which(import$sig == 1)]) %>%
        unique()
}
FCs <- cbind(FCs, import$GeneID)
colnames(FCs)[16] <- "ensembl_gene_id_version"

geneNames <- getBM(filters = "ensembl_gene_id_version",
                   attributes = c("ensembl_gene_id_version",
                                  "hgnc_symbol", "ensembl_gene_id"),
                   values = FCs$ensembl_gene_id_version,
                   mart = mart)
FCs <- merge(FCs, geneNames)
FCs$hgnc_symbol <- ifelse(FCs$hgnc_symbol == "", FCs$ensembl_gene_id,
                          FCs$hgnc_symbol)

FC_selected <- FCs[which(FCs$ensembl_gene_id %in% rownames(selectedz)),] %>%
    dplyr::select(-c(ensembl_gene_id, ensembl_gene_id_version))

# Melt dataframe for heatmap
forHeatmap <- FC_selected %>%
    rownames_to_column() %>%
    reshape::melt() %>%
    dplyr::select(-rowname)
colnames(forHeatmap) <- c("Gene", "Patient", "Log2FC")
forHeatmap <- forHeatmap %>%
    arrange(Gene)

# Rename those left unnamed
forHeatmap[which(forHeatmap$Gene == "ENSG00000224967"), 1] <- "AC009303.3"
forHeatmap[which(forHeatmap$Gene == "ENSG00000225854"), 1] <- "RP11-569G9.7"
forHeatmap[which(forHeatmap$Gene == "ENSG00000226989"), 1] <- "AL049758.2"
forHeatmap[which(forHeatmap$Gene == "ENSG00000233287"), 1] <- "AC009362.2"
forHeatmap[which(forHeatmap$Gene == "ENSG00000234648"), 1] <- "AL162151.3"
forHeatmap[which(forHeatmap$Gene == "ENSG00000237176"), 1] <- "RP11-813P10.2"
forHeatmap[which(forHeatmap$Gene == "ENSG00000238084"), 1] <- "RP3-469D22.1"
forHeatmap[which(forHeatmap$Gene == "ENSG00000254719"), 1] <- "RP11-351I24.3"
forHeatmap[which(forHeatmap$Gene == "ENSG00000266356"), 1] <- "RP11-578C11.2"

forHeatmap[which(forHeatmap$Gene == "ENSG00000231392"), 1] <- "ENSG00000231392"
forHeatmap[which(forHeatmap$Gene == "ENSG00000258111"), 1] <- "ATP5G1"
forHeatmap[which(forHeatmap$Gene == "ENSG00000260335"), 1] <- "CORO1A"
forHeatmap[which(forHeatmap$Gene == "ENSG00000260882"), 1] <- "ENSG00000260882"
# Remove GAGE2E
forHeatmap <- forHeatmap[-which(forHeatmap$Gene == "GAGE2E" | forHeatmap$Gene == "PHB1P2"),]

forHeatmap <- forHeatmap %>%
    arrange(Gene)

# Plot heatmap
ggplot(forHeatmap, aes(Patient, Gene, fill = Log2FC)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue4",
                         mid = "lightyellow1",
                         high = "red3") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 8),
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)) +
    scale_y_discrete(limits = rev) +
    labs(fill = TeX("$Log_{2} Change$")) +
    ggtitle("Log fold change in gene expression vs GTEx controls")
#ggsave("Heatmap.png", bg = "white")

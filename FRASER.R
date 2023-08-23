library(magrittr)
library(FRASER)
library(biomaRt)
setwd("LSS/lss_chndr/UDP_Research/RNAseq/BAMs")

# Folder for FRASER to use
workingDir <- ("FRASER")

# Import all BAM files from folder
bamFiles <- list.files(pattern = "*?\\.bam$")
sampleIDs <- stringr::str_split(bamFiles, pattern = "[.]") %>%
    unlist()
sampleIDs <- sampleIDs[seq(1, length(sampleIDs), 2)]

# Create table of sample files (must be >1 sample)
sampleTable <- data.table::data.table(sampleID = sampleIDs,
                                      bamFile = bamFiles,
                                      group = c(1:length(sampleIDs)),
                                      gene = rep(NA, length(sampleIDs)),
                                      pairedEnd = rep(T, length(sampleIDs)))

settings <- FraserDataSet(colData = sampleTable, workingDir = workingDir)

# Use multiple cores if available
if (.Platform$OS.type == "unix") {
    register(MulticoreParam(workers = min(10, multicoreWorkers())))
} else {
    register(SnowParam(workers = min(10, multicoreWorkers())))
}

# Create count data
fds <- countRNAData(settings)

# Calculate psi values
fds_psi <- calculatePSIValues(fds)

# Filter
# At least one sample has 20 reads
# 5% of the samples have at least 1 read
fds_filtered <- filterExpressionAndVariability(fds_psi,
                                               minDeltaPsi = 0.0,
                                               filter = FALSE)
plotFilterExpression(fds_filtered, bins = 100)
fds_filtered <- fds_filtered[mcols(fds_filtered, type = "j")[,"passed"],]

# Fit splicing model
fds_fraser <- FRASER(fds_filtered, q = c(psi5 = 3, psi3 = 5, theta = 2))
plotCountCorHeatmap(fds_filtered,
                    type = "psi5",
                    logit = TRUE,
                    normalized = FALSE)

# Annotate introns with the HGNC symbols of the corresponding gene
fds <- annotateRanges(fds)
# Retrieve unfiltered results
res <- results(fds_annotated, padjCutoff = NA)
res

plotVolcano(fds_annotated, type = "psi3", "UIRDB20230016")

# Sashimi plot doesn't seem to exist anymore?
#plotBamCoverageFromResultTable(fds_annotated, res[1,], control_samples = "UIRDB20230003")

# Export results
res_df <- as.data.frame(res)
readr::write_csv(res, "FRASER_Results.csv")
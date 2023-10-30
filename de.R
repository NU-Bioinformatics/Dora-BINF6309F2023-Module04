#!usr/bin/env Rscript
# de.R

library(tximport)
library(readr)
library(DESeq2)
library(knitr)

# TODO: update constants for your machine
# Define constants
TESTING <- FALSE # Change to FALSE if using entire Samples set
RESULTS_DIR <- "/home/beslic.i/BINF6309/module04/results"
AIPTASIA_DIR <- "/courses/BINF6309/AiptasiaMiSeq"

# for testing purposes - alternative samples table
testing_samples <- data.frame(Sample = c("Aip02", "Aip02", "Aip02", "Aip02"),
                              Menthol = c("Control", "Control", "Menthol", "Menthol"),
                              Vibrio = c("Control", "Vibrio", "Control", "Vibrio"))
head(testing_samples)

# True script begins
tx2gene <- read.csv(file.path(RESULTS_DIR, "tx2gene.csv"))
head(tx2gene)

if (TESTING) {
  print("***Running test with Aip02 only***")
  samples <- testing_samples
} else {
  samples <- read.csv(file.path(AIPTASIA_DIR, "Samples.csv"), header=TRUE)
}
head(samples)

files <- file.path(RESULTS_DIR, "quant", samples$Sample, "quant.sf")

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = samples, 
                                design = ~ Menthol + Vibrio)

dds$Vibrio <- relevel(dds$Vibrio, ref = "Control")
dds$Menthol <- relevel(dds$Menthol, ref = "Control")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

padj <- .05
minLog2FoldChange <- .5
dfAll <- data.frame()
# Get all DE results except Intercept, and "flatten" into a single file.
for (result in resultsNames(dds)){
  if(result != 'Intercept'){
    res <- results(dds, alpha=.05, name=result)
    dfRes <- as.data.frame(res)
    dfRes <- subset(subset(dfRes, select=c(log2FoldChange, padj)))
    dfRes$Factor <- result
    dfAll <- rbind(dfAll, dfRes)
  }
}
head(dfAll)

write.csv(dfAll, file=file.path(RESULTS_DIR, "dfAll.csv"))

# Read pathway information from the tx2gene.csv file in the RESULTS_DIR directory
pathway_info <- read.csv(file.path(RESULTS_DIR, "tx2gene.csv"))

# Filter for adjusted p-value < 0.05
dfFiltered <- dfAll[dfAll$padj < 0.05, ]
kable(head(dfFiltered))
print(dim(dfFiltered))

# Merge with pathway information
merged_df <- merge(dfFiltered, pathway_info, by.x = 0, by.y = "ko")
names(merged_df)[names(merged_df) == "Row.names"] <- "ko"
kable(head(merged_df))
print(dim(merged_df))

# Write the merged and filtered results to CSV
write.csv(merged_df, file=file.path(RESULTS_DIR, "deAnnotated.csv"))

# end of de.R script

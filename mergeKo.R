#!/usr/bin/env Rscript
# mergeKo.R

#install.packages("knitr")

# Load the kable library to display formatted tables
library(knitr)

# TODO: update constants for your machine

# Define output directory
OUTDIR <- "/home/beslic.i/BINF6309/module04/results"

# Load BLAST results as a table
blastFile <- "/courses/BINF6309/data_BINF6309/Module4/Annotation/transcriptBlast.txt"
keggFile <- "/courses/BINF6309/data_BINF6309/Module4/Annotation/kegg.txt"
koFile <- "/courses/BINF6309/data_BINF6309/Module4/Annotation/ko.txt"

blast <- read.table(blastFile, sep="\t", header=FALSE)
# Set column names to match fields selected in BLAST
colnames(blast) <- c("trans", "sp", "qlen", "slen", "bitscore", 
                     "length", "nident", "pident", "evalue", "ppos")
# Calculate the percentage of identical matches relative to subject length
blast$cov <- blast$nident/blast$slen
# Filter for at least 50% coverage of subject(SwissProt) sequence
blast <- subset(blast, cov > .5)
# Check the blast table
kable(head(blast))

# Load SwissProt to KEGG as a table
kegg <- read.table(keggFile, sep="\t", header=FALSE)
# Set the Swissprot to KEGG column names
colnames(kegg) <- c("sp", "kegg")
# Remove the up: prefix from sp column
kegg$sp <- gsub("up:", "", kegg$sp)
# Check the kegg table
kable(head(kegg))

# Merge BLAST and SwissProt-to-KEGG
blastKegg <- merge(blast, kegg)
# Check the merged table
kable(head(blastKegg))

# Load KEGG to KO as a table
ko <- read.table(koFile, sep="\t", header=FALSE)
# Set column names
colnames(ko) <- c("kegg", "ko")
# Check the ko table
kable(head(ko))

# Merge KOs
blastKo <- merge(blastKegg, ko)
# Check the blast ko table
kable(head(blastKo))

tx2gene <- unique(subset(blastKo, select=c(trans, ko)))
# Check the tx2gene table
kable(head(tx2gene))

# Write as a csv file, excluding row.names
write.csv(tx2gene, file=file.path(OUTDIR, "tx2gene.csv"), row.names=FALSE)
# end of mergeKo.R script


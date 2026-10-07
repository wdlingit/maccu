#!/usr/bin/env Rscript

suppressMessages(library(edgeR))

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: TMM.R <input.tsv> <rawName> <output.tsv>")
}
inFile <- args[1]
rowname <- args[2]
outFile <- args[3]


args = commandArgs(trailingOnly=TRUE)

x <- read.delim(inFile,row.names=rowname)
dge <- DGEList(counts=x)
dge <- calcNormFactors(dge)
v <- voom(dge,normalize="none")
write.table(v$E, file = outFile, sep = "\t", quote = FALSE, col.names = NA)

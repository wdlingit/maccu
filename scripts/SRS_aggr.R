#!/usr/bin/env Rscript

suppressMessages(library(rhdf5))

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript SRS_aggr.R <input.h5> <SRR_SRS> <output.tsv>")
}
h5file <- args[1]
srr_srs <- args[2]
outfile <- args[3]

# Load row/col names
rownames <- h5read(h5file, "/rownames")
colnames <- h5read(h5file, "/colnames")

# Load SRR -> SRS mapping
map <- read.table(srr_srs, header = FALSE, stringsAsFactors = FALSE)
colnames(map) <- c("SRR","SRS")

# Keep only SRRs present in mapping
keep_idx <- which(rownames %in% map$SRR)
sub_mat <- h5read(h5file, "/bigmatrix", index = list(keep_idx, NULL))

# Reorder mapping to match matrix rows
map_sub <- map[match(rownames[keep_idx], map$SRR), ]

# Aggregate counts by SRS (sum across SRRs belonging to same SRS)
agg_mat <- rowsum(sub_mat, group = map_sub$SRS)

# Assign proper colnames (genes) and rownames (SRS)
colnames(agg_mat) <- colnames

# Prepare output file
con <- file(outfile, "w")

# Write header: first column is "gene", then SRS names
writeLines(paste(c("Symbol", rownames(agg_mat)), collapse="\t"), con)

# Stream column by column: each gene becomes one line
for (j in seq_len(ncol(agg_mat))) {
    colvec <- agg_mat[, j]
    line   <- paste(c(colnames(agg_mat)[j], colvec), collapse = "\t")
    writeLines(line, con)
}

close(con)

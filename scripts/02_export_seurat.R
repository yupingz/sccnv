# Export a Seurat object (counts + metadata) for the Python CNV pipeline.
#   Rscript export_seurat.R <object.rds> <out_dir> [integer]
# Writes counts_i.bin / counts_p.bin / counts_x*.bin (CSC pieces, genes x cells), genes.txt, cells.txt, meta.csv.
suppressMessages({library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
o <- readRDS(args[1]); d <- args[2]; integer <- length(args) > 2 && args[3] == "integer"
dir.create(d, recursive = TRUE, showWarnings = FALSE)
cn <- GetAssayData(o, layer = "counts")                       # Seurat v5; use slot = "counts" on v4
wb <- function(v, f) { con <- file(file.path(d, f), "wb"); writeBin(v, con); close(con) }
wb(as.integer(cn@i), "counts_i.bin"); wb(as.integer(cn@p), "counts_p.bin")
if (integer) wb(as.integer(cn@x), "counts_x_int32.bin") else wb(as.numeric(cn@x), "counts_x_f64.bin")
writeLines(rownames(cn), file.path(d, "genes.txt")); writeLines(colnames(cn), file.path(d, "cells.txt"))
md <- o@meta.data; md$cell <- rownames(md); write.csv(md, file.path(d, "meta.csv"), row.names = FALSE)
cat("exported", dim(cn), "to", d, "\n")

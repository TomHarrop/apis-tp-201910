#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(GenomicFeatures)

gff_file <- snakemake@input[["gff"]]
regions_file <- snakemake@output[["regions"]]

# dev
gff_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.gff"

# make a TXDB
txdb <- makeTxDbFromGFF(file = gff_file,
                        format = "gff3",
                        dataSource = "Amel_HAv3.1",
                        organism = "Apis mellifera",
                        taxonomyId = 7460)
# get features
all_exons <- exonsBy(txdb, "gene")

# find the exon
hvr_exon <- all_exons$Csd[mcols(all_exons$Csd)$exon_name == "id75802"]
hvr_dt <- as.data.table(hvr_exon)
outlines <- hvr_dt[, paste0(seqnames, ":", start, "-", end)]

# write the output in samtools format
writeLines(outlines, regions_file, sep = " ")

# log
sessionInfo()

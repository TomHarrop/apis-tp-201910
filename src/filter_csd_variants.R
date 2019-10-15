#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(GenomicFeatures)
library(VariantAnnotation)

gff_file <- snakemake@input[["gff"]]
vcf_file <- snakemake@input[["vcf"]]
tbi_file <- snakemake@input[["tbi"]]
fa_file <- snakemake@input[["fa"]]

csd_vcf <- snakemake@output[["csd"]]
coding_file <- snakemake@output[["coding"]]

# gff_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.gff"
# fa_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
# vcf_file <- "output/030_freebayes/flongle/variants.vcf.gz"
# tbi_file <- "output/030_freebayes/flongle/variants.vcf.gz.tbi"
# csd_vcf <- "test/csd.vcf"

# read the txdb
txdb <- makeTxDbFromGFF(file = gff_file,
                        format = "gff3",
                        dataSource = "Amel_HAv3.1",
                        organism = "Apis mellifera",
                        taxonomyId = 7460)

# load the fasta
fa <- FaFile(fa_file, paste0(fa_file, ".fai"))

# read the vcf
tabix_file <- TabixFile(vcf_file, tbi_file)
vcf <- readVcf(tabix_file)

# get coding regions
coding <- predictCoding(vcf, txdb, fa)

# clean up / release mem
rm(vcf, fa, txdb)
gc(TRUE)

# subset by gene-level coding variants
goi_ids <- "Csd"
rng <- coding[coding$GENEID %in% goi_ids &
                  coding$QUAL > 30 &
                  coding$CONSEQUENCE == "nonsynonymous"]
vcf_rng <- readVcf(tabix_file, param = rng)

# write output
writeVcf(vcf_rng, csd_vcf, index = FALSE)
saveRDS(coding, coding_file)

# log
sessionInfo()

library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)

options(ucscChromosomeNames=FALSE)

gff_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.gff"
txdb_file <- "test/txdb.sqlite"
if(!file.exists(txdb_file)) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(gff_file)
    AnnotationDbi::saveDb(txdb, txdb_file)
} else {
    txdb <- AnnotationDbi::loadDb(txdb_file)
}

# set up transcripts, exons and coding sequences
tx <- transcriptsBy(txdb, "gene")
exons <- exonsBy(txdb, "gene")
cds <- cdsBy(txdb, "gene")

# find Csd in the txdb
csd_transcripts <- tx$Csd
csd_exons <- exons$Csd
csd_cds <- cds$Csd
csd_cds$cds_name <- "complementary sex determiner"
csd_cds$id <- csd_cds$cds_name

# define ranges to plot
csd_chr <- as.character(seqnames(csd_transcripts))[[1]]
csd_start <- min(start(csd_cds))
csd_end <- max(end(csd_cds))

# set up track schemes
set1 <- viridis::viridis_pal()(5)

my_scheme <- list(
    shape = "smallArrow",
    background.title = "transparent",
    fontface.title = 1,
    col.title = set1[1],
    col.frame = set1[1],
    thinBoxFeature = c("lincRNA"),
    fontcolor.group = set1[1],
    cex.title = 1,
    cex.group = 0.5,
    showTitle=FALSE
)

grt_scheme <- as.list(c(
    my_scheme,
    col.frame = set1[2],
    fontcolor.group = set1[2],
    fontcolor.title = set1[2],
    size = 1,
    col = set1[2],
    col.line = set1[2],
    fill = set1[2],
    fontface.group = 3,
    transcriptAnnotation = "symbol",
    cex.title = 0
))

aln_scheme <- as.list(c(
    my_scheme,
    fontcolor.title = set1[3],
    type = c("coverage", "pileup")))

# the genome
gat <- GenomeAxisTrack(showTitle = FALSE,
                       col = set1[4],
                       fontcolor = set1[4],
                       cex.title = 1,
                       col.title = set1[4],
                       name = csd_chr,
                       cex = 0.5)

# sequence track for reference
fasta_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
fa <- Biostrings::readDNAStringSet(fasta_file)
names(fa) <- gsub("^([^[:space:]]+).*", "\\1", names(fa)) 
st <- SequenceTrack(fa, chromosome = csd_chr)

# make tracks
csd_grt <- GeneRegionTrack(txdb,
                           chromosome = csd_chr,
                           start = csd_start,
                           end = csd_end,
                           name = "NCBI transcripts")
displayPars(csd_grt) <- grt_scheme

# alignment track (too big, have to subset)
at1 <- AlignmentsTrack(range = "output/030_process-aln/merged.bam",
                       referenceSequence = st,
                       isPaired = TRUE,
                       name = "All individuals",
                       chromosome = csd_chr,
                       start = csd_start,
                       end = csd_end)
displayPars(at1) <- aln_scheme

# join tracks and highlight hypervariable region
ht1 <- HighlightTrack(trackList = list(at1, csd_grt, gat),
                      start = 11771976,
                      end = 11772216,
                      chromosome = csd_chr,
                      col = set1[5],
                      fill = set1[5])

# output
pdf("test/csd_illumina.pdf",
    width = 5.5,
    height = 2.95,
    pointsize = 8)

plotTracks(list(ht1),
           cex.main = 1,
           fontface.main = 4,
           add = TRUE,
           col.mates = set1[5])

dev.off()
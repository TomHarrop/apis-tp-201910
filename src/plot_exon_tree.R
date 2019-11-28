#!/usr/bin/env Rscript

library(data.table)
library(ggtree)
library(ggplot2)

dist_dt <- fread("output/060_derived-alleles/all-indivs_aa.dist", header = FALSE)
setnames(dist_dt, names(dist_dt), c("V1", dist_dt[, V1]))

# turn it into a tree
d <- as.dist(as.matrix(data.frame(dist_dt, row.names = "V1")))
nj <- ape::bionj(d)
tree_data <- treeio::as.treedata(nj)

# annotate
annot <- data.table(tip_label = nj$tip.label)
annot[, type_code := substr(tip_label, 1, 2)]
annot[, indiv := gsub("_.*", "", tip_label)]
dup_indivs <- annot[duplicated(annot, by = "indiv"), unique(indiv)]
annot[indiv %in% dup_indivs, tech_dup := TRUE]
annot[!indiv %in% dup_indivs, tech_dup := NA]

type_order <- c("BB" = "Betta bees",
                "DR" = "Single drone",
                "TY" = "Taylor's Pass",
                "WS" = "Worker pool")

annot[, type := plyr::revalue(type_code, type_order)]

gp <- ggtree(tree_data,
       layout = "circular") %<+% annot +
    theme(legend.position="right",
          text = element_text(size = 8)) +
    xlim(c(-0.1, 0.3)) +
    scale_colour_viridis_d(guide = guide_legend(title = NULL)) +
    geom_tiplab2(mapping = aes(colour = type),
                 offset = 0.01,
                 size = 3)

ggsave("test/tree.pdf", gp, device = cairo_pdf, width = 140, height = 75, units = "mm")

saveRDS(gp, "test/tree.Rds")



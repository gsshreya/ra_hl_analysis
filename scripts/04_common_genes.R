library(annotate)
library(hgu133a.db)

ra <- read.csv("data/processed/ra_deg_results.csv")

ra$gene <- getSYMBOL(ra[,1], "hgu133a")
ra <- ra[!is.na(ra$gene), ]

ra <- ra[order(-abs(ra$logFC)), ]
ra <- ra[!duplicated(ra$gene), ]

ra$significant <- ifelse(ra$adj.P.Val < 0.05 & abs(ra$logFC) > 1, "yes", "no")

genes_ra <- ra$gene
deg_ra <- ra$gene[ra$significant == "yes"]


hl <- read.csv("data/processed/hl_deg_results.csv")

hl$gene <- getSYMBOL(hl[,1], "hgu133a")
hl <- hl[!is.na(hl$gene), ]

hl <- hl[order(-abs(hl$logFC)), ]
hl <- hl[!duplicated(hl$gene), ]

hl$significant <- ifelse(hl$adj.P.Val < 0.05 & abs(hl$logFC) > 1, "yes", "no")

genes_hl <- hl$gene
deg_hl <- hl$gene[hl$significant == "yes"]


common_genes <- intersect(deg_ra, deg_hl)

write.csv(deg_ra, "results/ra_genes.csv", row.names = FALSE)
write.csv(deg_hl, "results/hl_genes.csv", row.names = FALSE)
write.csv(common_genes, "results/common_genes.csv", row.names = FALSE)

#summary
cat("RA genes:", length(genes_ra), "\n")
cat("RA significant genes:", length(deg_ra), "\n\n")

cat("HL genes:", length(genes_hl), "\n")
cat("HL significant genes:", length(deg_hl), "\n\n")

cat("Common genes:", length(common_genes), "\n")
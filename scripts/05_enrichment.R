library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

common <- read.csv("results/common_genes.csv")   # load shared genes
genes <- common[,1]                              # extract gene symbols

# convert gene symbols → Entrez IDs (required for enrichment)
gene_df <- bitr(genes,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

entrez <- gene_df$ENTREZID

# GO Biological Process enrichment
ego <- enrichGO(gene          = entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05)

# KEGG pathway enrichment
ekegg <- enrichKEGG(gene         = entrez,
                    organism     = "hsa",
                    pvalueCutoff = 0.05)


write.csv(as.data.frame(ego), "results/go_enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "results/kegg_enrichment.csv", row.names = FALSE)

dotplot(ego, showCategory = 15)
ggsave("plots/go_dotplot.png")

dotplot(ekegg, showCategory = 15)
ggsave("plots/kegg_dotplot.png")

head(as.data.frame(ego))
head(as.data.frame(ekegg))
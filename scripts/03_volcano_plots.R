library(ggplot2)
library(annotate)
library(hgu133a.db)

ra <- read.csv("data/processed/ra_deg_results.csv")

ra$gene <- getSYMBOL(ra[,1], "hgu133a")   # map probes → genes
ra <- ra[!is.na(ra$gene), ]

ra <- ra[order(-abs(ra$logFC)), ]         # strongest probe per gene
ra <- ra[!duplicated(ra$gene), ]

ra$significant <- ifelse(ra$adj.P.Val < 0.05 & abs(ra$logFC) > 1,
                         "yes", "no")

sig_count_ra <- sum(ra$significant == "yes")

p1 <- ggplot(ra, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal() +
  ggtitle("RA Volcano (Gene-level)") +
  annotate("text",
           x = min(ra$logFC) + 0.5,
           y = max(-log10(ra$P.Value)) - 0.5,
           label = paste("Significant genes:", sig_count_ra),
           hjust = 0,
           vjust = 1,
           size = 4)

ggsave("plots/ra_volcano_gene.png", plot = p1, width = 6, height = 5)

hl <- read.csv("data/processed/hl_deg_results.csv")

hl$gene <- getSYMBOL(hl[,1], "hgu133a")
hl <- hl[!is.na(hl$gene), ]

hl <- hl[order(-abs(hl$logFC)), ]
hl <- hl[!duplicated(hl$gene), ]

hl$significant <- ifelse(hl$adj.P.Val < 0.05 & abs(hl$logFC) > 1,
                         "yes", "no")

sig_count_hl <- sum(hl$significant == "yes")

p2 <- ggplot(hl, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("no" = "grey", "yes" = "blue")) +
  theme_minimal() +
  ggtitle("HL Volcano (Gene-level)") +
  annotate("text",
           x = min(hl$logFC) + 0.5,
           y = max(-log10(hl$P.Value)) - 0.5,
           label = paste("Significant genes:", sig_count_hl),
           hjust = 0,
           vjust = 1,
           size = 4)

ggsave("plots/hl_volcano_gene.png", plot = p2, width = 6, height = 5)
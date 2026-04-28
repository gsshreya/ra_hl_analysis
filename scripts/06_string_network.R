library(STRINGdb)   # STRING database access
library(igraph)     # network analysis

common <- read.csv("results/common_genes.csv")   # load shared genes
genes <- common[,1]

string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,        # human
  score_threshold = 400  # medium confidence
)

mapped <- string_db$map(
  data.frame(gene = genes),
  "gene",
  removeUnmappedRows = TRUE
)

interactions <- string_db$get_interactions(mapped$STRING_id)   # get network edges

g <- graph_from_data_frame(
  interactions[, c("from", "to")],
  directed = FALSE
)

# map STRING IDs → gene symbols
id_map <- mapped[, c("STRING_id", "gene")]
V(g)$name <- id_map$gene[match(V(g)$name, id_map$STRING_id)]

g <- delete_vertices(g, which(is.na(V(g)$name)))   # remove unmapped nodes

print(g)   # basic network summary

# centrality metrics (CytoHubba-style)
deg <- degree(g)
bet <- betweenness(g)
clo <- closeness(g)

hub_df <- data.frame(
  gene = names(deg),
  degree = deg,
  betweenness = bet,
  closeness = clo
)

hub_df <- hub_df[order(-hub_df$degree), ]   # rank by degree

top_hubs <- head(hub_df, 10)

print(top_hubs)

write.csv(interactions, "results/string_interactions.csv", row.names = FALSE)
write.csv(hub_df, "results/hub_metrics.csv", row.names = FALSE)
write.csv(top_hubs, "results/top_hub_genes.csv", row.names = FALSE)

# simple network plot
png("plots/string_network.png", width = 800, height = 800)
plot(g,
     vertex.size = 5,
     vertex.label = NA,
     main = "STRING Network (RA + HL shared genes)")
dev.off()
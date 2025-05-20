library("igraph", lib.loc = "/Storage/data1/jorge.munoz/NRGSC.new/libraries")
table <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_current.triples", header = F, sep = " ")
colnames(table) <- c("V1", "V2", "weight")

G <- simplify(graph_from_data_frame(table, directed = F), remove.multiple = T)
names <- V(G)[order(degree(G), decreasing = TRUE)][1:110]
hub <- degree(G)[names]
names(hub) <- names(degree(G)[names])
hub <- as.data.frame(hub)
write.table(rownames(hub), "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/new/newhub_genes.txt", col.names = F, quote = F, row.names = F)


#Nid <- V(G)$name
#Nid <- as.data.frame(Nid)
#write.table(Nid, "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/TF/network.id", col.names = T, quote = F)



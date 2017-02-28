library(data.table)
library(arules)
library(dplyr)
library(corrr)

#otu2 <- read.csv("C:/Users/borre_000/Desktop/otu_table_psn_v13.csv", row.names = 1)
otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
#metadat <- read.csv("C:/Users/borre_000/Desktop/v13_map_uniquebyPSN.csv")
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")
colnames(otu1)[1:4]
dim(otu1)

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
metadat[stoolsamp,]
spptab <- colnames(otu2) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu3 <- otu2[-which(rowSums(otu2[,spptab]) == 0),spptab]
#rules1 <- apriori(data = t(otu3), parameter = list(support = 0.1, confidence = 0.8, maxlen = 3))

#inspect(sort(rules1, by = "confidence"))
#plot(rules1, method = NULL, measure = "support", shading = "lift", interactive = FALSE, data = NULL, control = NULL)
#plot(rules1, method = "graph")
#subrules1 <- head(sort(rules1, by = "support"), 50)
#plot(subrules1, method = "graph")

otu3.1 <- ceiling(apply(otu3, 2, function(x) x/sum(x)))
otu3.2 <- apply(otu3, 2, function(x) x/sum(x))

dotu3 <- dist(t(otu3.2), upper = T)
d3 <- as.matrix(dotu3)
dim(d3)
d3[d3 >= quantile(d3)[4]] <- 1
d3[d3 < quantile(d3)[4]] <- 0
sum(d3)
dg1 <- graph.adjacency(d3)
plot(degree(dg1)[order(degree(dg1), decreasing = T)])

plot(rowMeans(otu3.2)~apply(otu3.2, 1, function(x) sum(x > 0)))
cor.test(rowMeans(otu3.2),apply(otu3.2, 1, function(x) sum(x > 0)))
cor.test(apply(otu3.2, 1, median),apply(otu3.2, 1, function(x) sum(x > 0)))
otu4 <- apply(otu3[which(apply(otu3, 1, function(x) sum(x != 0))>125), 1:100], 2, function(x) x/sum(x))
library(bipartite)


otu2.1 <- apply(apply(as.matrix(otu2)[,-ncol(otu2)], 2, as.numeric), 2, function(x) x/sum(x))
plot(rowMeans(otu2.1)~apply(otu2.1, 1, function(x) sum(x > 0)))



###################

library(tidyverse)
library(igraph)
library(ggraph)
library(corrr)
otu2 <- otu2[,-ncol(otu2)]
tidy_cors <- otu3.2 %>% correlate() %>% stretch()
head(tidy_cors)
hist(tidy_cors$r)
sum(tidy_cors$r > .5, na.rm = T)
nrow(tidy_cors)

graph_cors <- tidy_cors %>% filter(abs(r) > 0.5) %>% graph_from_data_frame(directed = F)
plot(graph_cors, node_labels = NA)

sampls <- which(paste("X",as.character(metadat$SampleID), sep = "")%in%names(V(graph_cors)))
ggraph(graph_cors) + geom_edge_link() + geom_node_point(col = as.numeric(metadat$sex[sampls])) + theme_graph()

library(data.table)
library(arules)
otu1 <- fread(input = "~/Desktop/otu_table_psn_v13.csv", header = TRUE)
otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv")
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")
colnames(otu1)[1:4]
dim(otu1)
rules <- apriori(t(otu2[,-2911]))

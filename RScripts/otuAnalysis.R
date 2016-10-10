library(data.table)
library(arules)
#otu1 <- fread(input = "~/Desktop/otu_table_psn_v13.csv", header = TRUE)
otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv")
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")
colnames(otu1)[1:4]
dim(otu1)

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
metadat[stoolsamp,]
spptab <- colnames(otu2) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu3 <- otu2[-which(rowSums(otu2[,spptab]) == 0),spptab]
rules1 <- apriori(data = t(otu3), parameter = list(support = 0.1, confidence = 0.8, maxlen = 3))

inspect(sort(rules1, by = "confidence"))
plot(rules1, method = NULL, measure = "support", shading = "lift", interactive = FALSE, data = NULL, control = NULL)
plot(rules1, method = "graph")
subrules1 <- head(sort(rules1, by = "support"), 50)
plot(subrules1, method = "graph")

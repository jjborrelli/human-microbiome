devtools::install_github("joey711/biomformat")
library(biomformat)
library(data.table)
library(gambin)

system.time(test <- read_hdf5_biom("~/Downloads/ag_fecal.biom"))
names(test)
length(test$data[[2]])

sapply(test$data, sum)
test$data[[1]]

test <- do.call(rbind, test$data)
dim(test)

nsp <- (apply(test, 2, function(x) sum(x != 0)))
nrds <- apply(test, 2, sum)

gav <- apply(test[,nrds >= 2000], 2, function(x){rev(sort(get_abundvec(x[x!=0], N = 2000)))})
fzag <- lapply(gav, fzmod)
fzag <- do.call(rbind, fzag)
plot(fzag$s~fzag$N, ylim = c(0,3.5), xlim = c(0, 600))

fzag1 <- lapply(gav[sapply(gav, length) >= 200], function(x) fzmod(x[1:200]))
fzag1 <- do.call(rbind, fzag1)


############################################################################################################
############################################################################################################
############################################################################################################
system.time(
ag <- read_hdf5_biom("~/Documents/AmericanGut/ag_10k_fecal.biom")
)

ag.tax <- t(sapply(ag$rows, function(x) x$metadata$taxonomy))
ag <- as.data.table(do.call(rbind, ag$data))

ag.meta <- fread("~/Documents/AmericanGut/ag_10k_fecal.txt", header = T)
colnames(ag.meta)[1] <- "SampleID"
which(ag.meta$SampleID %in% colnames(ag)[1])

ag.meta <- (ag.meta[match(colnames(ag), ag.meta$SampleID),])

fzag <- apply(ag, 2, function(x) fzmod(rev(sort(x[x>0]))))
fzag <- rbindlist(fzag)

gbag <- apply(ag, 2, function(x) fitGambin(rev(sort(x[x>0]))))

#


table(ag.meta$PD_whole_tree_10k)
boxplot(fzag$s~toupper(ag.meta$SUBSET_BMI))
# quick comparisons
t.test(fzag$s~toupper(ag.meta$SUBSET_HEALTHY))
t.test(fzag$N~toupper(ag.meta$SUBSET_HEALTHY))
cvec <- vector(length = length(toupper(ag.meta$SUBSET_HEALTHY)))
cvec <- ifelse(toupper(ag.meta$SUBSET_HEALTHY), "blue", "green4")
plot(fzag$s,as.logical(toupper(ag.meta$SUBSET_HEALTHY)), pch = 20)


colnames(ag.meta)[(grep("VIOSCREEN", colnames(ag.meta)))]
df1 <- data.frame(fiber = ag.meta$VIOSCREEN_FIBER, mv = ag.meta$VIOSCREEN_MULTIVITAMIN, yog = ag.meta$VIOSCREEN_D_YOGURT, fat = ag.meta$VIOSCREEN_FAT, prot = ag.meta$VIOSCREEN_PROTEIN, glu = ag.meta$VIOSCREEN_GLUCOSE, fruc = ag.meta$VIOSCREEN_FRUCTOSE)
df2 <- apply(df1, 2, function(x){x[x == "Unknown" | x == "Unspecified" | x == "no_data"] <- NA;return(as.numeric(x))})

vio <- ag.meta[,(grep("VIOSCREEN", colnames(ag.meta))), with = F]
vio <- apply(vio, 2, function(x){x[x == "Unknown" | x == "Unspecified"| x == "no_data"] <- NA;return((x))})
vio <- as.data.frame(vio[!apply(vio, 1, function(x) all(is.na(x))),])

summary(lm(s~fiber+yog+fat+prot+glu+fruc, data = data.frame(s = fzag$N[!apply(df2, 1, function(x) all(is.na(x)))], df2[!apply(df2, 1, function(x) all(is.na(x))),])))
ag.meta$VIOSCREEN_







############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################

otu2 <- read.csv("~/Desktop/otu_table_psn_v13.csv", row.names = 1)
metadat <- read.csv("~/Desktop/v13_map_uniquebyPSN.csv")

stoolsamp <- which(metadat$HMPbodysubsite == "Stool")
spptab <- colnames(otu2) %in% paste0("X",metadat[stoolsamp,]$SampleID)
otu3 <- otu2[-which(rowSums(otu2[,spptab]) == 0),spptab]
rm(otu2)
rm(metadat)

gav1 <- apply(otu3[apply(otu3, 2, sum) >= 2000], 2, function(x) rev(sort(get_abundvec(x, N = 2000))))
gavfz1 <- lapply(gav1, fzmod)
gavfz1 <- do.call(rbind, gavfz1)

plot(fzag$s~fzag$N, ylim = c(0,3.5), xlim = c(0, 1000), pch = 20)
points(gavfz1$s~gavfz1$N, col = "blue", pch = 20)
points(simfz6$s~simfz6$N, col = "green4", pch = 20)
points(simfz1$s~simfz1$N, col = "green3", pch  = 20)

ir50.6 <- inrange(trunc = 50, gav = gavsim6, realdat = gav)
sum(ir50.6, na.rm = T)/length(ir50.6[!is.na(ir50.6)])

ir100.6 <- inrange(trunc = 100, gav = gavsim6, realdat = gav)
sum(ir100.6, na.rm = T)/length(ir100.6[!is.na(ir100.6)])

ir200.6 <- inrange(trunc = 200, gav = gavsim6, realdat = gav)
sum(ir200.6, na.rm = T)/length(ir200.6[!is.na(ir200.6)])

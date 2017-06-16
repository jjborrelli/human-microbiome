devtools::install_github("joey711/biomformat")
library(biomformat)
library(data.table)
library(gambin)

system.time(ag <- read_hdf5_biom("~/Downloads/ag_fecal.biom"))
names(test)
length(test$data[[2]])

sapply(test$data, sum)
test$data[[1]]

test <- do.call(rbind, test$data)
dim(test)

nsp <- (apply(test, 2, function(x) sum(x != 0)))
nrds <- apply(ag, 2, sum) >= 2000

gav <- apply(ag[,nrds, with = FALSE], 2, function(x){rev(sort(get_abundvec(x[x!=0], N = 2000)))})
fzag <- lapply(gav, fzmod)
fzag <- do.call(rbind, fzag)
plot(log10(fzag$s)~log10(fzag$N))#, ylim = c(0,3.5), xlim = c(0, 600))

fzag1 <- lapply(gav[sapply(gav, length) >= 200], function(x) fzmod(x[1:200]))
fzag1 <- do.call(rbind, fzag1)

fs <- lapply(gav, function(x) fitsad(x, sad = "poilog")@fullcoef)
N <- sapply(gav, length)

############################################################################################################
############################################################################################################
############################################################################################################
system.time(
ag <- read_hdf5_biom("~/Documents/AmericanGut/ag_10k_fecal.biom")
)

ag.tax <- t(sapply(ag$rows, function(x) x$metadata$taxonomy))
ag <- as.data.table(do.call(rbind, ag$data))
ma.ag2 <- apply(ag, 2, function(x) max(x))
hist(ma.ag)
ag.meta <- fread("~/Documents/AmericanGut/ag_10k_fecal.txt", header = T)
colnames(ag.meta)[1] <- "SampleID"
which(ag.meta$SampleID %in% colnames(ag)[1])

ag.meta <- (ag.meta[match(colnames(ag), ag.meta$SampleID),])

fzag <- apply(ag, 2, function(x){
  ga <- get_abundvec(rev(sort(x[x>0])), 2000)
  if(length(ga > 25)){smod <-fzmod(ga)}else{smod <- data.frame(N = NA, s = NA, nll = NA, r2 = NA)}
  return(smod)
}) 

fzag <- rbindlist(fzag)
plot(fzag[,1:2])


gbag <- apply(ag, 2, function(x) fitGambin(rev(sort(x[x>0]))))

#


table(ag.meta$PD_whole_tree_10k)
boxplot(fzag$s~toupper(agm$COUNTRY), las = 2)
# quick comparisons
t.test(fzag$s~toupper(ag.meta$SUBSET_HEALTHY))
t.test(fzag$N~toupper(ag.meta$SUBSET_HEALTHY))
cvec <- vector(length = length(toupper(ag.meta$SUBSET_HEALTHY)))
cvec <- ifelse(toupper(ag.meta$SUBSET_HEALTHY), "blue", "green4")
plot(fzag$s,as.logical(toupper(ag.meta$SUBSET_HEALTHY)), pch = 20)


colnames(ag.meta)[(grep("VIOSCREEN", colnames(ag.meta)))]
df1 <- data.frame(fiber = ag.meta$VIOSCREEN_FIBER, mv = ag.meta$VIOSCREEN_MULTIVITAMIN, yog = ag.meta$VIOSCREEN_D_YOGURT, fat = ag.meta$VIOSCREEN_FAT, prot = ag.meta$VIOSCREEN_PROTEIN, glu = ag.meta$VIOSCREEN_GLUCOSE, fruc = ag.meta$VIOSCREEN_FRUCTOSE)
df2 <- apply(df1, 2, function(x){x[x == "Unknown" | x == "Unspecified" | x == "no_data"] <- NA;return(as.numeric(x))})

vio <- ag.meta[,(grep("VIOSCREEN", colnames(agm)))]
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

gav <- apply(otu3[apply(otu3, 2, sum) >= 2000], 2, function(x) rev(sort(get_abundvec(x, N = 2000))))
gavfz1 <- lapply(gav, fzmod)
gavfz1 <- do.call(rbind, gavfz1)

gavfz2 <- do.call(rbind, lapply(1:ncol(otu3), function(x) fzmod(rev(sort(otu3[,x][otu3[,x]>0])))))
N1 <- apply(otu3, 2, function(x) sum(x > 0))

plot(fzag$s~fzag$N, ylim = c(.5,4), xlim = c(0, 600), pch = 20)
points(gavfz1$s~gavfz1$N, col = "blue", pch = 20)
points(simfz6$s~simfz6$N, col = "green4", pch = 20)
points(simfz1$s~simfz1$N, col = "green3", pch  = 20)

ir50.6 <- inrange(trunc = 50, gav = gavsim6, realdat = gav)
sum(ir50.6, na.rm = T)/length(ir50.6[!is.na(ir50.6)])

ir100.6 <- inrange(trunc = 100, gav = gavsim6, realdat = gav)
sum(ir100.6, na.rm = T)/length(ir100.6[!is.na(ir100.6)])

ir200.6 <- inrange(trunc = 200, gav = gavsim6, realdat = gav)
sum(ir200.6, na.rm = T)/length(ir200.6[!is.na(ir200.6)])

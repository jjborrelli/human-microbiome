devtools::install_github("joey711/biomformat")
library(biomformat)

system.time(test <- read_hdf5_biom("~/Downloads/ag_fecal.biom"))

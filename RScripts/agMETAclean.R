# Script to clean up the American Gut metadata
# metadata file: ag_10k_fecal.txt
library(dplyr)

dd <- read.csv("~/Documents/AmericanGut/AG data dictionary - Sheet1.csv")
ag.meta <- fread("~/Documents/AmericanGut/ag_10k_fecal.txt", header = T)
agm <- ag.meta

agm[agm == "Unknown" | agm == "Unspecified" | agm == "no_data"] <- NA

perNA <- sapply(1:ncol(agm), function(x) sum(is.na(agm[[colnames(agm)[x]]]))/(agm[,.N])) != 1

agm <- agm[,colnames(agm)[perNA], with = FALSE]
#agm <- agm2

agm <- as.data.frame(agm)
for(i in 1:ncol(agm)){
  if(sum(colnames(agm) %in% dd$Column.name[i]) == 0){next}
  if(dd$Data.type[i] == "bool"){
    agm[,colnames(agm) %in% dd$Column.name[i]] <- toupper(agm[,colnames(agm) %in% dd$Column.name[i]])
    agm[,colnames(agm) %in% dd$Column.name[i]][agm[,colnames(agm) %in% dd$Column.name[i]] == "NO"] <- FALSE
    agm[,colnames(agm) %in% dd$Column.name[i]][agm[,colnames(agm) %in% dd$Column.name[i]] == "YES"] <- TRUE
    agm[,colnames(agm) %in% dd$Column.name[i]] <- as.logical(agm[,colnames(agm) %in% dd$Column.name[i]])
  }
  

  if(dd$Data.type[i] == "float"){
    agm[,colnames(agm) %in% dd$Column.name[i]] <- as.numeric(agm[,colnames(agm) %in% dd$Column.name[i]])
  }
  
  if(dd$Data.type[i] == "int"){
    agm[,colnames(agm) %in% dd$Column.name[i]] <- as.numeric(agm[,colnames(agm) %in% dd$Column.name[i]])
  }
  print(i)
}

agm$ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL <- toupper(agm$ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL)
agm$ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL[agm$ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL == "NO"] <- FALSE
agm$ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL[agm$ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL == "YES"] <- TRUE
agm$ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL <- as.logical(agm$ALCOHOL_TYPES_SPIRITSHARD_ALCOHOL)

agm$NON_FOOD_ALLERGIES_BEESTINGS <- toupper(agm$NON_FOOD_ALLERGIES_BEESTINGS)
agm$NON_FOOD_ALLERGIES_BEESTINGS[agm$NON_FOOD_ALLERGIES_BEESTINGS == "NO"] <- FALSE
agm$NON_FOOD_ALLERGIES_BEESTINGS[agm$NON_FOOD_ALLERGIES_BEESTINGS == "YES"] <- TRUE
agm$NON_FOOD_ALLERGIES_BEESTINGS <- as.logical(agm$NON_FOOD_ALLERGIES_BEESTINGS)

agm$SUBSET_ANTIBIOTIC_HISTORY <- toupper(agm$SUBSET_ANTIBIOTIC_HISTORY)
agm$SUBSET_ANTIBIOTIC_HISTORY[agm$SUBSET_ANTIBIOTIC_HISTORY == "NO"] <- FALSE
agm$SUBSET_ANTIBIOTIC_HISTORY[agm$SUBSET_ANTIBIOTIC_HISTORY == "YES"] <- TRUE
agm$SUBSET_ANTIBIOTIC_HISTORY <- as.logical(agm$SUBSET_ANTIBIOTIC_HISTORY)



for(i in grep("VIOSCREEN_HEI_", colnames(agm))){agm[,i] <- as.numeric(agm[,i])}

str(agm[,100:200])
max(agm$HEIGHT_CM, na.rm  = T)

agm.s <- data.frame(agm, s = fzag$s)
fit <- lm(s~SALTED_SNACKS_FREQUENCY+EXERCISE_FREQUENCY+BOWEL_MOVEMENT_FREQUENCY+FERMENTED_PLANT_FREQUENCY+SUGARY_SWEETS_FREQUENCY+PREPARED_MEALS_FREQUENCY+HOMECOOKED_MEALS_FREQUENCY+VITAMIN_D_SUPPLEMENT_FREQUENCY+READY_TO_EAT_MEALS_FREQUENCY+FLOSSING_FREQUENCY+ONE_LITER_OF_WATER_A_DAY_FREQUENCY+COSMETICS_FREQUENCY+TEETHBRUSHING_FREQUENCY+FROZEN_DESSERT_FREQUENCY+ALCOHOL_FREQUENCY+POULTRY_FREQUENCY+SMOKING_FREQUENCY+VITAMIN_B_SUPPLEMENT_FREQUENCY+MILK_SUBSTITUTE_FREQUENCY+SEAFOOD_FREQUENCY+PROBIOTIC_FREQUENCY+OTHER_SUPPLEMENT_FREQUENCY+FRUIT_FREQUENCY+WHOLE_GRAIN_FREQUENCY+SUGAR_SWEETENED_DRINK_FREQUENCY+MILK_CHEESE_FREQUENCY+HIGH_FAT_RED_MEAT_FREQUENCY+RED_MEAT_FREQUENCY+MEAT_EGGS_FREQUENCY+VEGETABLE_FREQUENCY+POOL_FREQUENCY, data = agm.s)
summary(fit)

SUBSET_HEALTHY + AGE_CAT + 

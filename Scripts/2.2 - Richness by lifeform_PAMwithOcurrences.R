####Get richness by lifeform ####
library(dplyr)
library(data.table)
library(terra)
library(pbapply)

#Create directory to save results
dir.create("Data/PAM_indices_withOccurrences/")

#Load species information
spinfo <- fread("Data/SpeciesData_v2.gz")
spinfo$species <- gsub(" ", "_", spinfo$species)

#Import PAM
PAM <- fread("Data/PAM_with_occurrences.gz", data.table = FALSE)

#Check species
spp_pam <- PAM %>% dplyr::select(-x, -y) %>% colnames()
spp_pam_in <- intersect(spp_pam, spinfo$species)
setdiff(spp_pam, spinfo$species) #Should be 0

# Get richness for All species
richness_all <- PAM %>% dplyr::select(-x, -y) %>% rowSums()

#Calculate PAM indices by lifeform
table(spinfo$lifeForm)
#Get lifeforms
lf <- unique(spinfo$lifeForm)
lf

richness_lifeforms <- pblapply(lf, function(i){
  # Get species
  sp_i <- spinfo %>% filter(lifeForm == i) %>% pull(species)
  # Get species in PAM
  sp_i <- intersect(sp_i, colnames(PAM))
  #Get PAM of lifeforms i
  pam_i <- PAM %>% dplyr::select(sp_i)
  # Get richness
  rowSums(pam_i)
})
names(richness_lifeforms) <- lf

#Append results for all
richness_lifeforms[["All"]] <- richness_all

# Arrange list
names(richness_lifeforms) %>% dput()
lifeform_order <- c("All", "Tree", "Liana", "Shrub", "Subshrub",
                    "Terrestrial_herb", "Epiphytic_herb", 
                    "Bamboo", "Palm_tree")
richness_lifeforms <- richness_lifeforms[lifeform_order]

# Append xy
richness_lifeforms[["xy"]] <- PAM %>% select(x, y)

# Save
saveRDS(richness_lifeforms, "Data/Richness_by_lifeform_withOccurrences.rds")

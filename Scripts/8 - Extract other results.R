#Get informations
pam <- readRDS("Data/PAM.RDS")
spp <- pam %>% dplyr::select(-x, -y) %>% colnames()
head(spp)
length(spp)
#Get information of species
sp.info <- read.csv("Data/SpeciesData.csv")
head(sp.info$species)
sp.info$species <- gsub(" ", "_", sp.info$species)
sp.info <- sp.info %>% filter(species %in% spp)
table(sp.info$lifeForm)

####Ecoregions####
#Get richness of lifeforms
lf_indices <- list.files("Data/PAM_indices/", full.names = TRUE)
#Remove others
lf_indices <- lf_indices[!grepl("Other", lf_indices)]

#Read data
lf_indices <- pblapply(lf_indices, readRDS)

#Get lifeforms
lf_names <- sapply(lf_indices, function(x) x$lifeform)
names(lf_indices) <- lf_names
#Reorder lifeforms
lf_names <- c("All", "Tree", "Liana", "Shrub", "Subshrub", "Herb")

#Reorder lf_indices
lf_indices <- lapply(lf_names, function(x) lf_indices[[x]])
names(lf_indices) <- lf_names

#Rasterize indices
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_dissolved.gpkg")
r_base <- rast(ext = ext(af), res = 0.08333333)
#Test
#x <- lf_indices[[1]]
r <- pblapply(lf_indices, function(x){
  rx <- rasterize(x$xy, r_base,
                  values = x$Richness)
  rx[rx == 0] <- NA
  names(rx) <- "Richness"
  return(rx)
}) %>% rast()
plot(r[[1]])


af_eco <- vect("../spatial_files/Data/Ecoregions_Atlantic_Forest_simplified.gpkg")
plot(af_eco)

#Extract mean richness by ecoregions
mean_list <- pblapply(r, function(i){
  mean_eco <- extract(i, af_eco, fun = "mean", na.rm = TRUE, ID = F)
            })
mean_eco <- do.call("cbind", mean_list) %>% 
  mutate(Ecoregion = af_eco$ECO_NAME, .before = 1)
mean_eco

#Extract top3 ecoregions for each lifeform
top3 <- pblapply(lf_names, function(i){
  top_n(mean_eco[, c("Ecoregion", i)], n = 3) %>% 
    pull(Ecoregion)
})
names(top3) <- lf_names

#### Build PAM in Atlantic Forest ####
library(dplyr)
library(terra)
library(biosurvey)
library(pbapply)
library(tidyverse)
library(data.table)

#Set folder
setwd("C:/Users/wever/Desktop/kumodels")

####Import AF vector ####
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_limite_integrador.gpkg")
plot(af)
#Buffer of 10km
af <- buffer(af, width = 10*1000)
plot(af)

#Import any raster base to rasterize AF
#Current neotropic variables
var <- list.files("Current_Neotropic/", pattern = "\\.tif", full.names = T) %>%
  rast()
names(var)

#Get r base (with background 0)
r_base <- rasterize(af, var[[1]])
r_base <- crop(r_base, af, mask = TRUE)
r_base <- r_base * 0
plot(r_base)
#Rasterize and aggregate afr
f <- 0.08333333/res(r_base)[1]
afr <- aggregate(r_base, factor = f, fun = "min")
res(afr)
#Create directory to save
dir.create("Vetores")
writeRaster(afr, "Vetores/AF_raster.tiff", overwrite = TRUE)

#Pack spatial objects
# pack_afr <- terra::wrap(afr)
# pack_af <- terra::wrap(af)

#Get species
spp <- list.dirs("kuenm_models/", full.names = F, recursive = F)

#Create dataframe with path to models
#Tiff
tiff <- list.files(path = "kuenm_models/", pattern = "tif", recursive = TRUE,
                   full.names = TRUE)
#Get species and save as dataframe
names(tiff) <- dirname(tiff) %>% basename()
dt <- data.frame(species = names(tiff),
                 pathtiff = tiff, row.names = NULL)
any(duplicated(dt$species))

#RDS
rd <- list.files(path = "kuenm_models/", pattern = "itt.*", recursive = TRUE,
                 full.names = TRUE)
names(rd) <- dirname(rd) %>% basename()
dr <- data.frame(species = names(rd),
                 pathrds = rd, row.names = NULL)
any(duplicated(dr$species))
#Merge
df <- left_join(dt, dr, by = "species")

#Get species
spp <- df$species

# #Make cluster?
# library(parallel)
# cl <- makeCluster(45)
# clusterExport(cl, varlist= c("spp", "pack_afr", "pack_af"), #Get all objects created untill now
#               envir=environment())
# clusterEvalQ(cl, {
#   library(dplyr)
#   library(pbapply)
#   library(terra)
#   library(data.table)
# })

sp_l <- pblapply(seq_along(spp), function(i){
  #Erro com Tillandsia_xiphioides
  try({
  #Get specie
  sp <- spp[i]
  #print(sp)
  sp_dir <- file.path("kuenm_models/", sp)
  
  #Get species paths
  d_i <- df %>% filter(species == sp)
  
  #Unpack spatial objects
  # af <- unwrap(pack_af)
  # afr <- unwrap(pack_afr)
  #Get raster (continuous)
  sp_r <- rast(d_i$pathtiff)
  #Get threshold (10%)
  fm <- readRDS(d_i$pathrds)
  thr10 <- fm$thresholds$consensus$median
  r10 <- terra::app(sp_r, function(x) ifelse(x < thr10, 0, 1))
  
  #Sum rasters
  r10_res <- resample(r10, afr, method = "near")
  r10_res[is.na(r10_res)] <- 0
  r_final <- terra::trim(crop(r10_res, af, mask = T))
  names(r_final) <- sp
  # #Convert to dataframe
  # r_df <- as.data.frame(r_final, xy = F, cells = T) %>% 
  #   as.data.table()
  gc()
  return(r_final) })
})

# Remove errors
r_all <- sp_l[sapply(sp_l, function(x) class(x) == "SpatRaster")]
# Rasterize
rast_all <- rast(r_all)
# Convert to PAM
d <- terra::as.data.frame(rast_all, xy = TRUE)


#Save pam
library(data.table)
fwrite(d, "C:/Users/wever/Desktop/GitHub/Richness_patterns/Data/PAM.gzip", 
       compress = "gzip", row.names = FALSE)

####Rasterize some specie just to check####
d <- fread("C:/Users/wever/Desktop/GitHub/Richness_patterns/Data/PAM.gzip")

sp <- "Araucaria_angustifolia"

spp <- d %>% dplyr::select(x, y, sp) %>% as.matrix()

spr <- rasterize(spp[,c("x","y")],
                 afr, values = spp[,3])
plot(spr)


#### Get PAMS of other species included in the second round of review ####
# Species from napibio project
# Import first PAM
d <- fread("C:/Users/wever/Desktop/GitHub/Richness_patterns/Data/PAM.gzip")

# Get species
# 175 species
spp2 <- list.dirs("models_from_napibio/", full.names = FALSE, recursive = FALSE)

# Looping in the species
sp_l2 <- pblapply(seq_along(spp2), function(i){
  #Erro com Tillandsia_xiphioides
  try({
    #Get specie
    sp <- spp2[i]
    #print(sp)
    sp_dir <- file.path("models_from_napibio/", sp)
    
    #Get raster (binarized)
    r10 <- rast(file.path(sp_dir, "Current_Median_bin.tiff"))
   
    #Sum rasters
    r10_res <- resample(r10, afr, method = "near")
    r10_res[is.na(r10_res)] <- 0
    r_final <- terra::trim(crop(r10_res, af, mask = T))
    names(r_final) <- sp
    # #Convert to dataframe
    # r_df <- as.data.frame(r_final, xy = F, cells = T) %>% 
    #   as.data.table()
    gc()
    return(r_final) })
})

# Remove errors
r_all2 <- sp_l2[sapply(sp_l2, function(x) class(x) == "SpatRaster")]
# Rasterize
rast_all2 <- rast(r_all2)
# Convert to PAM
d2 <- terra::as.data.frame(rast_all2, xy = TRUE)

#Save pam
library(data.table)
fwrite(d2, "C:/Users/wever/Desktop/GitHub/Richness_patterns/Data/PAM_napibio.gzip", 
       compress = "gzip", row.names = FALSE)

# Merge PAMS
d_final <- left_join(d, d2, by = c("x", "y"))

# Save
fwrite(d_final, "C:/Users/wever/Desktop/GitHub/Richness_patterns/Data/PAM_final.gzip", 
       compress = "gzip", row.names = FALSE)

#### Get PAMS of species with more than 10 records to second round of review ####

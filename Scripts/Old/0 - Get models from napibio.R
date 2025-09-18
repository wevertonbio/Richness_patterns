#### Add more species ####
library(dplyr)
library(kuenm2)
library(terra)
library(data.table)
library(pbapply)

# Create directory to save new models from napibio
dir.create("C:/Users/wever/Desktop/kumodels/models_from_napibio")

# See all species
d <- fread("Data/SpeciesData.gz")

# See species with models
spp <- list.files(path = "C:/Users/wever/Desktop/kumodels/kuenm_models/",
                  pattern = "Current_Median",
                  recursive = TRUE, full.names = FALSE)
spp <- dirname(spp)
head(spp)

# Remove underline
spp <- gsub("_", " ", spp)

# See species without models
spp_out <- setdiff(d$species, spp)

# See species with models in napibio
spp_napi <- list.files("../napibio/Models/Plants/",
                       pattern = "FinalModels",
                       recursive = TRUE, full.names = FALSE)
spp_napi <- dirname(spp_napi)

# See species with models, but with models now in napibio
spp_in <- intersect(spp_out, spp_napi)

# See species withou m
spp_with_m <- list.files("../napibio/Models/Plants/",
                         pattern = "m.gpkg",
                         recursive = TRUE, full.names = FALSE)
spp_with_m <- dirname(spp_with_m)
no_m <- setdiff(spp_in, spp_with_m)

# For each species, rasterize occurrence and save in a folder
# Just to check a random species with models
r <- rast("C:/Users/wever/Desktop/kumodels/kuenm_models/Araucaria_angustifolia/Current_Median.tiff")
plot(r)
res(r)

####Import AF vector ####
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_limite_integrador.gpkg")
plot(af)
#Buffer of 10km
af <- buffer(af, width = 10*1000)
plot(af)

# Import variables
v <- rast("../napibio/Variables/Present/Variables.tif")
# Mask to Atlantic Forest
v_af <- crop(v, af, mask = TRUE)
plot(v_af[[1]])

#Loop in the species and predict to Atlantic forest
pblapply(spp_in, function(sp){
  print(sp)
  # sp <- spp_in[1]
  
  # Get species dir
  sp_dir <- file.path("../napibio/Models/Plants/", sp)
  
  # Get species fitted models
  fm <- readRDS(file.path(sp_dir, "FinalModels.rds"))
  
  # Get 10%
  fm <- fm[["10"]]
  
  # Predict models
  p <- predict_selected(fm, raster_variables = v_af)
  
  # Get M
  m <- vect(file.path(sp_dir, "m.gpkg"))
  
  # Mask predictions
  p <- mask(p$General_consensus$median, m, updatevalue = 0)
  
  # Binarize
  p_bin <- (p > fm$thresholds$consensus$median) * 1
  
  # Get occurrences
  occ <- fread(file.path(sp_dir, "Occurrences.gz"))
  
  # Save data
  sp_save <- file.path("C:/Users/wever/Desktop/kumodels/models_from_napibio", sp)
  dir.create(sp_save)
  
  writeRaster(p_bin,
              file.path(sp_save, "Current_Median_bin.tiff"), 
              overwrite = TRUE)
  fwrite(occ,
         file.path(sp_save, "Occurrences.gz"))
})

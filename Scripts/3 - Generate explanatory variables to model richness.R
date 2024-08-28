####Generate explanatory variables to model richness####

#Load packages
library(terra)
library(dplyr)
library(mapview)

#Load atlantic forest limits
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_limite_integrador.gpkg")
#Buffer of 50 km
af <- buffer(af, width = 10*1000)

#Calculate topographic heterogeneity based on altitude
alt <- rast("Current_Neotropic/Elevation_1km/elevation_1KMmd_GMTEDmd.tif")
alt_af <- crop(alt, af, mask = TRUE, touch = TRUE)
plot(alt_af)
mapview(alt_af)
window_size <- 5
weights <- matrix(1, nrow = window_size, ncol = window_size)
topo_het <- terra::focal(x = alt_af, w = weights, fun = sd,
                         na.policy = "all", pad = TRUE)
plot(topo_het)
#Aggregate to 5 arc-min
f <- 0.08333333/res(topo_het)
agg_het <- terra::aggregate(topo_het, fact = f, fun = mean)
names(agg_het) <- "Topo_het"
plot(agg_het)
mapview::mapview(agg_het, col.regions = pals::brewer.rdylgn(10))

#Salvar
dir.create("Others_variables")
writeRaster(agg_het, "Others_variables/Topo_heterogeneity.tiff", overwrite = T)

#Agregate to 5 arc-min considering range
f <- 0.08333333/res(alt_af)
range_alt_max <- terra::aggregate(alt_af, fact = f, fun = "max", na.rm = TRUE)
plot(range_alt_max)
range_alt_min <- terra::aggregate(alt_af, fact = f, fun = "min", na.rm = TRUE)
plot(range_alt_min)
range_alt <- range_alt_max - range_alt_min
plot(range_alt)


#Calculate topographic heterogeneity based on standard deviation of topoindex
#Download from ENVIREM
topoi <- rast("https://data.earthenv.org/topography/tpi_1KMmd_GMTEDmd.tif")
topoi_af <- crop(topoi, af, mask = TRUE, touch = TRUE)
plot(topoi_af)
#Aggregate to 5 arc-min and sum sd
f <- 0.08333333/res(topoi_af)
agg_topoi <- aggregate(topoi_af, fact = f, fun = sum)
names(agg_topoi) <- "Topoi_het"
plot(agg_topoi)
mapview::mapview(agg_topoi, col.regions = pals::brewer.rdylgn(10),
                 at = c(1, 10, 20, 30, 40, 50))

#Salvar
dir.create("Others_variables")
writeRaster(agg_topoi, "Others_variables/Topoindex_heterogeneity.tiff", overwrite = T)

####Calculate climate Stability ####
library(climateStability)

###Temperature
#After downloading the files, cut to AF extent, aggregate to 5 arc-min and save file
temp_dir <- "C:/Users/wever/Documents/Chelsa_Trace21/Bio01/"
dir.create(file.path(temp_dir, "Bio01_5_arc_min"))
all_temp <- list.files(path = temp_dir, full.names = T, recursive = F,
                       include.dirs = F, pattern = ".tif")
#Reorder data
names(all_temp) <- gsub("C:/Users/wever/Documents/Chelsa_Trace21/Bio01/CHELSA_TraCE21k_bio01_|_V1.0.tif", "", all_temp)
#Create sequence of layers (20 to -200)
time_slice <- as.numeric(names(all_temp))
#Get from 500 years
select_slice <- seq(-200, 20, 5) %>% sort(., decreasing = TRUE)
#Find indices
indices <- match(select_slice, names(all_temp))
select_temp <- all_temp[indices]
f_chelsa <- 0.08333333/0.008333333
#Looping to cut

library(pbapply)
pblapply(seq_along(select_temp), function(i){
  r <- rast(select_temp[i])
  name_r <- paste0("Temp", sprintf("%02d", i), "_",
                   2 - (names(select_temp[i]) %>% as.numeric() / 10),
                   "kBP")
  r_af <- crop(r, af, mask = TRUE)
  r_5 <- aggregate(r_af, fact = 10, FUN = mean)
  writeRaster(r_5, filename = paste0(temp_dir, "/Bio01_5_arc_min/", name_r, ".tiff"),
              overwrite = TRUE)
})
#Calculate deviation of temperature
temp_dev <- deviationThroughTime(variableDirectory = file.path(temp_dir, "/Bio01_5_arc_min/"),
                                 timeSlicePeriod = 500, fileExtension = ".tiff")
plot(temp_dev)
#Calculate stability
temp_stab <- climateStability::stabilityCalc(temp_dev)
names(temp_stab) <- "Temp_stab"
plot(temp_stab)
#Save
writeRaster(temp_stab, "Others_variables/Temperature_Stability.tiff",
            overwrite = TRUE)

###Precipitation
#After downloading the files, cut to AF extent, aggregate to 5 arc-min and save file
Prec_dir <- "C:/Users/wever/Documents/Chelsa_Trace21/Bio12/"
dir.create(file.path(Prec_dir, "Bio12_5_arc_min"))
all_Prec <- list.files(path = Prec_dir, full.names = T, recursive = F,
                       include.dirs = F, pattern = ".tif")
#Reorder data
names(all_Prec) <- gsub("C:/Users/wever/Documents/Chelsa_Trace21/Bio12/CHELSA_TraCE21k_bio12_|_V1.0.tif", "", all_Prec)
#Create sequence of layers (20 to -200)
time_slice <- as.numeric(names(all_Prec))
#Get from 500 years
select_slice <- seq(-200, 20, 5) %>% sort(., decreasing = TRUE)
#Find indices
indices <- match(select_slice, names(all_Prec))
select_Prec <- all_Prec[indices]
f_chelsa <- 0.08333333/0.008333333
#Looping to cut

library(pbapply)
pblapply(seq_along(select_Prec), function(i){
  r <- rast(select_Prec[i])
  name_r <- paste0("Prec", sprintf("%02d", i), "_",
                   2 - (names(select_Prec[i]) %>% as.numeric() / 10),
                   "kBP")
  r_af <- crop(r, af, mask = TRUE)
  r_5 <- aggregate(r_af, fact = 10, FUN = mean)
  writeRaster(r_5, filename = paste0(Prec_dir, "/Bio12_5_arc_min/", name_r, ".tiff"),
              overwrite = TRUE)
})
#Calculate deviation of Precipitation
Prec_dev <- deviationThroughTime(variableDirectory = file.path(Prec_dir, "/Bio12_5_arc_min/"),
                                 timeSlicePeriod = 500, fileExtension = ".tiff")
plot(Prec_dev)
#Calculate stability
Prec_stab <- climateStability::stabilityCalc(Prec_dev)
names(Prec_stab) <- "Prec_stab"
plot(Prec_stab)
#Save
writeRaster(Prec_stab, "Others_variables/Precipitation_Stability.tiff", 
            overwrite = TRUE)

####Aridity and PET####
#Downloaded from ENVIREM
ari <- rast("C:/Users/wever/Downloads/SAmerica_current_5arcmin_geotiff/current_5arcmin_aridityIndexThornthwaite.tif")
pet <- rast("C:/Users/wever/Downloads/SAmerica_current_5arcmin_geotiff/current_5arcmin_annualPET.tif")
aripet <- c(ari, pet)
ap_af <- crop(aripet, af, mask = T)
plot(ap_af)
cor(na.omit(as.data.frame(ap_af)))
names(ap_af) <- c("PET", "Aridity")
#Resample
ap_af <- resample(ap_af, agg_het)
#Save
writeRaster(ap_af, "Others_variables/PET_Aridity.tiff")

####MID-DOMAIN####
any_raster <- rast("Others_variables/Precipitation_Stability.tiff")
df_r <- as.data.frame(any_raster, xy = TRUE)
centroid_y <- df_r$y %>% median()
df_r$Distance_centroide <- abs(df_r$y - centroid_y)
r_mid <- rasterize(df_r %>% dplyr::select(x, y) %>% as.matrix(),
                   any_raster,
                   values = df_r$Distance_centroide)
plot(r_mid)
names(r_mid) <- "Mid_domain"
writeRaster(r_mid, "Others_variables/Mid_domain.tiff", overwrite = TRUE)

#####Eingenvector-based spatial filtering####
library(terra)
library(pbapply)
library(dplyr)
#Install older version of bigmds package
#remotes::install_version("bigmds", version = "2.0.1")
library(bigmds)
#Get coordinates
coords <- as.data.frame(r_mid, xy = TRUE) %>% dplyr::select(x, y) %>% as.matrix()

## calculating eigenvectors via principal coordinate analysis (K = numero de puntos - 1)
pcoa2 <- fast_mds(coords, l = 200, s_points = 5 * 10, r = 10, n_cores = 1)

eigs2 <- pcoa2$eig # all eigenvalues
pos_eigs2 <- 1:length(eigs2[eigs2 > 0]) # position of positive eigenvalues

eigens2 <- pcoa2$points[, pos_eigs2[1:100]] # getting only eigenvectors with positive eigenvalues


## selecting eigenvectors that should be considered in glms ********TWO EIGENVECTORS WERE SELECTED**********
### plots to select eigenvectors
plot(1:10, eigs2[1:10] / max(eigs2), type = "l", col = "black", las = 1, xlab = "", ylab = "") # all eigenvalues
abline(h = 0, col = "grey75", lwd = 1, lty = 2)
title(xlab = "Rank of eigenvalues", ylab = "", cex.lab = 1.2, line = 3)
title(xlab = "", ylab = "Eigenvalues (normalized)", cex.lab = 1.2, line = 3.3)
abline(v = 2, col = "red", lwd = 1, lty = 2)

#Get eingevectors
E1_2 <- eigens2[, 1]
E2_2 <- eigens2[, 2]
E3_2 <- eigens2[, 3]
#Rasterize eingevectors
resolution <- res(r_mid)
extension <- terra::ext(r_mid)
r_base <- rast(resolution = resolution, extent = extension)
r_res_1 <- rasterize(coords,
                     r_base,
                     values = E1_2)
plot(r_res_1) #Latitudinal variation
r_res_2 <- rasterize(coords,
                     r_base,
                     values = E2_2)
plot(r_res_2) #Longitudinal variation
r_res_3 <- rasterize(coords,
                     r_base,
                     values = E3_2)
plot(r_res_3) #??????variation
#Save variable sin object
spatial_variation <- c(r_res_1, r_res_2)
names(spatial_variation) <- c("latitudinal_eingevector", 
                              "longitudinal_eingevector")
plot(spatial_variation)
#save as environmental variables
writeRaster(spatial_variation, "Others_variables/Spatial_eingenvectors.tiff",
            overwrite = TRUE)



#### Build PAM in Atlantic Forest ####
library(dplyr)
library(terra)

# Import models binarized
m <- rast("Data/Binarized_models_v2.tiff")

# Convert to dataframe
d <- as.data.frame(m, xy = TRUE)

# Check nummer of species
d %>% select(-x, -y) %>% ncol()
# 10341 species

#Save
fwrite(d, "Data/PAM.gz")


#### Add occurrences of species with less than 10 records to PAM ####

# Import occurrences
occ <- fread("Data/Occurrences.gz")

# Get species info
sp <- fread("Data/SpeciesData_v2.gz")

# Get species with less than 10 records (not modeled)
sp10 <- sp %>% filter(n > 0 & n < 10)

# Subset occurrences
occ10 <- occ %>% filter(species %in% sp10$species)
occ10$species %>% unique() %>% length()

# Get raster base
b <- m[[1]] * 0
plot(b)

# Rasterize points
spp <- sp10$species
r_occ <- pblapply(spp, function(i){
  # i = spp[1]
  # Extract coordinates
  occ_i <- occ10 %>% filter(species == i) %>% select(x, y) %>% as.matrix()
  
  # Extract ID with records
  ids <- extract(b, occ_i, cells = TRUE)[["cell"]] %>% unique()
  b_sp <- b
  b_sp[ids] <- 1
  return(b_sp)
})
# Rasterize
r_occ <- rast(r_occ)
names(r_occ) <- spp
plot(r_occ[[1:5]])

# Mask to Atlantic Forest
r_occ_af <- mask(r_occ, b)

# Save
writeRaster(r_occ_af, "Data/Occurrences_rasterized.tiff", overwrite = TRUE)

# Get richness to check
richness <- app(r_occ_af, "sum", na.rm = TRUE)
plot(richness)

# Merge PAM
m2 <- c(m, r_occ_af)
# Standardize names (replace space with _)
names(m2) <- gsub(" ", "_", names(m2))

# Convert to PAM
d2 <- as.data.frame(m2, xy = TRUE)

# Check nummer of species
d2 %>% select(-x, -y) %>% ncol()
# 14238 species

#Save
fwrite(d2, "Data/PAM_with_occurrences.gz")

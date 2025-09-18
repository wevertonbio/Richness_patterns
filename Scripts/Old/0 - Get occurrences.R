
# Get number of records of species
library(dplyr)
library(data.table)
library(pbapply)
library(terra)
library(mapview)

# Importar dados de espécies
d <- fread("Data/SpeciesData_v2.gz")

# Standardize species
d$species <- gsub("_", " ", d$species)

# Get number of records for some modeled species
# First models
kuenm_models_occ1 <- list.files(path = "C:/Users/wever/Desktop/kumodels/kuenm_models/",
                               pattern = "Candidate_and_selected_models", 
                               full.names = TRUE, 
                               recursive = TRUE)

kuenm_models_occ1 <- pblapply(kuenm_models_occ1, function(i){
  # Get calibration data
  d_i <- readRDS(i)
  d_i <- d_i$data_xy[d_i$calibration_data$pr_bg == 1, ] %>% 
    mutate(species = d_i$species, .before = 1)
})
kuenm_models_occ1 <- rbindlist(kuenm_models_occ1)
# Update species
kuenm_models_occ1$species <- gsub("_", " ", kuenm_models_occ1$species)

# Second models
kuenm_models_occ2 <- list.files(path = "../Models_to_run/Models_to_run_done/",
                                pattern = "data.RDS", full.names = TRUE, 
                                recursive = TRUE)
kuenm_models_occ2 <- pblapply(kuenm_models_occ2, function(i){
  # Get calibration data
  d_i <- readRDS(i)
  d_i <- d_i$data_xy[d_i$calibration_data$pr_bg == 1, ] %>% 
    mutate(species = d_i$species, .before = 1)
})
kuenm_models_occ2 <- rbindlist(kuenm_models_occ2)
kuenm_models_occ2$species <- gsub("_", " ", kuenm_models_occ2$species)

# Merge
kuenm_models_occ <- bind_rows(kuenm_models_occ1, kuenm_models_occ2) %>% 
  distinct()
kuenm_models_occ$species %>% unique() %>% length()

# Get data from Missiones
missiones_occ <- fread("../Models_ARG_PAR/Occurrences_final.csv")
# Merge data
kuenm_missiones_occ <- bind_rows(kuenm_models_occ, missiones_occ) %>% 
  distinct()
kuenm_missiones_occ$species %>% unique() %>% length()

# Get data from napibio
napibio_occ <- list.files(path = "C:/Users/wever/Desktop/kumodels/models_from_napibio/",
                          pattern = "Occurrences", full.names = TRUE, recursive = TRUE)
napibio_occ <- pblapply(napibio_occ, fread) %>% rbindlist()

# Merge data
kuenm_missiones_napibio_occ <- bind_rows(kuenm_missiones_occ, napibio_occ) %>% 
  distinct()
kuenm_missiones_napibio_occ$species %>% unique() %>% length()

# Get data from models fit to second round of review
r2_occ <- list.files(path = "C:/Users/wever/Desktop/kumodels/more_models_richness/Models/",
                     pattern = "Occurrences.gz", recursive = TRUE, full.names = TRUE)
r2_occ <- pblapply(r2_occ, fread)
r2_occ <- rbindlist(r2_occ)

# Merge all data
all_occ <- bind_rows(kuenm_missiones_napibio_occ, r2_occ) %>% 
  distinct()
all_occ$species %>% unique() %>% length()

# Get number of records from Get and filter records
gfp <- fread("../The_Invisible_Species/Check_Points/8-Rescued_from_Sea.gz", 
             select = c("species", "decimalLongitude.new1", 
                        "decimalLatitude.new1"))
# Get only species not in kuenm/missiones/napibio/more models
spp_ready <- unique(all_occ$species)
spp_ready <- gsub("_", " ", spp_ready)
spp_not_ready <- setdiff(gfp$species, spp_ready)
#Subset
gfp <- gfp %>% filter(species %in% spp_not_ready)
# Fix colnames
colnames(gfp) <- c("species", "x", "y")

# Merge data
occ <- bind_rows(gfp, all_occ) %>% distinct()

# Standardize species
occ$species <- gsub("_", " ", occ$species)

# Get number of species
length(unique(occ$species))

# Get number of occurrences per specie
n_occ <- occ %>% count(species)

# Get species
setdiff(n_occ$species, d$species)
setdiff(d$species, n_occ$species)

# Join data
d2 <- left_join(d, n_occ)
# Number of species with NA
d2 %>% filter(is.na(n)) %>% nrow()

# Fill with NA
d2$n[is.na(d2$n)] <- 0

# Check if species is in Atlantic Forest
pts <- vect(occ, geom = c(x = "x", y = "y"), crs = "epsg:4326")
#Import AF
af <- vect("https://github.com/wevertonbio/spatial_files/raw/refs/heads/main/Data/AF_limite_integrador.gpkg")
# Rasterize
v <- rast("Current_Neotropic/Bio01.tif")
af_r <- rasterize(af, v) %>% trim()
plot(af_r)
# Extract values
in_af <- extract(af_r, pts)
occ$in_af <- !is.na(in_af$layer)

# Count number of records in AF
n_in_af <- occ %>% count(species, in_af)
# Arrange columns
n_in_af <- n_in_af %>% filter(in_af) %>%  select(species, n_in_af = n)

# Merge data
d2 <- left_join(d2, n_in_af)

# Update NA
d2 %>% filter(is.na(n_in_af)) %>% nrow()
d2$n_in_af[is.na(d2$n_in_af)] <- 0

# Add information on modeled species
pam <- fread("Data/PAM_final_v2.gzip", data.table = FALSE)
spp_modeled <- colnames(pam[ ,-c(1:2)])
spp_modeled <- gsub("_", " ", spp_modeled)
d2$modeled <- NA
d2$modeled[d2$species %in% spp_modeled] <- TRUE
d2$modeled[!(d2$species %in% spp_modeled)] <- FALSE
table(d2$modeled)

# Check inconsistencies
d2 %>% filter(modeled & n < 10) %>% nrow()
d2 %>% filter(modeled & n < 10) %>% View()
# Change to not modeles
d2$modeled[d2$modeled & d2$n < 10] <- FALSE

# Check species with more than 10 records and not modeles
d2 %>% filter(!modeled & n >= 10 & n_in_af >= 1) %>% nrow()
d2 %>% filter(!modeled & n >= 10& n_in_af >= 1) %>% View()

# See if we have these species in napibio
spp_weird <- d2 %>% filter(!modeled & n >= 10) %>% pull(species) %>% unique()
spp_napibio <- list.dirs("../napibio/Models/Plants/", recursive = FALSE, 
                         full.names = FALSE)
spp_weird_in <- intersect(spp_napibio, spp_weird)
# We don't ????

# Update number of records based on unique records in one pixel
v <- rast("Current_Neotropic/Bio01.tif")
res(v) * 111 #Resolution in km

# Get data of species with >10 records and not modelled
occ_weird <- occ %>% filter(species %in% spp_weird) %>% select(-kingdom)
# Extract cell number
occ_cell <- extract(v, occ_weird[, 2:3], cell = TRUE)[["cell"]]
occ_weird <- cbind(occ_weird, "cell" = occ_cell)
# Get unique records
unique_weird <- occ_weird %>% distinct(species, cell, .keep_all = TRUE)
# Get number of records
n_occ_weird <- unique_weird %>% count(species)

# Update this number
d3 <- rows_update(x = d2, y = n_occ_weird, by = "species")
# Check inconsistencies
d3 %>% filter(modeled & n < 10) %>% nrow()
d3 %>% filter(!modeled & n >= 10 & n_in_af >= 1) %>% nrow()
d3 %>% filter(!modeled & n >= 10 & n_in_af >= 1) %>% View()


# Save species information
fwrite(d3, "Data/SpeciesData_v2.gz")

# Import again
d3 <- fread("Data/SpeciesData_v2.gz")

# Import PAM again
pam <- fread("Data/PAM_final_v2.gzip", data.table = FALSE)

# Update PAM with modeled species
spp_in_pam <- colnames(pam) %>% gsub(" ", "_", .)
spp_modeled <- d3 %>% filter(modeled) %>% pull(species) %>% gsub(" ", "_", .)
spp_to_remove <- setdiff(spp_in_pam, spp_modeled)
d3 %>% filter(species %in% 
                gsub("_", " ", spp_to_remove)) %>% View()
#Update PAM
pam_v3 <- pam
colnames(pam_v3) <- gsub(" ", "_", colnames(pam_v3))
pam_v3 <- pam_v3 %>% select(-spp_to_remove)
# Check
ncol(pam) - length(spp_to_remove)
ncol(pam_v3)
d3 %>% filter(modeled) %>% nrow()
#Save
fwrite(pam_v3, "Data/PAM_final_v3.gzip")

# Subset species modeled
r <- rast("Data/Binarized_models.tiff")
names(r)
# Standardize names
names(r) <- gsub(" ", "_", names(r))
# Get species to keep
to_keep <- d3 %>% filter(modeled) %>% pull(species) %>% 
  gsub(" ", "_", .)
to_keep <- intersect(names(r), to_keep)
r_updated <- r[[to_keep]]

# Save
writeRaster(r_updated, "Data/Binarized_models_v2.tiff", overwrite = TRUE)

# Save occurrences
fwrite(occ, "Data/Occurrences.gz")

# # Save species to model
# to_model <- d3 %>% filter(!modeled & n >= 10 & n_in_af >= 1)
# # Get occurrences
# occ_to_model <- occ %>% filter(species %in% to_model$species)
# unique(occ_to_model$species) %>% length()
# # Create directory to save
# dir.create("C:/Users/wever/Desktop/kumodels/more_models_richness")
# fwrite(occ_to_model, "C:/Users/wever/Desktop/kumodels/more_models_richness/Occurrences.gz")
# 
# # Plot some weird specie
# sp <- spp_weird[1]
# sp <- "Aphelandra nemoralis"
# occ_sp <- occ %>% filter(species == sp)
# pts <- vect(occ_sp, geom = c(x = "x", y = "y"), crs = "epsg:4326")
# mapview(pts)

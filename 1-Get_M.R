#### Get M of species ####
library(grinnell)
library(dplyr)
library(terra)
library(pbapply)
#library(mapview)
library(data.table)
library(parallel)
library(doSNOW)
library(foreach)

# #Join all records in a single data.frame (plants and animals)
# p <- fread("Data/Plants/Check_points/K-Final_Occurrences.gz")
# p$kingdom <- "Plants"
# a <- fread("Data/Animals/Check_points/J - Final_Occurrences.gz")
# a$kingdom <- "Animals"
# pa <- bind_rows(a, p)
# fwrite(pa, "Data/Occurrence_plants_animals.gzip", row.names = FALSE)

#Set folder to napibio
setwd("more_models_richness")

#Import folder with the environmental variables at 2.5, 5 and 10 arc-minutes, and wrap variables
res_var <- c("2.5", "5", "10")

current_variables <- pblapply(res_var, function(x){
  dir_x <- paste0("paleoclim_", x, "//cur.tif")
  #terra::wrap(rast(dir_x))
})
names(current_variables) <- res_var

variables <- c("bio_1", "bio_10", "bio_11", "bio_12", "bio_13", "bio_14",
               "bio_15", "bio_16", "bio_17","bio_2", "bio_3",
               "bio_4", "bio_5", "bio_6", "bio_7")

projection_dir <- paste0("paleoclim_", res_var)
names(projection_dir) <- res_var

# Capturar a especie da linha de comando
args <- commandArgs(trailingOnly = TRUE)
sp <- args[1]  # A espécie será o primeiro argumento passado ao script

#Import occurrences
occ <- fread("Thinned_Occurrences.gz")

#Subset species
occ <- occ %>% filter(species == sp)

#Get variables according to resolution
res_sp <- ifelse(nrow(occ) < 500, "5", ifelse(
    nrow(occ) > 500, "10", NA))
#Current variables
current_variables_sp <- terra::rast(current_variables[[res_sp]])
#Past variables
projection_dir_sp <- projection_dir[res_sp]

#Keep one record per grid
occ_id <- extract(current_variables_sp[[1]], occ[,c("x","y")], cell = T)$cell
occ <- occ[!duplicated(occ_id), ]

#Remove occurences manually
#occ$id <- 1:nrow(occ)
# to_remove <- c(47, 48, 49)
# occ <- occ[-to_remove, ]

#Get distances between records to set initial buffer
pts <- vect(occ, geom = c(x = "x", y = "y"),
            crs="+init=epsg:4326")

#mapview(pts)
all.dist <- pts %>% #Get distance matrix
  terra::distance(unit = "km") %>% as.numeric()
di <- quantile(all.dist, 0.75) %>% as.numeric() + 100
#Set minoimum and maximum of distance
if(di < 500)
  di <- 500
if(di > 2500)
  di <- 2500

# Draw mcp and add buffer
b <- hull(pts) %>% buffer(width = di * 1000)
current_variables_sp <- crop(current_variables_sp, b, mask = TRUE) %>%
  trim()

#Grid of combinations to test M
comb_grid <- data.frame(dispersal_events = c(5, 10, 15, 20, 25, 30, 35, 40),
                        kernel_spread = c(2, 2, 2, 2, 3, 4, 6, 8))

periods <- c("lh", "mh", "eh", "yds", "ba", "hs1", "lgm", "30kya", "40kya",
             "50kya", "60kya", "70kya", "80kya", "90kya",
             "100kya", "110kya", "120kya","lig")

#Run M simulation
m <- m_simulations(data = occ, long = "x", lat = "y",
                   current_variables = current_variables_sp,
                   variables = variables,
                   scale = TRUE,
                   center = TRUE,
                   projection_dir = projection_dir_sp,
                   periods = periods,
                   pattern = ".tif",
                   initial_m_buffer = NULL,
                   suitability_threshold = 5,
                   starting_proportion = 0.75,
                   proportion_to_disperse = 1,
                   sampling_rule = "random",
                   dispersal_kernel = "normal",
                   kernel_spread = 2,
                   max_dispersers = 4,
                   dispersal_events = 25,
                   comb_grid = comb_grid,
                   replicates = 3,
                   threshold = 5,
                   set_seed = 1,
                   skip_extinction = TRUE,
                   results_by_event = FALSE,
                   results_by_scenario = FALSE,
                   remove_m_without_records = TRUE,
                   extra_buffer = 50,
                   progress_bar = TRUE,
                   verbose = TRUE)
#Get best M
# m$summary %>% View()
# mapview(m$Combination_1$m_final) + mapview(pts) +
#         mapview(m$Combination_5$m_final) +
#         mapview(m$Combination_7$m_final)

#mapview(m$m[[12]]$m) + mapview(pts)
#Get minimum number of records in surplus m = 0
mn <- min(m$summary$occ_outside[m$summary$surplus_m == 0])

#Number of records allowed to be outside M: 1%
occ_outside_limit <- floor(nrow(occ) * 0.05)
# #Set limitof records allowed to be outside M: 1, 2 or 3
# if(occ_outside_limit > 3) {occ_outside_limit <- 3 }
# if(occ_outside_limit < 1) {occ_outside_limit <- 1 }
# }

#Check if minimum records is available
if(mn > occ_outside_limit){
  occ_outside_limit <- mn
}

#Get best combination
best_comb <- which(m$summary$surplus_m <= 0 &
                     m$summary$occ_outside <= occ_outside_limit)

if(length(best_comb) == 0){
  best_comb <- which(m$summary$surplus_m == min(m$summary$surplus_m) &
                       m$summary$occ_outside <= occ_outside_limit)
}

bc<- paste0("Combination_", min(best_comb))

#Get M
best_m <- m[[bc]]
#mapview(best_m$m_final) + mapview(pts)

#Get records to save (only inside M)
occ_m <- best_m$occ %>%
  filter(inside_m) %>% dplyr::select(-inside_m)

#Identify best m in m summary
m$summary$best_m <- FALSE
m$summary$best_m[min(best_comb)] <- TRUE

#Save results
sp_path <- file.path("Models/", sp)
dir.create(sp_path, recursive = TRUE)
writeVector(best_m$m_final, filename = paste0(sp_path, "/m.gpkg"), overwrite = TRUE)
fwrite(m$summary, paste0(sp_path, "/m_report.gz"), compress = "gzip",
       row.names = FALSE)
fwrite(occ_m, paste0(sp_path, "/Occurrences.gz"), compress = "gzip",
       row.names = FALSE)

#Exit R and run this
# cat more_models_richness/species_list.txt | parallel -j 15 "Rscript 1-Get_M.R {} > more_models_richness/log/resultado_{}.log 2> more_models_richness/log/erro_{}.log"

#Check species ready
library(dplyr)
library(data.table)

# #Set folder
# setwd("more_models_richness")
# 
# d <- fread("Thinned_Occurrences.gz")
# #Get species
# spp <- unique(d$species)
# 
# #Get species ready
# spp_ready <- list.dirs("Models/", recursive = FALSE, full.names = F)
# spp_not_ready <- setdiff(spp, spp_ready)
# #Save as species list
# writeLines(spp_not_ready, "Data/species_list_3.txt")
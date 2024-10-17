####PAM indices by lifeform ####
library(dplyr)
library(data.table)
library(terra)
library(GPUmatrix)

#Create directory to save results
dir.create("Data/PAM_indices/")

#Load species information
spinfo <- fread("Data/SpeciesData.gz")
spinfo$species <- gsub(" ", "_", spinfo$species)

#Import PAM
PAM <- fread("Data/PAM.gzip")

#Check species
spp_pam <- PAM %>% dplyr::select(-x, -y) %>% colnames()
spp_pam_in <- intersect(spp_pam, spinfo$species)
setdiff(spp_pam, spinfo$species) #Should be 0

# #Remove species
# PAM2 <- PAM %>% dplyr::select(x, y, spp_pam_in)
# fwrite(PAM2, "Data/PAM.gzip", row.names = FALSE, compress = "gzip")

#Remove coordinates
pam_herbs <- PAM %>% dplyr::select(-x, -y) %>% as.matrix()

####Import function####
#pam_indices_gpu modified from biosurvey to run in GPU (faster!)
source("Scripts/Functions/PAM_indices_gpu.R")

#Calculate PAM indice from all species together
ind_all <- PAM_indices_gpu(PAM = pam_all, indices = "basic")
#Include in the list x and y
ind_all$xy <- PAM %>% dplyr::select(x, y) %>% as.matrix()
#Include in the list the lifeform
ind_all$lifeform <- "All"
#Empty cuda cache
torch::cuda_empty_cache()
#Save
saveRDS(ind_all, "Data/PAM_indices/Indices_All.rds")
rm(ind_all)

#Calculate PAM indices by lifeform
table(spinfo$lifeForm)
#Get lifeforms
lf <- unique(spinfo$lifeForm)
lf

pblapply(seq_along(lf), function(i){
  #Get species in lifeform i
  lf_i <- lf[i] 
  sp_i <- spinfo %>% filter(lifeForm == lf_i)
  spp_i <- subset(colnames(pam_herbs), colnames(pam_herbs) %in% sp_i$species)
  #Get PAM of lifeforms i
  pam_i <- PAM %>% dplyr::select(spp_i) %>% as.matrix()
  #Calculate index
  ind_i <- PAM_indices_gpu(PAM = pam_i, indices = "basic")
  #Include in the list x and y
  ind_i$xy <- PAM %>% dplyr::select(x, y) %>% as.matrix()
  #Include in the list the lifeform
  ind_i$lifeform <- lf_i
  #Save
  saveRDS(ind_i, paste0("Data/PAM_indices/Indices_", lf_i, ".RDS"))
  #Empty cache
  rm(ind_i)
  gc()
  torch::cuda_empty_cache()
})

#Calculate PAM indice merging herbs
spp_herbs <- spinfo %>% filter(grepl("herb", lifeForm))
spp_herbs <- intersect(spp_herbs$species, colnames(PAM))
pam_herbs <- PAM %>% dplyr::select(spp_herbs) %>% as.matrix()
ind_herbs <- PAM_indices_gpu(PAM = pam_herbs, indices = "basic")
#Include in the list x and y
ind_herbs$xy <- PAM %>% dplyr::select(x, y) %>% as.matrix()
#Include in the list the lifeform
ind_herbs$lifeform <- "Herb"
#Empty cuda cache
torch::cuda_empty_cache()
#Save
saveRDS(ind_herbs, "Data/PAM_indices/Indices_Herb.rds")
rm(ind_herbs)


####Not used####
# ####Run Null_matrix_in_cluster and Null_matrix_in_cluster_by_lifeform to get the null Matrix####
# dir.create("PAM_indices/Dispersion_sign/")
# ####Calculate indices from null matrix####
# lf <- list.dirs("PAM_indices/Null_Matrix/", recursive = F, full.names = F)
# 
# lapply(seq_along(lf), function(x){
#   lf_i <- lf[x]
#   print(lf_i)
#   lf_dir <- file.path("PAM_indices/Null_Matrix", lf_i)
#   #Get random matrix
#   random_l <- list.files(lf_dir, pattern = "m", full.names = T)
#   pblapply(seq_along(random_l), function(i){
#     #Check if indice has already been calculated
#     has_file <- file.exists(paste0(lf_dir, "/r_ind", i, ".rds"))
#     if(!has_file){
#     random_m <- readRDS(paste0(lf_dir, "/m", i, ".rds"))
#     ind_i <- PAM_indices_gpu(PAM = random_m, indices = "DF")
#     saveRDS(ind_i, paste0("PAM_indices/Null_Matrix/", lf_i, "/r_ind", i, ".rds"))
#     rm(ind_i)
#     gc();torch::cuda_empty_cache();gc();torch::cuda_empty_cache()
#     }
#     })
#   })
# 
# ####Get percentis and calculate dispersion field sign####
# lf <- list.dirs("PAM_indices/Null_Matrix/", recursive = F, full.names = F)
# pblapply(seq_along(lf), function(i){
#   lf_i <- lf[i]
#   print(lf_i)
#   ind_i <- readRDS(paste0("PAM_indices/Indices_", lf_i, ".rds"))
#   ind_null_l <- list.files(path = paste0("PAM_indices/Null_Matrix/",
#                                          lf_i), 
#                            pattern = "r_ind", full.names = T)
#   ind_null_i <- lapply(seq_along(ind_null_l), function(i){
#     return(readRDS(ind_null_l[i]))
#   })
#   #Get dispersion field normalized by site and richness
#   ds <- (ind_i$Dispersion_field/ind_i$One_value_indices["Sites_Cells",])
#   ds <- ds/ind_i$One_value_indices["Species",]
#   ds[which(is.na(ds))] <- 0
#   
#   ds_null <- sapply(seq_along(ind_null_i), function(i){
#     ind_null_i_i <- ind_null_i[[i]]
#     n_sites <- ind_null_i_i$One_value_indices["Sites_Cells",]
#     S <- ind_null_i_i$One_value_indices["Species",]
#     ds_null_i <- ind_null_i_i$Dispersion_field/n_sites
#     ds_null_i <- ds_null_i/S
#     ds_null_i[which(is.na(ds_null_i))] <- 0
#     return(ds_null_i)
#   })
#   
#   ind_pos <- sapply(1:nrow(ds_null), function(i){
#     ds_i <- ds[i]
#     ds_null_i <- ds_null[i,]
#     pos5 <- quantile(ds_null_i, 0.05)
#     pos95 <- quantile(ds_null_i, 0.95)
#     pos_i <- ifelse(ds_i < pos5, -1,
#                     ifelse(ds_i > pos95, 1, 0))
#     return(pos_i)
#   })
#   
#   #Join with coordinates to plot
#   disp_sign <- PAM %>% as.data.frame() %>% dplyr::select(x, y) %>% 
#     mutate(DispersionSign = ind_pos)
#   #Join with Richness, dispersal field normalized and richness normalized
#   disp_sign <- disp_sign %>% 
#     mutate(Richness = ind_i$Richness,
#       NormalizedRichness = ind_i$Richness_normalized,
#            DispersedFieldNormalized = ds)
#   
#   # #Rasterize to see
#   afr <- rast("Vetores/AF_raster.tiff")
#   disp_sign_r <- rasterize(x = disp_sign %>% dplyr::select(x, y) %>% as.matrix(),
#                            y = afr,
#                            values = disp_sign$DispersionSign)
# 
#   plot(disp_sign_r, main = lf_i)
#   
#   #Salvar dataframe
#   write.csv(disp_sign,
#             paste0("PAM_indices/Dispersion_sign/", lf_i, ".csv"),
#             row.names = F)
#   
# })
# 
# 




# #Get NULL metrics of trees
# #Get random matrix in paralelle
# library(parallel)
# cl <- makeCluster(40)
# clusterExport(cl, varlist = "pam_herbs", #Get PAM2
#               envir=environment())
# clusterEvalQ(cl, {
#   library(picante)
#   library(pbapply)
# })
# 
# 
# dir.create("PAM_indices/Null_Matrix/All/")
# pblapply(1:999, function(i){
#   set.seed(i)
#   random_matrix_i <- picante::randomizeMatrix(samp = pam_herbs, 
#                                               null.model = "trialswap",
#                                               iterations = length(pam_herbs)*10)
#   saveRDS(random_matrix_i, paste0("PAM_indices/Null_Matrix/All/m", i, ".rds"))
# }, cl = cl)
# stopCluster(cl)
# 
# #get random matrix
# random_l <- list.files("PAM_indices/Null_Matrix/All/", full.names = T, 
#                        pattern = "m")
# 
# #Calculate PAM indices from NULL matrix
# pblapply(seq_along(random_l), function(i){
#   random_m <- readRDS(random_l[i])
#   ind_i <- PAM_indices_gpu(PAM = random_m, indices = "DF")
#   saveRDS(ind_i, paste0("PAM_indices/Null_Matrix/All/", "r_ind", i, ".rds"))
#   return(ind_i) })
# 
# 
# ####Get PAM of trees####
# sp_i <- spinfo %>% filter(lifeForm == "Tree")
# spp_i <- subset(colnames(PAM), colnames(PAM) %in% sp_i$species)
# pam_i <- PAM %>% dplyr::select(spp_i) %>% as.matrix()
# 
# #Get indices from trees
# ind_is <- PAM_indices_gpu(PAM = pam_i, indices = "all")
# #Save
# saveRDS(ind_is, "PAM_indices/Indices_is.rds")
# rm(ind_is)
# 
# #Get NULL metrics of trees
# #Get random matrix in paralelle
# library(parallel)
# cl <- makeCluster(40)
# clusterExport(cl, varlist = "pam_i", #Get PAM2
#               envir=environment())
# clusterEvalQ(cl, {
#   library(picante)
#   library(pbapply)
# })
# 
# #dir.create("PAM_indices/Null_Matrix")
# dir.create("PAM_indices/Null_Matrix/Trees")
# pblapply(1:999, function(i){
#   set.seed(i)
#   random_matrix_i <- picante::randomizeMatrix(samp = pam_i, 
#                                               null.model = "trialswap",
#                                               iterations = length(pam_i)*10)
#   saveRDS(random_matrix_i, paste0("PAM_indices/Null_Matrix/Trees/m", i, ".rds"))
# }, cl = cl)
# 
# #get random matrix
# random_l <- list.files("PAM_indices/Null_Matrix/Trees/", full.names = T)
# 
# #Calculate PAM indices from NULL matrix
# pblapply(seq_along(random_l), function(i){
#   random_m <- readRDS(random_l[i])
#   ind_i <- PAM_indices_gpu(PAM = random_m, indices = "DF")
#   saveRDS(ind_i, paste0("PAM_indices/Null_Matrix/Trees/", "r_ind", i, ".rds"))
#   return(ind_i) })
# stopCluster(cl)
# 
# 
# ####Get only herbs####
# sp_herb <- spinfo %>% filter(lifeForm == "Herb")
# #Get PAM of herbs
# #PAM <- fread("PAM.csv")
# spp_herb <- subset(colnames(PAM), colnames(PAM) %in% sp_herb$species)
# pam_herb <- PAM %>% dplyr::select(spp_herb) %>% as.matrix()
# 
# #Get indices from herbs
# ind_herbs <- PAM_indices_gpu(PAM = pam_herb, indices = "all")
# #Save
# saveRDS(ind_herbs, "PAM_indices/Indices_Herbs.rds")
# 
# #Get NULL metrics of herbs
# #Get random matrix in paralelle
# library(parallel)
# cl <- makeCluster(40)
# clusterExport(cl, varlist = "pam_herb", #Get PAM2
#               envir=environment())
# clusterEvalQ(cl, {
#   library(picante)
#   library(pbapply)
# })
# 
# #dir.create("PAM_indices/Null_Matrix")
# dir.create("PAM_indices/Null_Matrix/Herbs")
# pblapply(1:999, function(i){
#   set.seed(i)
#   random_matrix_i <- picante::randomizeMatrix(samp = pam_herb, 
#                                               null.model = "trialswap",
#                                               iterations = length(pam_herb)*10)
#   saveRDS(random_matrix_i, paste0("PAM_indices/Null_Matrix/Herbs/m", i, ".rds"))
# }, cl = cl)
# stopCluster(cl)
# #get random matrix
# random_l <- list.files("PAM_indices/Null_Matrix/Herbs/", full.names = T,
#                        pattern = "m")
# 
# #Calculate PAM indices from NULL matrix
# pblapply(seq_along(random_l), function(i){
#   random_m <- readRDS(random_l[i])
#   ind_i <- PAM_indices_gpu(PAM = random_m, indices = "DF")
#   saveRDS(ind_i, paste0("PAM_indices/Null_Matrix/Herbs/", "r_ind", i, ".rds"))
#   return(ind_i) })
# 
# ####Get only shrubs####
# sp_shrub <- spinfo %>% filter(lifeForm == "Shrub")
# #Get PAM of shrubs
# #PAM <- fread("PAM.csv")
# spp_shrub <- subset(colnames(PAM), colnames(PAM) %in% sp_shrub$species)
# pam_shrub <- PAM %>% dplyr::select(spp_shrub) %>% as.matrix()
# 
# #Get indices from shrubs
# ind_shrubs <- PAM_indices_gpu(PAM = pam_shrub, indices = "all")
# #Save
# saveRDS(ind_shrubs, "PAM_indices/Indices_Shrubs.rds")
# 
# #Get NULL metrics of shrubs
# #Get random matrix in paralelle
# library(parallel)
# cl <- makeCluster(40)
# clusterExport(cl, varlist = "pam_shrub", #Get PAM2
#               envir=environment())
# clusterEvalQ(cl, {
#   library(picante)
#   library(pbapply)
# })
# 
# #dir.create("PAM_indices/Null_Matrix")
# dir.create("PAM_indices/Null_Matrix/Shrubs")
# pblapply(1:999, function(i){
#   set.seed(i)
#   random_matrix_i <- picante::randomizeMatrix(samp = pam_shrub, 
#                                               null.model = "trialswap",
#                                               iterations = length(pam_shrub)*10)
#   saveRDS(random_matrix_i, paste0("PAM_indices/Null_Matrix/Shrubs/m", i, ".rds"))
# }, cl = cl)
# stopCluster(cl)
# #get random matrix
# random_l <- list.files("PAM_indices/Null_Matrix/Shrubs/", full.names = T,
#                        pattern = "m")
# 
# #Calculate PAM indices from NULL matrix
# pblapply(seq_along(random_l), function(i){
#   random_m <- readRDS(random_l[i])
#   ind_i <- PAM_indices_gpu(PAM = random_m, indices = "DF")
#   saveRDS(ind_i, paste0("PAM_indices/Null_Matrix/Shrubs/", "r_ind", i, ".rds"))
#   return(ind_i) })
# 
# ####Get only subshrubs####
# sp_subshrub  <- spinfo %>% filter(lifeForm == "Subshrub")
# #Get PAM of subshrubs
# #PAM <- fread("PAM.csv")
# spp_subshrub  <- subset(colnames(PAM), colnames(PAM) %in% sp_subshrub$species)
# pam_subshrub  <- PAM %>% dplyr::select(spp_subshrub ) %>% as.matrix()
# 
# #Get indices from subshrubs
# ind_subshrubs <- PAM_indices_gpu(PAM = pam_subshrub , indices = "all")
# #Save
# saveRDS(ind_subshrubs, "PAM_indices/Indices_Subshrubs.rds")
# 
# #Get NULL metrics of subshrubs
# #Get random matrix in paralelle
# library(parallel)
# cl <- makeCluster(40)
# clusterExport(cl, varlist = "pam_subshrub", #Get PAM2
#               envir=environment())
# clusterEvalQ(cl, {
#   library(picante)
#   library(pbapply)
# })
# 
# #dir.create("PAM_indices/Null_Matrix")
# dir.create("PAM_indices/Null_Matrix/Subshrubs")
# pblapply(1:999, function(i){
#   set.seed(i)
#   random_matrix_i <- picante::randomizeMatrix(samp = pam_subshrub , 
#                                               null.model = "trialswap",
#                                               iterations = length(pam_subshrub )*10)
#   saveRDS(random_matrix_i, paste0("PAM_indices/Null_Matrix/Subshrubs/m", i, ".rds"))
# }, cl = cl)
# 
# #get random matrix
# random_l <- list.files("PAM_indices/Null_Matrix/Subshrubs/", full.names = T,
#                        pattern = "m")
# 
# #Calculate PAM indices from NULL matrix
# pblapply(seq_along(random_l), function(i){
#   random_m <- readRDS(random_l[i])
#   ind_i <- PAM_indices_gpu(PAM = random_m, indices = "DF")
#   saveRDS(ind_i, paste0("PAM_indices/Null_Matrix/Subshrubs/", "r_ind", i, ".rds"))
#   return(ind_i) })
# 
# 
# ####Wood species - Trees, shrubs and subshrubs####
# sp_wood  <- spinfo %>% filter(lifeForm == "Subshrub" |
#                                  lifeForm == "Shrub" | lifeForm == "Tree")
# unique(sp_wood$lifeForm)
# #Get PAM of woods
# #PAM <- fread("PAM.csv")
# spp_wood  <- subset(colnames(PAM), colnames(PAM) %in% sp_wood$species)
# pam_wood  <- PAM %>% dplyr::select(spp_wood ) %>% as.matrix()
# 
# #Get indices from woods
# ind_woods <- PAM_indices_gpu(PAM = pam_wood , indices = "all")
# #Save
# saveRDS(ind_woods, "PAM_indices/Indices_Woods.rds")
# 
# #Get NULL metrics of woods
# #Get random matrix in paralelle
# library(parallel)
# cl <- makeCluster(40)
# clusterExport(cl, varlist = "pam_wood", #Get PAM2
#               envir=environment())
# clusterEvalQ(cl, {
#   library(picante)
#   library(pbapply)
# })
# 
# #dir.create("PAM_indices/Null_Matrix")
# dir.create("PAM_indices/Null_Matrix/Woods")
# pblapply(1:999, function(i){
#   set.seed(i)
#   random_matrix_i <- picante::randomizeMatrix(samp = pam_wood , 
#                                               null.model = "trialswap",
#                                               iterations = length(pam_wood )*10)
#   saveRDS(random_matrix_i, paste0("PAM_indices/Null_Matrix/Woods/m", i, ".rds"))
# }, cl = cl)
# 
# #get random matrix
# random_l <- list.files("PAM_indices/Null_Matrix/Woods/", full.names = T,
#                        pattern = "m")
# 
# #Calculate PAM indices from NULL matrix
# pblapply(seq_along(random_l), function(i){
#   random_m <- readRDS(random_l[i])
#   ind_i <- PAM_indices_gpu(PAM = random_m, indices = "DF")
#   saveRDS(ind_i, paste0("PAM_indices/Null_Matrix/Woods/", "r_ind", i, ".rds"))
#   return(ind_i) })
# 
# ####Non_wood species - Trees, shrubs and subshrubs####
# sp_non_wood  <- spinfo %>% filter(lifeForm == "Herb" |
#                                  lifeForm == "Liana")
# unique(sp_non_wood$lifeForm)
# #Get PAM of non_woods
# #PAM <- fread("PAM.csv")
# spp_non_wood  <- subset(colnames(PAM), colnames(PAM) %in% sp_non_wood$species)
# pam_non_wood  <- PAM %>% dplyr::select(spp_non_wood ) %>% as.matrix()
# 
# #Get indices from non_woods
# ind_non_woods <- PAM_indices_gpu(PAM = pam_non_wood , indices = "all")
# #Save
# saveRDS(ind_non_woods, "PAM_indices/Indices_Non_woods.rds")
# 
# #Get NULL metrics of non_woods
# #Get random matrix in paralelle
# library(parallel)
# cl <- makeCluster(40)
# clusterExport(cl, varlist = "pam_non_wood", #Get PAM2
#               envir=environment())
# clusterEvalQ(cl, {
#   library(picante)
#   library(pbapply)
# })
# 
# #dir.create("PAM_indices/Null_Matrix")
# dir.create("PAM_indices/Null_Matrix/Non_woods")
# pblapply(1:999, function(i){
#   set.seed(i)
#   random_matrix_i <- picante::randomizeMatrix(samp = pam_non_wood , 
#                                               null.model = "trialswap",
#                                               iterations = length(pam_non_wood )*10)
#   saveRDS(random_matrix_i, paste0("PAM_indices/Null_Matrix/Non_woods/m", i, ".rds"))
# }, cl = cl)
# 
# #get random matrix
# random_l <- list.files("PAM_indices/Null_Matrix/Non_woods/", full.names = T,
#                        pattern = "m")
# 
# #Calculate PAM indices from NULL matrix
# pblapply(seq_along(random_l), function(i){
#   random_m <- readRDS(random_l[i])
#   ind_i <- PAM_indices_gpu(PAM = random_m, indices = "DF")
#   saveRDS(ind_i, paste0("PAM_indices/Null_Matrix/Non_woods/", "r_ind", i, ".rds"))
#   return(ind_i) })


####Get percentis and calculate dispersion field sign####
#Create folder to store results
dir.create("PAM_indices/Dispersion_sign")

#Import PAM
PAM <- fread("PAM.csv")

#By lifeform
lf <- list.dirs("PAM_indices/Null_Matrix/", full.names = F, recursive = F)
lf <- lf[2]
lf

# pblapply(seq_along(lf), function(i){
#   lf_i <- lf[i]
#   print(lf_i)
#   ind_i <- readRDS(paste0("PAM_indices/Indices_", lf_i, ".rds"))
#   ind_null_l <- list.files(path = paste0("PAM_indices/Null_Matrix/",
#                                          lf_i), 
#                            pattern = "r_ind", full.names = T)
#   ind_null_i <- lapply(seq_along(ind_null_l), function(i){
#     return(readRDS(ind_null_l[i]))
#   })
#   #Get dispersion field normalized by site and richness
#   ds <- (ind_i$Dispersion_field/ind_i$One_value_indices["Sites_Cells",])
#   ds <- ds/ind_i$One_value_indices["Species",]
#   ds[which(is.na(ds))] <- 0
#   
#   ds_null <- sapply(seq_along(ind_null_i), function(i){
#     ind_null_i_i <- ind_null_i[[i]]
#     n_sites <- ind_null_i_i$One_value_indices["Sites_Cells",]
#     S <- ind_null_i_i$One_value_indices["Species",]
#     ds_null_i <- ind_null_i_i$Dispersion_field/n_sites
#     ds_null_i <- ds_null_i/S
#     ds_null_i[which(is.na(ds_null_i))] <- 0
#     return(ds_null_i)
#   })
#   
#   ind_pos <- sapply(1:nrow(ds_null), function(i){
#     ds_i <- ds[i]
#     ds_null_i <- ds_null[i,]
#     pos5 <- quantile(ds_null_i, 0.05)
#     pos95 <- quantile(ds_null_i, 0.95)
#     pos_i <- ifelse(ds_i < pos5, -1,
#                     ifelse(ds_i > pos95, 1, 0))
#     return(pos_i)
#   })
#   
#   #Join with coordinates to plot
#   disp_sign <- PAM %>% as.data.frame() %>% dplyr::select(x, y) %>% 
#     mutate(DispersionSign = ind_pos)
#   #Join with dispersal field normalized and richness normalized
#   disp_sign <- disp_sign %>% 
#     mutate(NormalizedRichness = ind_i$Richness_normalized,
#            DispersedFieldNormalized = ds)
#   
#   # #Rasterize to see
#   # afr <- rast("Vetores/AF_raster.tiff")
#   # disp_sign_r <- rasterize(x = disp_sign %>% dplyr::select(x, y) %>% as.matrix(),
#   #                          y = afr,
#   #                          values = disp_sign$DispersionSign)
#   # 
#   # plot(disp_sign_r)
#   
#   #Salvar dataframe
#   write.csv(disp_sign,
#             paste0("PAM_indices/Dispersion_sign/", lf_i, ".csv"),
#             row.names = F)
#   
# })
# 
# 




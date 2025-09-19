#### INLA script 2 ####
library(dplyr)
library(pbapply)
library(terra)
library(Matrix)
library(INLA)
library(data.table)
library(hier.part)

# Create directory to save models
dir.create("Data/INLA_models/Models", recursive = TRUE)

# Configure formulas for Partition importance
# Get vars
covs <- c("Aridity","PET","Bio06","Bio02","Prec_stab","Temp_stab","Topo_het")

#Create combination
ragged_comb <- combos(length(covs))$ragged
combs <- sapply(1:nrow(ragged_comb), function(i){
  comb_i <- as.vector(ragged_comb[i, ][ragged_comb[i, ] > 0])
  current_comb <- paste(covs[comb_i], collapse = " + ")
  current_comb <- paste0("y ~ -1 + intercept + ", current_comb, " + f(spatial, model = spde)")
})
#Add null model in the first line
combs <- c("y ~ -1 + f(spatial, model = spde)", combs)


#### Analyzes with modeles species ####
#Import data
all_var <- rast("Data/Variables/Explanatory_Variables.tiff")
my_var <- names(all_var)

#Get richness of lifeforms
lf_indices <- readRDS("Data/Richness_by_lifeform.rds")
# Get xy
xy <- lf_indices$xy
lf_indices$xy <- NULL

# # Merge data
# d <- lapply(names(lf_indices), function(x){
#   data.frame(lifeform = x, Richness = lf_indices[[x]])
# }) %>% rbindlist()

#Extract values of environmental variables 
all_data <- pblapply(names(lf_indices), function(i){
  #Create dataframe with lifeform, coordinates and richness
  df_i <- data.frame(xy,
                     lifeform = i,
                     Richness = lf_indices[[i]])
  #Extract values from explanatory variabkes
  pt_var <- terra::extract(all_var, as.matrix(xy)) %>%
    cbind(df_i, .) %>% na.omit()
  return(pt_var)
})

#Rename
names(all_data) <- names(lf_indices)

# Looping in lifeforms
lf <- c("All", "Tree", "Liana", "Shrub", "Subshrub", "Terrestrial_herb", 
        "Epiphytic_herb")
# To test
i <- "Tree"


res <- pblapply(lf, function(i){
  
  message("Fitting model for ", i)
  
  # Get data
  d_i <- all_data[[i]]
  
  # Remove pixels with less than 1 species
  d_i <- d_i %>% filter(Richness >= 1)
  
  
  # -------------------------------
  # 1?????? Preparar dados
  # -------------------------------
  covs <- c("Aridity","PET","Bio06","Bio02","Prec_stab","Temp_stab","Topo_het")
  d_i[covs] <- scale(d_i[covs])  # padronizar covariáveis
  coords <- as.matrix(d_i[,c("x","y")])
  
  # -------------------------------
  # 2?????? Definir mesh automaticamente
  # -------------------------------
  # Amostra para estimar distância média (economiza tempo)
  set.seed(42)
  sample_idx <- sample(1:nrow(coords), min(1000, nrow(d_i)))
  coords_sample <- coords[sample_idx, ]
  dists_sample <- as.matrix(dist(coords_sample))
  mean_dist <- mean(dists_sample[upper.tri(dists_sample)])
  
  # Extensão da área
  x_range <- max(coords[,1]) - min(coords[,1])
  y_range <- max(coords[,2]) - min(coords[,2])
  area_range <- max(x_range, y_range)
  
  # Parâmetros automáticos
  cutoff_val <- 0.2 * mean_dist #0.05
  min_edge <- area_range / 20 #20
  max_edge <- area_range / 5 #5
  
  # cat("cutoff =", cutoff_val, "\n")
  # cat("max.edge = c(", min_edge, ",", max_edge, ")\n")
  
  # Criar mesh
  mesh <- inla.mesh.2d(loc = coords,
                       max.edge = c(min_edge, max_edge),
                       cutoff = cutoff_val)
  # plot(mesh)
  # points(coords, col="red", cex=0.3)
  
  # -------------------------------
  # 3?????? Definir SPDE
  # -------------------------------
  spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
  s.index <- inla.spde.make.index("spatial", n.spde = spde$n.spde)
  A <- inla.spde.make.A(mesh = mesh, loc = coords)
  
  # -------------------------------
  # 4?????? Criar stack
  # -------------------------------
  stk <- inla.stack(
    data = list(y = d_i$Richness),
    A = list(A, 1),
    effects = list(
      spatial = s.index$spatial,
      data.frame(intercept = 1, d_i[, covs])
    ),
    tag = "est"
  )
  
  # Set formula
  # 6. Fórmula
  form <- y ~ -1 + intercept + Aridity + PET + Bio06 + Bio02 + Prec_stab + Temp_stab + Topo_het +
    f(spatial, model = spde)
  
  # Fit model
  m <- inla(form,
              family = "nbinomial",
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
              control.compute = list(dic = TRUE, waic = TRUE),
              verbose = FALSE)
  
  message("Calculating importance...")
  
  comb_m <- pblapply(combs, function(z){
    # message("Fitting ", z)
    m_x <- inla(as.formula(z),
                family = "nbinomial",
                data = inla.stack.data(stk),
                control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
                control.compute = list(dic = TRUE, waic = TRUE),
                verbose = FALSE)
    
    #mlik
    mlik <- m_x$mlik[1]
    
  #Create dataframe
    df_z <- data.frame(variable.combination = z,
                       mlik = mlik)
  })
  
  #Join results in a unique dataframe
  comb_df <- bind_rows(comb_m)
  
  #Partition
  var_p <- partition(comb_df$mlik, pcan = length(covs),
                     var.names = covs)
  
  var_imp <- var_p$I.perc %>% mutate(Variable = row.names(.), .before = 1) %>%
    rename(Importance = ind.exp.var)
  var_imp <- var_imp %>%
    mutate(lifeForm = i, .before = 1)
 
  return(list("data" = d_i,
              "model" = m,
              "stk" = stk,
              "mesh_matrice" = A,
              "var_imp" = var_imp))
})
names(res) <- lf

# Salvar resultados
saveRDS(res, "Data/INLA_models/models.rds")

# Fast check of results
r <- pblapply(names(res), function(x){
  # For test
  # x <- res[[1]]
  # Efeitos fixos
  res[[x]]$model$summary.fixed %>% 
    mutate(Variable = row.names(.), lifeform = x,
           .before = 1)
}) %>% rbindlist()

# Residuals
residuos <- pblapply(names(res), function(x){
  # For test
  # x <- res[[1]]
  
  # Get data
  # Get data
  d_i <- all_data[[x]]
  
  # Remove pixels with less than 1 species
  d_i <- d_i %>% filter(Richness >= 1)
  
  # Check residuals
  index <- inla.stack.index(res[[x]]$stk, tag="est")$data
  d_i$predicted <- res[[x]]$model$summary.fitted.values[index, "mean"]
  d_i$residuals <- d_i$predicted - d_i$Richness
  r_residuals <- rasterize(as.matrix(d_i[,c("x", "y")]),
                           all_var, values = d_i$residuals)
  plot(r_residuals, main = x)
  # Check moran
  moran <- moranfast::moranfast(d_i$residuals, d_i$x, d_i$y)
  moran
})
residuos


importance <- pblapply(names(res), function(x){
  res[[x]]$var_imp
}) %>% rbindlist()


#### Analyzes including rare and undersamples species ####
#Import data
all_var <- rast("Data/Variables/Explanatory_Variables.tiff")
my_var <- names(all_var)

#Get richness of lifeforms
lf_indices_wo <- readRDS("Data/Richness_by_lifeform_withOccurrences.rds")
# Get xy
xy_wo <- lf_indices_wo$xy
lf_indices_wo$xy <- NULL


#Extract values of environmental variables 
all_data_wo <- pblapply(names(lf_indices_wo), function(i){
  #Create dataframe with lifeform, coordinates and richness
  df_i <- data.frame(xy_wo,
                     lifeform = i,
                     Richness = lf_indices_wo[[i]])
  #Extract values from explanatory variabkes
  pt_var <- terra::extract(all_var, as.matrix(xy_wo)) %>%
    cbind(df_i, .) %>% na.omit()
  return(pt_var)
})

#Rename
names(all_data_wo) <- names(lf_indices_wo)

# Looping in lifeforms
lf <- c("All", "Tree", "Liana", "Shrub", "Subshrub", "Terrestrial_herb", 
        "Epiphytic_herb")
# To test
i <- "Tree"


res_wo <- pblapply(lf, function(i){
  
  message("Fitting model for ", i)
  
  # Get data
  d_i <- all_data_wo[[i]]
  
  # Remove pixels with less than 1 species
  d_i <- d_i %>% filter(Richness >= 1)
  
  
  # -------------------------------
  # 1?????? Preparar dados
  # -------------------------------
  covs <- c("Aridity","PET","Bio06","Bio02","Prec_stab","Temp_stab","Topo_het")
  d_i[covs] <- scale(d_i[covs])  # padronizar covariáveis
  coords <- as.matrix(d_i[,c("x","y")])
  
  # -------------------------------
  # 2?????? Definir mesh automaticamente
  # -------------------------------
  # Amostra para estimar distância média (economiza tempo)
  set.seed(42)
  sample_idx <- sample(1:nrow(coords), min(1000, nrow(d_i)))
  coords_sample <- coords[sample_idx, ]
  dists_sample <- as.matrix(dist(coords_sample))
  mean_dist <- mean(dists_sample[upper.tri(dists_sample)])
  
  # Extensão da área
  x_range <- max(coords[,1]) - min(coords[,1])
  y_range <- max(coords[,2]) - min(coords[,2])
  area_range <- max(x_range, y_range)
  
  # Parâmetros automáticos
  cutoff_val <- 0.2 * mean_dist #0.05
  min_edge <- area_range / 20 #20
  max_edge <- area_range / 5 #5
  
  # cat("cutoff =", cutoff_val, "\n")
  # cat("max.edge = c(", min_edge, ",", max_edge, ")\n")
  
  # Criar mesh
  mesh <- inla.mesh.2d(loc = coords,
                       max.edge = c(min_edge, max_edge),
                       cutoff = cutoff_val)
  # plot(mesh)
  # points(coords, col="red", cex=0.3)
  
  # -------------------------------
  # 3?????? Definir SPDE
  # -------------------------------
  spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
  s.index <- inla.spde.make.index("spatial", n.spde = spde$n.spde)
  A <- inla.spde.make.A(mesh = mesh, loc = coords)
  
  # -------------------------------
  # 4?????? Criar stack
  # -------------------------------
  stk <- inla.stack(
    data = list(y = d_i$Richness),
    A = list(A, 1),
    effects = list(
      spatial = s.index$spatial,
      data.frame(intercept = 1, d_i[, covs])
    ),
    tag = "est"
  )
  
  # Set formula
  # 6. Fórmula
  form <- y ~ -1 + intercept + Aridity + PET + Bio06 + Bio02 + Prec_stab + Temp_stab + Topo_het +
    f(spatial, model = spde)
  
  # Fit model
  m <- inla(form,
            family = "nbinomial",
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE),
            verbose = FALSE)
  
  message("Calculating importance...")
  
  comb_m <- pblapply(combs, function(z){
    # message("Fitting ", z)
    m_x <- inla(as.formula(z),
                family = "nbinomial",
                data = inla.stack.data(stk),
                control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
                control.compute = list(dic = TRUE, waic = TRUE),
                verbose = FALSE)
    
    #mlik
    mlik <- m_x$mlik[1]
    
    #Create dataframe
    df_z <- data.frame(variable.combination = z,
                       mlik = mlik)
  })
  
  #Join results in a unique dataframe
  comb_df <- bind_rows(comb_m)
  
  #Partition
  var_p <- partition(comb_df$mlik, pcan = length(covs),
                     var.names = covs)
  
  var_imp <- var_p$I.perc %>% mutate(Variable = row.names(.), .before = 1) %>%
    rename(Importance = ind.exp.var)
  var_imp <- var_imp %>%
    mutate(lifeForm = i, .before = 1)
  
  return(list("data" = d_i,
              "model" = m,
              "stk" = stk,
              "mesh_matrice" = A,
              "var_imp" = var_imp))
})
names(res_wo) <- lf

# Salvar resultados
saveRDS(res_wo, "Data/INLA_models/models_with_Ocurrences.rds")

# Fast check of results
r_wo <- pblapply(names(res_wo), function(x){
  # For test
  # x <- res[[1]]
  # Efeitos fixos
  res_wo[[x]]$model$summary.fixed %>% 
    mutate(Variable = row.names(.), lifeform = x,
           .before = 1)
}) %>% rbindlist()

# Residuals
residuos_wo <- pblapply(names(res_wo), function(x){
  # For test
  # x <- res[[1]]
  
  # Get data
  # Get data
  d_i <- res_wo[[x]]$data
  
  # Remove pixels with less than 1 species
  d_i <- d_i %>% filter(Richness >= 1)
  
  # Check residuals
  index <- inla.stack.index(res_wo[[x]]$stk, tag="est")$data
  d_i$predicted <- res_wo[[x]]$model$summary.fitted.values[index, "mean"]
  d_i$residuals <- d_i$predicted - d_i$Richness
  r_residuals <- rasterize(as.matrix(d_i[,c("x", "y")]),
                           all_var, values = d_i$residuals)
  plot(r_residuals, main = x)
  # Check moran
  moran <- moranfast::moranfast(d_i$residuals, d_i$x, d_i$y)
  moran
})
residuos_wo


importance <- pblapply(names(res_wo), function(x){
  res_wo[[x]]$var_imp
}) %>% rbindlist()



### Try to plot response curve
# Use something like kuenm2


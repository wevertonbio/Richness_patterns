#### INLA script 2 ####
library(dplyr)
library(pbapply)
library(terra)
library(Matrix)
library(INLA)
library(data.table)

# Create directory to save models
dir.create("Data/INLA_models/Models", recursive = TRUE)


#Import data
all_var <- rast("Data/Variables/Explanatory_Variables.tiff")
my_var <- names(all_var)

#Get richness of lifeforms
lf_indices <- readRDS("Data/Richness_by_lifeform.rds")
# Get xy
xy <- lf_indices$xy
lf_indices$xy <- NULL

# Merge data
d <- lapply(names(lf_indices), function(x){
  data.frame(lifeform = x, Richness = lf_indices[[x]])
}) %>% rbindlist()



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


# Lista de covariáveis contínuas
covs <- c("Aridity","PET","Bio06","Bio02","Prec_stab","Temp_stab","Topo_het")

# Gerar todas as combinações com pelo menos 5 variáveis
all_combs <- unlist(lapply(5:length(covs), function(n) combn(covs, n, simplify = FALSE)), recursive = FALSE)

# Conferir quantas combinações existem
length(all_combs)  # Para saber quantos modelos vamos testar

# Exemplo: mostrar primeiras 5 combinações
all_combs[1:5]

# Para cada combinação, gerar a fórmula com SPDE
formulas <- lapply(all_combs, function(vars){
  as.formula(
    paste0("y ~ -1 + intercept + ",
           paste(vars, collapse = " + "),
           " + f(spatial, model = spde)")
  )
})

formulas_text <- lapply(all_combs, function(vars){
  paste0("y ~ -1 + intercept + ",
         paste(vars, collapse = " + "),
         " + f(spatial, model = spde)")
})

# Looping in lifeforms
lf <- names(all_data)
# To test
i <- "Terrestrial_herb"

pblapply(lf, function(i){
  # Get data
  d_i <- all_data[[i]]
  
  
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
  
  cat("cutoff =", cutoff_val, "\n")
  cat("max.edge = c(", min_edge, ",", max_edge, ")\n")
  
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
  
  # Testar modelos
  res <- pblapply(formulas, function(x){
    res_x <- inla(x,
         family = "nbinomial",
         data = inla.stack.data(stk),
         control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
         control.compute = list(dic = TRUE, waic = TRUE),
         verbose = FALSE)
    return(res_x)
  })
  
  # To test
  # bm <- res_x
  # Get dic and waic
  dic_values <- sapply(res, function(x) x$dic$dic)
  waic_values <- sapply(res, function(x) x$waic$waic)
  # Create dataframe
  df <- data.frame("formulas" = unlist(formulas_text),
                   dic = dic_values,
                   waic = waic_values)
  # Gest best model
  bm <- which(df$dic == min(df$dic))
  bm <- res[[bm]]
  bm$summary.fixed
  
  # p <- autoplot(bm)
  # p$fixed.marginals
  # p$hyper.marginals
  # p$random.effects.line
  # p$marginal.fitted
  
  # Efeitos fixos
  fixed <- bm$summary.fixed
  fixed$Variable <- rownames(fixed)
  
  # Plotando coeficientes com IC 95%
  ggplot(fixed[-1, ], aes(x = Variable, y = mean)) +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(y = "Estimativa", x = "Covariáveis", title = "Efeitos fixos do modelo INLA")
  
  # Check residuals
  index <- inla.stack.index(stk, tag="est")$data
  d_i$predicted <- bm$summary.fitted.values[index, "mean"]
  d_i$residuals <- d_i$predicted - d_i$Richness
  r_residuals <- rasterize(as.matrix(d_i[,c("x", "y")]), 
                           all_var, values = d_i$residuals)
  
  ggplot() +
    geom_spatraster(data = r_residuals) +
    scale_fill_whitebox_c()
  
  # Check moran
  moran <- moranfast::moranfast(d_i$residuals, d_i$x, d_i$y)
  moran
  
})



### Try to plot response curve
# Use something like kuenm2

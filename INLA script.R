#### Fit best models and partitioning the variance ####
#We didn't scale the variables. See: https://stackoverflow.com/questions/25481331/how-to-use-scale-in-logistic-regression-correctly
library(MASS)
library(dplyr)
library(pbapply)
library(terra)
library(sjPlot)
library(performance)
library(hier.part)
library(parallel)
library(Matrix)
library(INLA)

#Create directory to save data to plot
dir.create("Data/Models/Partitioning", recursive = T)
dir.create("Data/Models/Predictions")
dir.create("Data/Models/Best_models")


#Import data
my_var <- names(all_var)

#Get richness of lifeforms
lf_indices <- readRDS("Data/Richness_by_lifeform.rds")
# Get xy
xy <- lf_indices$xy
lf_indices$xy <- NULL

#Plot rasters to check
rich <- rast(pblapply(lf_indices, function(x){
  r_x <- rasterize(x = as.matrix(xy),
                   y = all_var, values = x, fun = "mean",
                   background = 0)
}))
plot(rich[[1]])

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

#Import candidate models
cm_l <- list.files("Data/Models/Candidate_models//", pattern = ".RDS",
                   full.names = TRUE, recursive = FALSE)

#Remove bamboos and palms
cm_l <- cm_l[!grepl("Bamboo|Palm", cm_l)]

#Generate data to expand grid with NA values for variables that are not selected
names_var <- names(all_var)
#Remove spatial predictors
names_var <- names_var[!grepl("EV", names_var)]

#To test
i <- cm_l[7]
i

cm_i <- readRDS(i) #Read i data
#Get lifeform and endemism
lf_i <- unique(cm_i$lifeForm)
message("Getting results for ", lf_i)

#Select data
d_i <- all_data[[lf_i]] %>% filter(Richness > 20)

# Salvar d_i
dir.create("INLA_teste")
write.csv(d_i, "INLA_teste/trees.csv")


# 1. padronizar covariáveis
covs <- c("Aridity","PET","Bio06","Bio02","Bio15","Prec_stab","Temp_stab","Topo_het")
d_i[covs] <- scale(d_i[covs])

# 2. mesh (controle do tamanho via max.edge)
coords <- as.matrix(d_i[,c("x","y")])
mesh <- inla.mesh.2d(loc=coords,
                     max.edge = c(0.5, 2),   # ajustar: valores em graus; menor -> mais nós -> mais custo
                     cutoff = 0.01)          # evita multi-pontos muito próximos

plot(mesh); points(coords, col="red", cex=0.3)

# 3. SPDE object
spde <- inla.spde2.matern(mesh = mesh, alpha = 2)

# 4. Índices e projeção
A <- inla.spde.make.A(mesh = mesh, loc = coords)
s.index <- inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)

# 5. Stack
stk <- inla.stack(
  data = list(y = d_i$Richness),
  A = list(A, 1),
  effects = list(spatial = s.index$spatial, data.frame(intercept = 1, d_i[,covs])),
  tag = "est"
)

# 6. Fórmula
form <- y ~ -1 + intercept + Aridity + PET + Bio06 + Bio02 + Bio15 +
  Prec_stab + Temp_stab + Topo_het +
  f(spatial, model = spde)


# 7. Rodar INLA (com NegBinomial)
inla.setOption(num.threads = 8)   # ajuste conforme cpu
res_spde_nb <- inla(form,
                    family = "nbinomial",
                    data = inla.stack.data(stk),
                    control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
                    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                    verbose = TRUE)

summary(res_spde_nb)
summary(res_spde_nb)$fixed

res_spde_nb$residuals

spatial_field <- inla.spde2.result(res_spde_nb, "spatial", spde)
plot(spatial_field)


# índice dentro do stack
index <- inla.stack.index(stk, tag="est")$data

# valores preditos de Richness
d_i$Richness_pred <- res_spde_nb$summary.fitted.values[index, "mean"]
d_i %>% select(Richness, Richness_pred) %>% View()
d_i$residuals <- d_i$Richness_pred - d_i$Richness

#Check spatial autocorrelation
residuals_r <- rasterize(d_i %>% select(x, y) %>% as.matrix(),
                         all_var, values = d_i$residuals)

plot(residuals_r)

ggplot(d_i, aes(x, y, colour = residuals)) +
  scale_color_gradient2() +
  geom_point(size = 0.5)


# campo espacial latente (residual espacial)
d_i$spatial_effect <- res_spde_nb$summary.random$spatial$mean

plot_model(res_spde_nb)

# 100 valores igualmente espaçados entre min e max
aridity_seq <- seq(min(d_i$Aridity), max(d_i$Aridity), length.out = 100)

# Criar dataframe com valores médios das outras covariáveis
cov_means <- as.list(colMeans(d_i[, covs[-which(covs=="Aridity")]]))
df_grid <- data.frame(intercept = 1,
                      Aridity = aridity_seq,
                      PET = cov_means$PET,
                      Bio06 = cov_means$Bio06,
                      Bio02 = cov_means$Bio02,
                      Bio15 = cov_means$Bio15,
                      Prec_stab = cov_means$Prec_stab,
                      Temp_stab = cov_means$Temp_stab,
                      Topo_het = cov_means$Topo_het)
A_grid <- inla.spde.make.A(mesh = mesh, loc = as.matrix(d_i[1, c("x","y")])) 

# Fitted = X %*% beta (efeito fixo)
beta <- res_spde_nb$summary.fixed$mean
X_grid <- as.matrix(df_grid)
pred_fixed <- X_grid %*% beta
pred_fixed <- exp(pred_fixed)  # NegBinomial usa log-link

Sigma <- res_spde_nb$cov.fixed
pred_sd <- sqrt(diag(X_grid %*% Sigma %*% t(X_grid)))

# Intervalos 95% (em log-scale)
pred_lower <- exp(pred_fixed - 1.96 * pred_sd)
pred_upper <- exp(pred_fixed + 1.96 * pred_sd)

library(ggplot2)

# Efeitos fixos
fixed <- res_spde_nb$summary.fixed
fixed$Variable <- rownames(fixed)

# Plotando coeficientes com IC 95%
ggplot(fixed[-1, ], aes(x = Variable, y = mean)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(y = "Estimativa", x = "Covariáveis", title = "Efeitos fixos do modelo INLA")

# Ver autocorrelação
moran <- moranfast::moranfast((d_i$residuals), d_i$x, d_i$y)
moran


#### Novo script ####
# -------------------------------
# 0?????? Pacotes
# -------------------------------
library(INLA)
library(sp)
library(ggplot2)
library(dplyr)
library(sf)
library(ggeffects)   # Para curvas marginais

# Impotar dados
d_i <- read.csv("INLA_teste/trees.csv")

# -------------------------------
# 1?????? Preparar dados
# -------------------------------
# Supondo que d_i já existe com colunas: x, y, Richness, Aridity, PET, Bio06, Bio02, Bio15, Prec_stab, Temp_stab, Topo_het
covs <- c("Aridity","PET","Bio06","Bio02","Prec_stab","Temp_stab","Topo_het")
d_i[covs] <- scale(d_i[covs])  # padronizar covariáveis
coords <- as.matrix(d_i[,c("x","y")])

# -------------------------------
# 2?????? Definir mesh automaticamente
# -------------------------------
# Amostra para estimar distância média (economiza tempo)
set.seed(123)
sample_idx <- sample(1:nrow(d_i), min(1000, nrow(d_i)))
coords_sample <- coords[sample_idx, ]
dists_sample <- as.matrix(dist(coords_sample))
mean_dist <- mean(dists_sample[upper.tri(dists_sample)])

# Extensão da área
x_range <- max(coords[,1]) - min(coords[,1])
y_range <- max(coords[,2]) - min(coords[,2])
area_range <- max(x_range, y_range)

# Parâmetros automáticos
cutoff_val <- 0.05 * mean_dist
min_edge <- area_range / 20
max_edge <- area_range / 5

cat("cutoff =", cutoff_val, "\n")
cat("max.edge = c(", min_edge, ",", max_edge, ")\n")

# Criar mesh
mesh <- inla.mesh.2d(loc = coords,
                     max.edge = c(min_edge, max_edge),
                     cutoff = cutoff_val)
plot(mesh)
points(coords, col="red", cex=0.3)

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

# -------------------------------
# 5?????? Rodar modelo NegBinomial
# -------------------------------
form <- y ~ -1 + intercept + Aridity + PET + Bio06 + Bio02 + 
  Prec_stab + Temp_stab + Topo_het +
  f(spatial, model = spde)


form <- y ~ -1 + intercept +
  Aridity + 
  PET + I(PET^2) +
  Bio06 + 
  Bio02 + 
  Prec_stab + 
  Temp_stab + 
  Topo_het + 
  f(spatial, model = spde)

res <- inla(form,
            family = "nbinomial",
            data = inla.stack.data(stk),
            control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
            verbose = TRUE)

# -------------------------------
# 6?????? Extrair predições
# -------------------------------
index <- inla.stack.index(stk, tag="est")$data
d_i$Richness_pred <- res$summary.fitted.values[index, "mean"]
# A matriz A projeta o campo SPDE nos pontos originais
spatial_mean <- as.vector(A %*% res$summary.random$spatial$mean)

# Agora adiciona ao dataframe
d_i$spatial_effect <- spatial_mean


# -------------------------------
# 7?????? Transformar em sf para plot
# -------------------------------
d_sf <- sf::st_as_sf(d_i, coords = c("x","y"), crs = 4326)

# -------------------------------
# 8?????? Mapas com ggplot2
# -------------------------------
# Riqueza prevista
ggplot(d_sf) +
  geom_sf(aes(color = Richness_pred), size = 1.5) +
  scale_color_viridis_c(option = "C") +
  theme_minimal() +
  labs(title = "Riqueza prevista de espécies", color = "Riqueza")

# Efeito espacial residual
ggplot(d_sf) +
  geom_sf(aes(color = spatial_effect), size = 1.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Efeito espacial residual", color = "Residual")

# Ver residuos 2
d_i$residuals2 <- d_i$Richness_pred - d_i$Richness
ggplot(d_i, aes(x, y, colour = residuals2)) +
  scale_color_gradient2() +
  geom_point(size = 0.5)



# -------------------------------
# 9?????? Curvas marginais das covariáveis (fixas)
# -------------------------------
for(var in covs){
  pred <- ggpredict(res, terms = "Aridity")
  plot(pred) + ggtitle(paste("Curva marginal:", var))
}

# -------------------------------
# 10?????? Resumo rápido do modelo
# -------------------------------
summary(res$summary.fixed)
summary(res$summary.hyperpar)
res$dic$dic
res$waic$waic

# Efeitos fixos
# Efeitos fixos
fixed <- res$summary.fixed
fixed$Variable <- rownames(fixed)

# Plotando coeficientes com IC 95%
ggplot(fixed[-1, ], aes(x = Variable, y = mean)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(y = "Estimativa", x = "Covariáveis", title = "Efeitos fixos do modelo INLA")


# Grid de Bio06
bio6 <- seq(min(d_i$Bio06), max(d_i$Bio06), length.out = 100)

# Médias das outras covariáveis
cov_means <- colMeans(d_i[, covs[-which(covs=="Bio06")]])

# Criar dataframe incluindo termo quadrático de Bio06
df_grid <- data.frame(
  intercept = 1,
  Aridity = cov_means["Aridity"],
  PET = cov_means["PET"],
  Bio06 = bio6,
  `I(Bio06^2)` = bio6^2,   # termo quadrático
  Bio02 = cov_means["Bio02"],
  Bio15 = cov_means["Bio15"],
  Prec_stab = cov_means["Prec_stab"],
  Temp_stab = cov_means["Temp_stab"],
  Topo_het = cov_means["Topo_het"]
)

# Obter nomes dos coeficientes do modelo
coef_names <- names(res$summary.fixed)

# Criar dataframe de predição com colunas correspondentes
df_grid <- data.frame(matrix(NA, nrow=100, ncol=length(coef_names)))
colnames(df_grid) <- coef_names

# Grid de Bio06
bio6 <- seq(min(d_i$Bio06), max(d_i$Bio06), length.out = 100)

# Preencher valores
df_grid$intercept <- 1
df_grid$Aridity <- mean(d_i$Aridity)
df_grid$`I(Aridity^2)` <- mean(d_i$Aridity)^2
df_grid$PET <- mean(d_i$PET)
df_grid$`I(PET^2)` <- mean(d_i$PET)^2
df_grid$Bio06 <- bio6
df_grid$`I(Bio06^2)` <- bio6^2
df_grid$Bio02 <- mean(d_i$Bio02)
df_grid$`I(Bio02^2)` <- mean(d_i$Bio02)^2
df_grid$Bio15 <- mean(d_i$Bio15)
df_grid$`I(Bio15^2)` <- mean(d_i$Bio15)^2
df_grid$Prec_stab <- mean(d_i$Prec_stab)
df_grid$`I(Prec_stab^2)` <- mean(d_i$Prec_stab)^2
df_grid$Temp_stab <- mean(d_i$Temp_stab)
df_grid$`I(Temp_stab^2)` <- mean(d_i$Temp_stab)^2
df_grid$Topo_het <- mean(d_i$Topo_het)
df_grid$`I(Topo_het^2)` <- mean(d_i$Topo_het)^2


X_grid <- as.matrix(df_grid)
beta <- res$summary.fixed$mean

pred_fixed <- exp(X_grid %*% beta)
plot(bio6, pred_fixed, type="l", 
     main="Curva marginal Bio06", 
     ylab="Riqueza prevista", xlab="Bio06")

# Predição na escala da resposta (NegBinomial)
pred_fixed <- exp(X_grid %*% beta)

# Plot
plot(bio6, pred_fixed, type="l", 
     main="Curva marginal Bio06", 
     ylab="Riqueza prevista", xlab="Bio06")


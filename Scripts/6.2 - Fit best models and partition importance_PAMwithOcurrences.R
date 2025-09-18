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

#Create directory to save data to plot
dir.create("Data/Models_withOccurrences//Partitioning", recursive = T)
dir.create("Data/Models_withOccurrences/Predictions")
dir.create("Data/Models_withOccurrences/Best_models")


#Import data
my_f <- readRDS("Data/Variables/Formulas.RDS")
my_f %>% data.frame() %>% View()
all_var <- rast("Data/Variables/Explanatory_Variables.tiff")
my_var <- names(all_var)

#Get richness of lifeforms
lf_indices <- readRDS("Data/Richness_by_lifeform_withOccurrences.rds")
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
cm_l <- list.files("Data/Models_withOccurrences//Candidate_models//", pattern = ".RDS",
                   full.names = TRUE, recursive = FALSE)

#Remove bamboos and palms
cm_l <- cm_l[!grepl("Bamboo|Palm", cm_l)]

#Generate data to expand grid with NA values for variables that are not selected
names_var <- names(all_var)
#Remove spatial predictors
names_var <- names_var[!grepl("EV", names_var)]

#To test
i <- cm_l[6]
i

pblapply(cm_l, function(i){
  cm_i <- readRDS(i) #Read i data
  #Get lifeform and endemism
  lf_i <- unique(cm_i$lifeForm)
  message("Getting results for ", lf_i)
  
  #Select data
  d_i <- all_data[[lf_i]] %>% filter(Richness > 20)
  
  #Select best model
  bm <- cm_i %>%
    # Remove MID Domain
    filter(!grepl("Mid_domain", Formula)) %>%
    filter(!grepl("Bio15", Formula)) %>%
    # Remove overdispersed models
    filter(p_dispersion_ration > 0.05) %>%
    #calculate delta AIC
    mutate(dAIC = AIC - min(AIC, na.rm = T), .after = AIC) %>%
    filter(dAIC <= 2) #AIC < 2
  bm$Formula
  
  if(nrow(bm) == 0){
    bm <- cm_i %>%
      # Remove MID Domain
      filter(!grepl("Mid_domain", Formula)) %>%
      filter(!grepl("Bio15", Formula)) %>%
      # Remove overdispersed models
      # filter(p_dispersion_ration > 0.05) %>%
      #calculate delta AIC
      mutate(dAIC = AIC - min(AIC, na.rm = T), .after = AIC) %>%
      filter(dAIC <= 2) #AIC < 2
  }
  
  if(nrow(bm) > 1){ #Get more comples model
    n_var <- sapply(1:nrow(bm), function(i){
      bm$Formula[i] %>% as.formula() %>% all.vars() %>% length()
    })
    to_keep <- n_var == max(n_var)
    bm <- bm[to_keep,]
  }
  
  #Get best formulas
  bf <- bm$Formula
  # bf <- "Richness ~ Aridity + PET + Bio06 + Bio02 + Prec_stab + Temp_stab + Topo_het + Lat_EV + Long_EV + I(Lat_EV^2) + I(Long_EV^2)"
  #Fit model
  m_i <- try(glm.nb(formula = bf,
                    data = d_i), silent = T)
  # Get summary
  s <- summary(m_i)$coefficients %>% as.data.frame() %>% 
    mutate(Variable = row.names(.), .before = 1)
  
  
  # plot_model(m_i, type = "pred")
  
  # performance::check_overdispersion(m_i)
  # 
  # # sr2 <- DHARMa::simulateResiduals(m_i)
  # # plot(sr2)
  #  p <- performance::check_model(m_i, check = c("outliers", "pp_check", 
  #                                               "overdispersion", "normality"))
  # # p
  # # 
  # 
  # simulated_residuals <- simulate_residuals(m_i)
  # check_residuals(simulated_residuals)
  # check_overdispersion(simulated_residuals)
  # check_outliers(simulated_residuals, type = 'bootstrap')
  # check_model(simulated_residuals, dot_size = 1.5)
  
  #summary(m_i)
  #plot_model(m_i, type = "pred")
  
  ####Variable importance####
  #Get variables
  vars <- gsub("Richness ~ ", "", bf) %>% strsplit(., " \\+ ") %>% unlist()
  
  # Get unique vars
  unique_vars <- gsub("I\\(|\\^2\\)", "", vars) %>% strsplit(., ":") %>%
    unlist() %>% unique()
  unique_vars
  
  # Group linear and quadratic variables
  agrupar_variaveis <- function(var) {
    if (grepl("^I\\(", var) & !grepl("EV|Mid_domain|Topoi_het", var)) {
      # If there is quadratic, find linear version of the variable
      nome_linear <- sub("^I\\(([^\\^]+)\\^2\\)", "\\1", var)
      # Check if linear is in the original vector
      if (nome_linear %in% vars) {
        return(paste(nome_linear, "+", var))
      } else {
        return(NULL)
      }
    }
    
    if(!grepl("^I\\(", var) & !grepl("EV|Mid_domain|Topoi_het", var)){
      # Check if there is quadratic version of the linear predictor
      var_quadratica <- paste("I(", var, "^2)", sep = "")
      if (var_quadratica %in% vars) {
        return(paste(var, "+", var_quadratica))
      } else {return(var)}
    }
    
    if(grepl("EV|Mid_domain|Topoi_het", var)) {
      return(var) }
  }
  
  # Appply function to group variable
  vars_agrupadas <- lapply(vars, agrupar_variaveis)
  vars_agrupadas <- Filter(Negate(is.null), vars_agrupadas) %>% unlist() %>% unique()
  
  # Group spatial variables
  vars_no_spatial <- vars_agrupadas[!grepl("EV", vars_agrupadas)]
  vars_spatial <- paste(vars_agrupadas[grepl("EV", vars_agrupadas)],
                        collapse = " + ")
  comb_formulas <- c(vars_no_spatial, vars_spatial)
  
  #Create combination
  ragged_comb <- combos(length(comb_formulas))$ragged
  combs <- sapply(1:nrow(ragged_comb), function(i){
    comb_i <- as.vector(ragged_comb[i, ][ragged_comb[i, ] > 0])
    current_comb <- paste(comb_formulas[comb_i], collapse = " + ")
    current_comb <- paste0("Richness ~ ", current_comb)
  })
  #Add null model in the first line
  combs <- c("Richness ~ 1", combs)
  
  #Fit combination of models
  #Make cluster
  cl <- parallel::makeCluster(10)
  parallel::clusterExport(cl, varlist= c("combs", "d_i"),
                          envir = environment()) #Send objects to nodes
  parallel::clusterEvalQ(cl, {  #Send packages to nodes
    library(MASS) #To fit negative binomial
    library(performance) #To check performance
    library(moranfast) #To calculate moran index
  })
  
  message("Calculating importance...")
  
  comb_m <- pblapply(combs, function(x){
    m_x <- glm.nb(formula = x, data = d_i)
    #RMSE
    rmse_x <- performance::rmse(m_x)
    #loglik
    loglig_x <- logLik(m_x)
    #r2
    r2_x <- performance::r2(m_x)
    
    #Create dataframe
    df_x <- data.frame(variable.combination = x,
                       rmse = rmse_x,
                       logLik = loglig_x,
                       twologlik = m_x$twologlik,
                       r2 = r2_x,
                       deviance = m_x$deviance)
  }, cl = cl)
  parallel::stopCluster(cl)
  #Join results in a unique dataframe
  comb_df <- bind_rows(comb_m)
  
  #Get importance based on twologlik:
  #See: https://www.certara.com/knowledge-base/what-is-the-2ll-or-the-log-likelihood-ratio/
  
  #Partition
  var_p <- partition(comb_df$twologlik, pcan = length(comb_formulas),
                     var.names = comb_formulas)
  
  var_imp <- var_p$I.perc %>% mutate(Variable = row.names(.), .before = 1) %>%
    rename(Importance = ind.exp.var)
  var_imp$Variable <- sub("^(\\S+).*", "\\1", var_imp$Variable)
  var_imp$Variable[grepl("EV", var_imp$Variable)] <- "Spatial"
  var_imp <- var_imp %>%
    mutate(lifeForm = lf_i, .before = 1)
  #Expand grid with NA values for variables that are not selected
  data_fill <- expand.grid(lifeForm = lf_i,
                           Variable = setdiff(names_var,
                                              unique(var_imp$Variable)),
                           Importance = 0)
  var_imp <- rbind(var_imp, data_fill)
  
  #Get data to plot
  pm <- plot_model(m_i, type = "pred")
  data_pm <- lapply(seq_along(pm), function(i){
    d_pm <- pm[[i]]$data %>% as.data.frame()
    d_pm <- d_pm %>% mutate(lifeForm = lf_i, .before = 1)
  }) %>% bind_rows()
  
  # Join model coefficientes
  data_pm <- left_join(data_pm, s, by = join_by("group" == "Variable"))
  
  #Save results
  saveRDS(data_pm,
          paste0("Data/Models_withOccurrences///Predictions/", lf_i,".RDS"))
  saveRDS(var_imp,
          paste0("Data/Models_withOccurrences//Partitioning/", lf_i,".RDS"))
  saveRDS(m_i,
          paste0("Data/Models_withOccurrences//Best_models/", lf_i, ".RDS"))
})

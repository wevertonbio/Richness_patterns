#### Fit best models and partitioning the variance ####
#We didn't scale the variables. See: https://stackoverflow.com/questions/25481331/how-to-use-scale-in-logistic-regression-correctly
library(MASS)
library(dplyr)
library(pbapply)
library(terra)
library(sjPlot)
library(performance)
library(hier.part)

#Com MID ou sem MID?
#Com MID fica mais legal, mas é uma hipótese polêmica :(

#Create directory to save data to plot
dir.create("Data/Models/Partitioning")
dir.create("Data/Models/Predictions")
dir.create("Data/Models/Best_models")

#Import data
my_f <- readRDS("Data/Variables/Formulas.RDS")
all_var <- rast("Data/Variables/Explanatory_Variables.tiff")
my_var <- names(all_var)

#Get richness of lifeforms
lf_indices <- list.files("Data/PAM_indices/", full.names = TRUE)
lf_indices <- pblapply(lf_indices, readRDS)

#Extract values of environmental variables 
all_data <- pblapply(lf_indices, function(i){
  lf_i <- i
  #Create dataframe with lifeform, coordinates and richness
  df_i <- data.frame(lifeform = lf_i$lifeform,
                     lf_i$xy,
                     Richness = lf_i$Richness,
                     Richness_normalized = lf_i$Richness_normalized)
  #Extract values from explanatory variabkes
  pt_var <- terra::extract(all_var, df_i %>% dplyr::select(x, y), ID = FALSE) %>%
    cbind(df_i, .) %>% na.omit()
  return(pt_var)
})
#Rename
names(all_data) <- sapply(all_data, function(x) {unique(x$lifeform)})

#Import candidate models
cm_l <- list.files("Data/Models/Candidate_models/", pattern = ".RDS",
                   full.names = TRUE, recursive = FALSE)

#Generate data to expand grid with NA values for variables that are not selected
names_var <- names(all_var)
names_var <- names_var[!grepl("EV", names_var)]

#To test
i <- cm_l[3]
i

pblapply(cm_l, function(i){
  cm_i <- readRDS(i) #Read i data
  #Get lifeform and endemism
  lf_i <- unique(cm_i$lifeForm)
  
  #Select data
  d_i <- all_data[[lf_i]]
  
  #### Do not include PET ####
  
  #Select best model - do not include PET
  bm <- cm_i %>%
    #Remove PET - Autocorrelation with Bio15
    #filter(!grepl("PET", Formula)) %>% 
    #calculate delta AIC
    mutate(dAIC = AIC - min(AIC, na.rm = T), .after = AIC) %>% #Recalculate AIC
    filter(dAIC == 0) #AIC < 2
  #Get best formulas
  bf <- bm$Formula
  bf
  #Fit model
  m_i <- try(glm.nb(formula = bf,
                    data = d_i), silent = T)
  
  #summary(m_i)
  #plot_model(m_i, type = "pred")
  
  ####Variable importance####
  #Get variables
  vars <- gsub("Richness ~ ", "", bf) %>% strsplit(., " \\+ ") %>% unlist()
  
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
  comb_formulas <- c(vars_no_spatial,  vars_spatial)
  
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
  })
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
    d_pm <- d_pm %>% mutate(lifeForm = lf_i, before = 1)
  }) %>% bind_rows()
  
  
  #Save results
  saveRDS(data_pm,
          paste0("Data/Models/Predictions/", lf_i,".RDS"))
  saveRDS(var_imp,
          paste0("Data/Models/Partitioning/", lf_i,".RDS"))
  saveRDS(m_i,
          paste0("Data/Models/Best_models/", lf_i, ".RDS"))
  
  # #### Do not include Bio15 ####
  # 
  # #Select best model - do not include Bio15
  # bm <- cm_i %>%
  #   #Remove PET - Autocorrelation with Bio15
  #   filter(!grepl("Bio15", Formula)) %>% 
  #   #calculate delta AIC
  #   mutate(dAIC = AIC - min(AIC, na.rm = T), .after = AIC) %>% #Recalculate AIC
  #   filter(dAIC == 0) #AIC < 2
  # #Get best formulas
  # bf <- bm$Formula
  # bf
  # #Fit model
  # m_i <- try(glm.nb(formula = bf,
  #                   data = d_i), silent = T)
  # 
  # #summary(m_i)
  # #plot_model(m_i, type = "pred")
  # 
  # ####Variable importance####
  # #Get variables
  # vars <- gsub("Richness ~ ", "", bf) %>% strsplit(., " \\+ ") %>% unlist()
  # 
  # # Group linear and quadratic variables
  # agrupar_variaveis <- function(var) {
  #   if (grepl("^I\\(", var) & !grepl("EV|Mid_domain|Topoi_het", var)) {
  #     # If there is quadratic, find linear version of the variable
  #     nome_linear <- sub("^I\\(([^\\^]+)\\^2\\)", "\\1", var)
  #     # Check if linear is in the original vector
  #     if (nome_linear %in% vars) {
  #       return(paste(nome_linear, "+", var))
  #     } else {
  #       return(NULL)
  #     }
  #   }
  #   
  #   if(!grepl("^I\\(", var) & !grepl("EV|Mid_domain|Topoi_het", var)){
  #     # Check if there is quadratic version of the linear predictor
  #     var_quadratica <- paste("I(", var, "^2)", sep = "")
  #     if (var_quadratica %in% vars) {
  #       return(paste(var, "+", var_quadratica))
  #     } else {return(var)}
  #   }
  #   
  #   if(grepl("EV|Mid_domain|Topoi_het", var)) {
  #     return(var) }
  # }
  # 
  # # Appply function to group variable
  # vars_agrupadas <- lapply(vars, agrupar_variaveis)
  # vars_agrupadas <- Filter(Negate(is.null), vars_agrupadas) %>% unlist() %>% unique()
  # 
  # # Group spatial variables
  # vars_no_spatial <- vars_agrupadas[!grepl("EV", vars_agrupadas)]
  # vars_spatial <- paste(vars_agrupadas[grepl("EV", vars_agrupadas)],
  #                       collapse = " + ")
  # comb_formulas <- c(vars_no_spatial,  vars_spatial)
  # 
  # #Create combination
  # ragged_comb <- combos(length(comb_formulas))$ragged
  # combs <- sapply(1:nrow(ragged_comb), function(i){
  #   comb_i <- as.vector(ragged_comb[i, ][ragged_comb[i, ] > 0])
  #   current_comb <- paste(comb_formulas[comb_i], collapse = " + ")
  #   current_comb <- paste0("Richness ~ ", current_comb)
  # })
  # #Add null model in the first line
  # combs <- c("Richness ~ 1", combs)
  # #Fit combination of models
  # comb_m <- pblapply(combs, function(x){
  #   m_x <- glm.nb(formula = x, data = d_i)
  #   #RMSE
  #   rmse_x <- performance::rmse(m_x)
  #   #loglik
  #   loglig_x <- logLik(m_x)
  #   #r2
  #   r2_x <- performance::r2(m_x)
  #   
  #   #Create dataframe
  #   df_x <- data.frame(variable.combination = x,
  #                      rmse = rmse_x,
  #                      logLik = loglig_x,
  #                      twologlik = m_x$twologlik,
  #                      r2 = r2_x,
  #                      deviance = m_x$deviance)
  # })
  # #Join results in a unique dataframe
  # comb_df <- bind_rows(comb_m)
  # 
  # #Get importance based on twologlik:
  # #See: https://www.certara.com/knowledge-base/what-is-the-2ll-or-the-log-likelihood-ratio/
  # 
  # #Partition
  # var_p <- partition(comb_df$twologlik, pcan = length(comb_formulas),
  #                    var.names = comb_formulas)
  # 
  # var_imp <- var_p$I.perc %>% mutate(Variable = row.names(.), .before = 1) %>% 
  #   rename(Importance = ind.exp.var)
  # var_imp$Variable <- sub("^(\\S+).*", "\\1", var_imp$Variable)
  # var_imp$Variable[grepl("EV", var_imp$Variable)] <- "Spatial"
  # var_imp <- var_imp %>%
  #   mutate(lifeForm = lf_i, .before = 1)
  # #Expand grid with NA values for variables that are not selected
  # data_fill <- expand.grid(lifeForm = lf_i,
  #                          Variable = setdiff(names_var,
  #                                             unique(var_imp$Variable)),
  #                          Importance = 0)
  # var_imp <- rbind(var_imp, data_fill)
  # 
  # #Get data to plot
  # pm <- plot_model(m_i, type = "pred")
  # data_pm <- lapply(seq_along(pm), function(i){
  #   d_pm <- pm[[i]]$data %>% as.data.frame()
  #   d_pm <- d_pm %>% mutate(lifeForm = lf_i, before = 1)
  # }) %>% bind_rows()
  # 
  # 
  # #Save results
  # saveRDS(data_pm,
  #         paste0("Data/Models/Predictions/", lf_i,"_pet.RDS"))
  # saveRDS(var_imp,
  #         paste0("Data/Models/Partitioning/", lf_i,"_pet.RDS"))
  # saveRDS(m_i,
  #         paste0("Data/Models/Best_models/", lf_i, "_pet.RDS"))
  
})

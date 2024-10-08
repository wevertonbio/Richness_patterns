#### Fit candidate models ####
#Load packages
library(enmpa)
library(terra)
library(dplyr)
library(sjPlot)
library(pbapply)
library(parallel)
library(MASS)
library(performance)
#Calculare I moran
library(moranfast)

####FIT CANDIDATE MODELS####
dir.create("Data/Models/Candidate_models", recursive = TRUE)
#Import data
my_f <- readRDS("Data/Variables/Formulas.RDS") #Formulas
all_var <- rast("Data/Variables/Explanatory_Variables.tiff") #Explanatory Variables
my_var <- names(all_var) #Get names of variables

#Get richness of lifeforms
lf_indices <- list.files("Data/PAM_indices/", full.names = TRUE)
lf_indices <- pblapply(lf_indices, readRDS)

#Plot rasters to check
rich <- rast(pblapply(lf_indices, function(x){
  r_x <- rasterize(x = x$xy, y = all_var, values = x$Richness, fun = "mean",
                   background = 0)
}))
plot(rich[[1]])

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

#Run models
#Test model
x <- all_data[[10]]

#Run models by lifeform
pblapply(all_data, function(x){
  dt <- x
  message("Fitting models for ", unique(dt$lifeform))
  #Normalize data
  #dt[names(all_var)] <- scale(dt[names(all_var)])
  
  #### Run in parallel ####
  cl <- makeCluster(8)
  clusterExport(cl, varlist= c("my_f", "dt"), #Send objects to nodes
                envir=environment())
  clusterEvalQ(cl, {  #Send packages to nodes
    library(MASS) #To fit negative binomial
    library(performance) #To check performance
    library(moranfast) #To calculate moran index
  })
  all_m <- pblapply(seq_along(my_f), function(i){
    tryCatch({
      #Test formula
 
      m_i <- try(glm.nb(formula = my_f[i],
                        data = dt), silent = T)
      #plot_model(m_i, type = "pred")
      
      # #Calculate imoran index on residuals
      moranI <- moranfast(m_i$residuals, dt$x, dt$y)
      
      #Calculate performance
      #Overdispersion - Must be > 0.05
      overdisp_i <- performance::check_overdispersion(m_i)
      disp_ratio_i <- overdisp_i$dispersion_ratio
      disp_p <- overdisp_i$p_value
      
      # #VIF - don't need to calculate with quadratic terms
      # vif_i <- try(performance::check_collinearity(m_i))
      # highest_vif <- try(max(vif_i$VIF, na.rm = TRUE))
      #RMSE
      rmse_i <- performance::rmse(m_i)
      
      #Generate dataframe
      df <- data.frame(lifeForm = dt$lifeform[1],
                       Formula = my_f[i],
                       AIC = m_i$aic,
                       loglik = m_i$twologlik,
                       Moran_I_obs = moranI$observed,
                       Moran_I_p = moranI$p.value,
                       Dispersion_ratio = disp_ratio_i,
                       p_dispersion_ration = disp_p,
                       #Highest_VIF = highest_vif,
                       rmse = rmse_i)
      return(df)
      gc()
    },
    error=function(e) NULL) #Avoid errors
  }, cl = cl)
  
  all_m2 <- Filter(Negate(is.null), all_m)
  df_m <- bind_rows(all_m2)
  # write.csv(df_m,
  #           paste0("Data/Models/Candidate_models/",
  #                  dt$lifeform[1], #Lifeform
  #                  ".csv"),
  #           row.names = F)
  saveRDS(df_m,
          paste0("Data/Models/Candidate_models/",
                 dt$lifeform[1], #Lifeform
                 ".RDS"))
  
  stopCluster(cl)
}) #End of looping

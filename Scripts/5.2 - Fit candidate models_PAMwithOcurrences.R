#### Fit candidate models ####
#Load packages
library(enmpa)
library(terra)
library(dplyr)
#library(sjPlot)
library(pbapply)
library(parallel)
library(MASS)
library(performance)
#Calculare I moran
library(moranfast)

####FIT CANDIDATE MODELS####
dir.create("Data/Models_withOccurrences//Candidate_models", recursive = TRUE)
#Import data
my_f <- readRDS("Data/Variables/Formulas.RDS") #Formulas
all_var <- rast("Data/Variables/Explanatory_Variables.tiff") #Explanatory Variables
my_var <- names(all_var) #Get names of variables

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

#Run models
#Test model
x <- all_data[[2]]

#Run models by lifeform
pblapply(all_data, function(x){
  dt <- x
  # Remove sites with less than 20 species
  dt <- dt %>% filter(Richness >= 20)
  
  message("Fitting models for ", unique(dt$lifeform))
  #Normalize data
  #dt[names(all_var)] <- scale(dt[names(all_var)])
  
  #### Run in parallel ####
  cl <- makeCluster(18)
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
  
  saveRDS(df_m,
          paste0("Data/Models_withOccurrences/Candidate_models/",
                 unique(dt$lifeform), #Lifeform
                 ".RDS"))
  
  stopCluster(cl)
}) #End of looping

#### Prepare variables to model ####
#Load packages
library(enmpa)
library(terra)
library(dplyr)
library(sjPlot)
library(pbapply)
library(parallel)
library(MASS)
library(performance)
library(ape)

####PREPARE VARIABLES####
dir.create("Richness_Models")
#Import AF extent
af <- vect("Vetores/AF_dissolved.shp")

####Get variables####
var <- list.files("Current_Neotropic/", pattern = ".tif", full.names = T) %>%
  rast()
var <- var[[c(#"Bio01", "Bio12",
  "Bio07", "Bio15", "Bio02", #Stability
  "Bio06", "Bio14")]]
var <- crop(var, af, mask = TRUE)
#Aggregate
f_ag <- round(0.08333333/res(var)[1],0)
var <- terra::aggregate(var, fact = f_ag, fun = mean)
res(var)
#Others variables
ol_var <- list.files("Others_variables/", pattern = ".tif", full.names = T)
o_var <- pblapply(seq_along(ol_var), function(i){
  o_i <- rast(ol_var[i])
  print(names(o_i))
  o_i <- crop(o_i, af, mask = TRUE)
  o_i <- terra::resample(o_i, var[[1]])
  return(o_i)
})
o_var <- rast(o_var)
names(o_var)
#All variables
all_var <- c(var, o_var)
#Subset and reorder variables
names(all_var) %>% dput()
all_var <- all_var[[c("Aridity", #Energy 
                      "Bio06", "Bio14", #Tolerance
                      "Bio02", "Bio07", "Bio15", #Seasonality
                      "Mid_domain", #Mid Domain
                      "Prec_stab", "Temp_stab", #Stability
                      "Topoi_het",#Topographic heterogeneity
                      "EV1", "EV2")]] #Spatial
#Save variables
writeRaster(all_var, "Data/Variables/Explanatory_Variables.tiff", overwrite = TRUE)

####See correlation between variables####
df_cor <- as.data.frame(all_var) %>% na.omit() %>% scale() %>% cor()
# correlation_finder(cor_mat = df_cor ,threshold = 0.7,verbose = T)
#Save variables correlation
write.csv(df_cor, "Data/Variables/Correlation_variables.csv")

#Create combinations of variables without correlation
#Identify correlated variables
cor_groups <- function(cor_mat, th = 0.7){
  var_names <- colnames(cor_mat)
  cor_mat<- abs(cor_mat)
  diag(cor_mat) <- 0
  cor_var <- cor_mat > th
  cor_var <- apply(cor_var, 2, function(x) colnames(cor_mat)[x])
  if (length(cor_var) == 0) {
    cor_var <- "No pair of variables reached the specified correlation threshold."
    var_list <- var_names
  } else{
    cor_var <- Filter(function(x) length(x) > 0, cor_var)
    #Get variables with correlation
    cor_names <- names(cor_var)
    #Get all combinations of order of variables
    l_var <- gtools::permutations(n = length(cor_names), r = length(cor_names),
                                  v = cor_names,
                                  repeats.allowed = FALSE)
    
    # Save in list
    l_var <- lapply(1:nrow(l_var), function(i){
      q_i <- as.character(unlist(l_var[i,]))
      return(q_i)
    })
    
    #Iterate on minha_lista
    var_list <- pblapply(seq_along(l_var), function(x){
      l_x <- l_var[[x]]
      the_var <- var_names
      # Inicialize o índice da iteração
      i <- 1
      while(i <= length(l_x)) {
        f_i <- l_x[[i]]
        if(f_i %in% the_var) {
          cor_var_i <- cor_var[f_i][[1]]
          the_var <- the_var[!(the_var %in% cor_var_i)] } else {
            the_var <- the_var
          }
        # Avance para a próxima iteração
        i <- i + 1
      }
      return(the_var)
    })
    return(unique(var_list))
  }
}

#Get var comb
var_comb <- cor_groups(df_cor, th = 0.7)
#Save variables combination
saveRDS(var_comb, "Data/Variables/var_comb.RDS")

#Get formulas: at least 4 variables
my_f <- pblapply(seq_along(var_comb), function(i){
  ind_i <- var_comb[[i]]
  f_i <- enmpa::get_formulas(dependent = "Richness",
                             independent = ind_i,
                             type = c("lq"), minvar = 4)
  return(f_i)
})
my_f2 <- unique(unlist(my_f))

#Remove models without any spatial variable
my_f2_spt <- my_f2[grepl("EV", my_f2)]

# Drop off formulas with quadratic predictor without its linear function
filter_quadratic_without_linear <- function(ff) {
  #Get variables
  f1 <- gsub("Richness ~ ", "", ff)
  f1 <- strsplit(f1, " \\+ ") %>% unlist()
  #Remove spatial predictors from the list
  f1 <- f1[!grepl("EV", f1)]
  
  # Has quadratic?
  q <- grepl("\\^2)", ff)
  #If has quadratic, search for linear version
  if(q) {
    quadratic_predictors <- f1[grepl("\\^2", f1)] %>% gsub("I\\(|\\^2)", "", .)
    linear_predictors <- f1[!grepl("\\^2", f1)]
    return(all(quadratic_predictors %in% linear_predictors))
  } else {
    return(TRUE)
  }
}
# Filter
my_f3 <- my_f2_spt[pbsapply(my_f2_spt, filter_quadratic_without_linear)]

#Filter formulas: remove quadratic models with Mid_domain and Topoi_het
my_f4 <- my_f3[which(!grepl("Mid_domain\\^2|Topo_het\\^2|Topoi_het\\^2", my_f3))]

#Save formulas
saveRDS(my_f4, "Data/Variables/Formulas.RDS")



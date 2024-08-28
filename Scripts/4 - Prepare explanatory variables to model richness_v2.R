#### Prepare variables to model ####
#Load packages
library(enmpa)
library(terra)
library(dplyr)
library(pbapply)
library(RcppAlgos)



####PREPARE VARIABLES####
dir.create("Richness_Models")
#Load atlantic forest limits
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_limite_integrador.gpkg")
#Buffer of 10km
af <- buffer(af, width = 10*1000)
plot(af)


####Get variables####
var <- list.files("Current_Neotropic", pattern = ".tif", full.names = T) %>%
  rast()
var <- var[[c(#"Bio01", "Bio12",
  "Bio02", "Bio07", "Bio15",  #Seasonality
  "Bio06", "Bio11", "Bio17", "Bio14")]] #Tolerance
var <- crop(var, af, mask = TRUE)
#Aggregate
f_ag <- round(0.08333333/res(var)[1],0)
var <- terra::aggregate(var, fact = f_ag, fun = mean)
res(var)
#Others variables
ol_var <- list.files("Others_variables", pattern = ".tif", full.names = T)

o_var <- pblapply(seq_along(ol_var), function(i){
  o_i <- rast(ol_var[i])
  print(names(o_i))
  o_i <- crop(o_i, af, mask = TRUE)
  o_i <- terra::resample(o_i, var[[1]])
  return(o_i)
})
o_var <- rast(o_var)
names(o_var)
#Rename
names(o_var) <- c("Mid_domain", "PET", "Aridity", "Prec_stab", "Lat_EV", 
                  "Long_EV", "Temp_stab", "Topo_het", "Topoi_het")

#All variables
all_var <- c(var, o_var)
plot(all_var)
#Get correlation
df_cor <- as.data.frame(all_var) %>% na.omit() %>% cor()

#Subset and reorder variables
names(all_var) %>% dput()
all_var <- all_var[[c("Aridity", "PET", #Energy 
                      "Bio06", "Bio14", #Tolerance
                      "Bio02", "Bio07", "Bio15", #Seasonality
                      "Mid_domain", #Mid Domain
                      "Prec_stab", "Temp_stab", #Stability
                      "Topo_het",#Topographic heterogeneity
                      "Lat_EV", "Long_EV")]] #Spatial
#Save variables
writeRaster(all_var, "Data/Variables/Explanatory_Variables.tiff", overwrite = TRUE)
all_var <- terra::rast("Data/Variables/Explanatory_Variables.tiff")

####See correlation between variables####
df_cor <- as.data.frame(all_var) %>% na.omit() %>% cor()
# z <- ntbox::correlation_finder(cor_mat = df_cor ,threshold = 0.7,verbose = T)
#Save variables correlation
write.csv(df_cor, "Data/Variables/Correlation_variables.csv")

#Create combinations of variables without correlation
#Get variables (remove spatial predictors)
v <- names(all_var)
v

all_comb <- pbsapply(v, function(x){
  c_x <- df_cor %>% as.data.frame() %>% dplyr::select(all_of(x)) %>% pull %>% abs()
  to_rm <- which(c_x > 0.7)
  to_keep <- c(x, row.names(df_cor[-to_rm,]))
  to_keep <- to_keep[to_keep %in% v]
  #Append spatial eingevectors
  v_x <- sort(to_keep)
  return(v_x)
})

#Get formulas: at least 4 variables
my_f <- pblapply(seq_along(all_comb), function(i){
  ind_i <- all_comb[[i]]
  f_i <- enmpa::get_formulas(dependent = "Richness",
                             independent = ind_i,
                             type = "lq", minvar = 4,
                             mode = "moderate")
  #Convert to formulas
  f_i <- sapply(f_i, as.formula)
  
  #Make sure there is at least 4 variables
  f_i <- f_i[sapply(f_i, function(x) length(attr(terms(x), "term.labels")) >= 4)]

  return(as.character(f_i))
})
my_f2 <- unique(unlist(my_f))
my_f2 %>% as.data.frame() %>% View()

#Remove models without any spatial variable
my_f2_spt <- my_f2[grepl("_EV", my_f2)]
my_f2_spt %>% as.data.frame() %>% View()

# Drop off formulas with quadratic predictor without its linear function
filter_quadratic_without_linear <- function(ff) {
  #Get variables
  f1 <- gsub("Richness ~ ", "", ff)
  f1 <- strsplit(f1, " \\+ ") %>% unlist()
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
my_f3 %>% as.data.frame() %>% View()

#Filter formulas: remove quadratic models with Mid_domain, Topo_het and Lat
my_f4 <- my_f3[which(!grepl("Mid_domain\\^2|Topo_het\\^2|Lat_EV\\^2", my_f3))]
#See formulas
my_f4 %>% as.data.frame() %>% View()

#Save formulas
saveRDS(my_f4, "Data/Variables/Formulas.RDS")



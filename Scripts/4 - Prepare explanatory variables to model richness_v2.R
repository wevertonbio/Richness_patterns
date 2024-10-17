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
#Get correlation only linear
df_cor <- as.data.frame(all_var) %>% na.omit() %>% cor()
#Get correlation (linear and quadratic)
d2 <- as.data.frame(all_var) %>% na.omit()
dq <- apply(d2, 2, function(x) x^2)
colnames(dq) <- paste0("I(", colnames(dq), ")^2")
dq <- cbind(d2, dq)
df_cor2 <- cor(dq)

#Save variables correlation
write.csv(df_cor, "Data/Variables/Correlation_variables.csv")
#Remove some variables with high correlation

#Subset and reorder variables
names(all_var) %>% dput()
all_var <- all_var[[c("Aridity",  #Energy 
                      "PET", 
                      "Bio06",  #Tolerance
                      #"Bio14", #Remove Bio14 because it has high correlation with Bio15
                      "Bio02", "Bio15", #Seasonality
                      #"Bio07", #Remove Bio07 because it has high correlation with Bio02 and Bio06
                      "Mid_domain", #Mid Domain
                      "Prec_stab", "Temp_stab", #Stability
                      "Topo_het",#Topographic heterogeneity
                      "Lat_EV", "Long_EV")]] #Spatial
#Save variables
writeRaster(all_var, "Data/Variables/Explanatory_Variables.tiff", overwrite = TRUE)
all_var <- terra::rast("Data/Variables/Explanatory_Variables.tiff")

#### Get formulas ####
z <- enmpa::get_formulas(dependent = "Richness",
                         independent = names(all_var),
                         type = "lq",
                         minvar = 4, #Minimum of 4 variables
                         mode = "intensive")
length(z)
head(z)

#Remove models without any spatial variable
my_f<- z[grepl("_EV", z)]
my_f %>% as.data.frame() %>% View()

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
my_f2 <- my_f[pbsapply(my_f, filter_quadratic_without_linear)]
my_f2 %>% as.data.frame() %>% View()

#Filter formulas: remove quadratic models with Mid_domain^2, Topo_het^2 and Lat^2
my_f3 <- my_f2[which(!grepl("Mid_domain\\^2|Topo_het\\^2|Lat_EV\\^2", my_f2))]
#See formulas
my_f3 %>% as.data.frame() %>% View()

#Remove formulas with correlated variables
library(tidyr)
head(df_cor2)
#Get pairwise correlations
df_cor2 <- df_cor2 %>% as.data.frame() %>% mutate(var1 = row.names(.), .before = 1)
par_cor <- gather(data = df_cor2, key = "var2", value = "correlation", -var1)
#Get only correlation > 0.7
cor7 <- par_cor %>% filter(correlation > 0.7 | correlation < -0.7)
#Remove formulas with PET and Bio15
pet15 <- pbsapply(my_f3, function(x){
  is_pet <- grepl("PET", x)
  is_15 <- grepl("Bio15", x)
  all(is_pet, is_15)
  })
my_f4 <- my_f3[!pet15]
as.data.frame(my_f4) %>% View()


#Save formulas
saveRDS(my_f4, "Data/Variables/Formulas.RDS")



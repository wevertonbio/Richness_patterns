#### Prepare variables to model ####
#Load packages
library(enmpa)
library(terra)
library(dplyr)
library(pbapply)



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
saveRDS(df_cor2, "Data/Variables/Correlation_variables.rds")
#Remove some variables with high correlation

#Subset and reorder variables
names(all_var) %>% dput()
all_var <- all_var[[c("Aridity",  #Energy 
                      "PET", 
                      "Bio06",  #Tolerance
                      #"Bio14", #Remove Bio14 because it has high correlation with Bio15
                      "Bio02", "Bio15", #Seasonality
                      #"Bio07", #Remove Bio07 because it has high correlation with Bio02 and Bio06
                      #"Mid_domain", #Mid Domain
                      "Prec_stab", "Temp_stab", #Stability
                      "Topo_het",#Topographic heterogeneity
                      "Lat_EV", "Long_EV")]] #Spatial
#Save variables
writeRaster(all_var, "Data/Variables/Explanatory_Variables.tiff", overwrite = TRUE)
# Import again
all_var <- terra::rast("Data/Variables/Explanatory_Variables.tiff")

#### Get formulas ####
# Do not consider spatial for now
z <- enmpa::get_formulas(dependent = "Richness",
                         independent = c("Aridity", "PET", "Bio06", "Bio02",
                                         "Bio15", "Prec_stab", "Temp_stab", 
                                         "Topo_het"),
                         type = "lqp",
                         minvar = 6, #Minimum of 6 variables
                         mode = "moderate")
length(z)
head(z)
data.frame(z) %>% View()

# # Remove models with less than 6 predictors
np <- sapply(z, function(i){
  all.vars(i %>% as.formula)[-1] %>% length()
})
z <- z[np >= 6]
data.frame(z) %>% View()

# Add spatial variables
z_spt <- enmpa::get_formulas(dependent = "Richness",
                             independent = c("Lat_EV", "Long_EV"),
                             type = "lq",
                             minvar = 1, 
                             mode = "intensive")
my_f <- pblapply(z, function(i){
  sapply(z_spt, function(x){
    x_i<- gsub("Richness ~ ", "", x)
    paste0(i, " + ", x_i)
  })
}) %>% unlist()
my_f %>% data.frame() %>% View()

# # #Remove models without any spatial variable
# my_f<- z[grepl("_EV", z2)]
# my_f %>% as.data.frame() %>% View()

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

#Filter formulas: remove quadratic models with Mid_domain^2, Topo_het^2, Prec_stab^2 e Temp_stab^2
my_f3 <- gsub(" + I(Topo_het^2)", "", my_f2, fixed = TRUE)
my_f3 <- gsub(" + I(Prec_stab^2)", "", my_f3, fixed = TRUE)
my_f3 <- gsub(" + I(Temp_stab ^2)", "", my_f3, fixed = TRUE)
my_f3 <- unique(my_f3)

#See formulas
my_f3 %>% as.data.frame() %>% View()

#Remove formulas with correlated variables
library(tidyr)
# df_cor2 <- read.csv("Data/Variables/Correlation_variables.csv")
head(df_cor2)
#Get pairwise correlations
df_cor2 <- df_cor2 %>% as.data.frame() %>% 
  mutate(var1 = row.names(.), .before = 1)
par_cor <- gather(data = df_cor2, key = "var2", value = "correlation", -var1)
#Get only correlation > 0.7
cor7 <- par_cor %>% filter(correlation > 0.7 | correlation < -0.7)

# Ignore correlation with spatial
cor7 <- cor7 %>% filter(!grepl("Lat|Long", var1)) %>% 
  filter(!grepl("Lat|Long", var2))

# Remove reverse duplicates
cor7 <- cor7 %>%
   mutate(
    par_ordenado_min = pmin(var1, var2),
    par_ordenado_max = pmax(var1, var2)
  ) %>%
  filter(var1 != var2) %>%
  distinct(par_ordenado_min, par_ordenado_max, .keep_all = TRUE) %>%
  dplyr::select(var1, var2, correlation)

# Remove correlation between linear and its quadratic term
cor7 <- cor7 %>% mutate(var_A = gsub("I\\(|\\)|\\^2", "", var1),
                        var_B = gsub("I\\(|\\)|\\^2", "", var2)) %>% 
  filter(var_A != var_B)

#Remove formulas with predictors with correlation > |0.7|
f_with_correlation <- pblapply(1:length(my_f3), function(i){
  f_i <- my_f3[[i]] %>% as.formula() %>% all.vars()
  # Check if correlated variables are part of the formula
  var_in_formula <- sapply(1:nrow(cor7), function(x){
    v_x1 <- cor7$var1[x] %in% f_i
    v_x2 <- cor7$var2[x] %in% f_i
    all(v_x1, v_x2)
  })
  any(var_in_formula)
})
f_with_correlation <- unlist(f_with_correlation)
#Subset and keep only formulas without correlation
my_f4 <- my_f3[!f_with_correlation]
as.data.frame(my_f4) %>% View()

#Save formulas
saveRDS(my_f4, "Data/Variables/Formulas_v2.RDS")


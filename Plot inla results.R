#### Plot results ####
library(dplyr)
library(INLA)
library(terra)
library(tidyterra)
library(ggplot2)
library(data.table)
library(pbapply)

#### Main analyzes with modeled species ####

# Import model results
res <- readRDS("Data/INLA_models/models.rds")

# Plot fixed effects
fixed <- lapply(names(res), function(x) {
  f_x <- res[[x]]$model$summary.fixed
  f_x <- f_x %>% mutate(lifeform = x,
                        Variable = row.names(.), 
                        .before = 1)
  }) %>% rbindlist()

# Remove intercept
fixed <- fixed %>% filter(Variable != "intercept") %>% 
  filter(lifeform != "All")
str(fixed)

# Identify positive and negative effects
fixed <- fixed %>% mutate(Effect = ifelse(mean < 0, "Negative", "Positive"))
# Identify non-significant variables
fixed$Effect[fixed$`0.025quant` < 0 & fixed$`0.975quant` > 0] <- "Non-significant"

# Set factors
# Variables
v <- factor(x = fixed$Variable,
            levels = c("Topo_het", "Temp_stab", "Prec_stab", "Bio02", "Bio06", "PET", 
                       "Aridity"),
            labels = c("Topographic\nheterogeneity", "Temperature\nstability", 
                       "Precipitation\nstability", "Temperature\nseasonality", 
                       "Cold extremes", "Energy availability", 
                       "Water stress"))
# Lifeforms
fixed$lifeform <- factor(x = fixed$lifeform,
             levels = c("Tree", "Liana", "Shrub", "Subshrub", 
                        "Terrestrial_herb", "Epiphytic_herb"),
             labels = c("Tree", "Liana", "Shrub", "Subshrub", 
                        "Terrestrial herb", "Epiphytic herb"))

g_fixed <- ggplot(fixed, aes(x = v, y = mean, colour = Effect)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "#0072B2") +
  scale_color_manual(values = c("#D55E00", "gray20", "#009E73")) +
  coord_flip() +
  theme_minimal() +
  labs(y = "Estimates", x = "Variables") +
  facet_wrap(lifeform ~ ., scales = "fixed") + 
  ggpubr::theme_pubclean() +
  theme(legend.position = "bottom",
        panel.border=element_rect(colour="black",size=0.5, fill = NA))
g_fixed
ggsave("Figures/Model_results.png",
       g_fixed,
       dpi = 600, units = "px",
       width = 5350,
       height = 4500, limitsize = F)

#### Variable importance ####
imp <- lapply(names(res), function(x) {
  res[[x]]$var_imp
}) %>% rbindlist()

# Remove results for all
imp <- imp %>% filter(lifeForm != "All")

#Get colors
unique(imp$Variable) %>% as.character() %>% sort()
myc <- c("#D55E00", #Aridity
         "#800020", #PET
         "blue",  #Bio06
         "#E69F00", #Bio02
         #"#F0E442", #Bio07
         "#3C5488FF", #Bio15
         "#F39B7FFF", #Temp_stab
         "#56B4E9", #Prec_stab
         "#7E6148FF") #Topo_het

#Rename colors
names(myc) <- c("Water stress", "Energy availability", 
                "Cold extremes", "Temperature\nseasonality",
                "Precipitation\nstability", 
                "Temperature\nstability", 
                "Topographic\nheterogeneity")
#Make factors
imp$Variable <-factor(imp$Variable,
                      levels = c("Aridity",
                                 "PET",
                                 "Bio06",
                                 "Bio02",
                                 "Temp_stab", "Prec_stab",
                                 "Topo_het"),
                      labels = c("Water stress", "Energy availability", 
                                 "Cold extremes", "Temperature\nseasonality",
                                 "Precipitation\nstability", 
                                 "Temperature\nstability", 
                                 "Topographic\nheterogeneity"))
unique(imp$lifeForm)
imp$lifeForm <- factor(imp$lifeForm,
                       levels = c("Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Terrestrial_herb", "Epiphytic_herb"
                       ),
                       labels = c("Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Terrestrial\nherb", "Epiphytic\nherb"
                       ))

#Plot
g_imp <- ggplot(data = imp, aes(x = Variable, y = Importance, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = myc) +
  facet_wrap(.~lifeForm, nrow = 2) +
  ggpubr::theme_pubclean() +
  #theme_stata() +
  #theme_wsj()+ 
  theme(axis.text.x= element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.border=element_rect(colour="black",size=1, fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 14),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 16),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))
g_imp
#Save
ggsave("Figures/VariableImportance.png",
       g_imp, dpi = 600, units = "px", width = 2550,
       height = 1450, scale = 3)

#### Complementary nalyzes including rare and undersampled species ####

# Import model results
res_wo <- readRDS("Data/INLA_models/models_with_Ocurrences.rds")

# Plot fixed effects
fixed_wo <- lapply(names(res_wo), function(x) {
  f_x <- res_wo[[x]]$model$summary.fixed
  f_x <- f_x %>% mutate(lifeform = x,
                        Variable = row.names(.), 
                        .before = 1)
}) %>% rbindlist()

# Remove intercept
fixed_wo <- fixed_wo %>% filter(Variable != "intercept") %>% 
  filter(lifeform != "All")
str(fixed_wo)

# Identify positive and negative effects
fixed_wo <- fixed_wo %>% mutate(Effect = ifelse(mean < 0, "Negative", "Positive"))
# Identify non-significant variables
fixed_wo$Effect[fixed_wo$`0.025quant` < 0 & fixed_wo$`0.975quant` > 0] <- "Non-significant"

# Set factors
# Variables
v_wo <- factor(x = fixed_wo$Variable,
            levels = c("Topo_het", "Temp_stab", "Prec_stab", "Bio02", "Bio06", "PET", 
                       "Aridity"),
            labels = c("Topographic\nheterogeneity", "Temperature\nstability", 
                         "Precipitation\nstability", "Temperature\nseasonality", 
                         "Cold extremes", "Energy availability", 
                         "Water stress"))
# Lifeforms
fixed_wo$lifeform <- factor(x = fixed_wo$lifeform,
                         levels = c("Tree", "Liana", "Shrub", "Subshrub", 
                                    "Terrestrial_herb", "Epiphytic_herb"),
                         labels = c("Tree", "Liana", "Shrub", "Subshrub", 
                                    "Terrestrial herb", "Epiphytic herb"))

g_fixed_wo <- ggplot(fixed_wo, aes(x = v, y = mean, colour = Effect)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "#0072B2") +
  scale_color_manual(values = c("#D55E00", "gray20", "#009E73")) +
  coord_flip() +
  theme_minimal() +
  labs(y = "Estimates", x = "Variables") +
  facet_wrap(lifeform ~ ., scales = "fixed") + 
  ggpubr::theme_pubclean() +
  theme(legend.position = "bottom",
        panel.border=element_rect(colour="black",size=0.5, fill = NA))
g_fixed_wo
ggsave("Figures/Model_results_with_Undersampled.png",
       g_fixed_wo,
       dpi = 600, units = "px",
       width = 5350,
       height = 4500, limitsize = F)

#### Variable importance ####
imp_wo <- lapply(names(res_wo), function(x) {
  res_wo[[x]]$var_imp
}) %>% rbindlist()

# Remove res_woults for all
imp_wo <- imp_wo %>% filter(lifeForm != "All")

#Get colors
unique(imp_wo$Variable) %>% as.character() %>% sort()
myc <- c("#D55E00", #Aridity
         "#800020", #PET
         "blue",  #Bio06
         "#E69F00", #Bio02
         #"#F0E442", #Bio07
         "#3C5488FF", #Bio15
         "#F39B7FFF", #Temp_stab
         "#56B4E9", #Prec_stab
         "#7E6148FF") #Topo_het

#Rename colors
names(myc) <- c("Water stress", "Energy availability", 
                "Cold extremes", "Temperature\nseasonality",
                "Precipitation\nstability", 
                "Temperature\nstability", 
                "Topographic\nheterogeneity")
#Make factors
imp_wo$Variable <-factor(imp_wo$Variable,
                      levels = c("Aridity",
                                 "PET",
                                 "Bio06",
                                 "Bio02",
                                 "Temp_stab", "Prec_stab",
                                 "Topo_het"),
                      labels = c("Water stress", "Energy availability", 
                                 "Cold extremes", "Temperature\nseasonality",
                                 "Precipitation\nstability", 
                                 "Temperature\nstability", 
                                 "Topographic\nheterogeneity"))
unique(imp_wo$lifeForm)
imp_wo$lifeForm <- factor(imp_wo$lifeForm,
                       levels = c("Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Terrestrial_herb", "Epiphytic_herb"
                       ),
                       labels = c("Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Terrestrial\nherb", "Epiphytic\nherb"
                       ))

#Plot
g_imp_wo <- ggplot(data = imp_wo, aes(x = Variable, y = Importance, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = myc) +
  facet_wrap(.~lifeForm, nrow = 2) +
  ggpubr::theme_pubclean() +
  #theme_stata() +
  #theme_wsj()+ 
  theme(axis.text.x= element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.border=element_rect(colour="black",size=1, fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 14),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 16),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill = guide_legend(byrow = TRUE))
g_imp_wo
#Save
ggsave("Figures/VariableImportance_with_Undersampled.png",
       g_imp_wo, dpi = 600, units = "px", width = 2550,
       height = 1450, scale = 3)

####Plot maps####
#Load packages
library(terra)
library(geobr)
library(tidyterra)
library(ggspatial)
library(scales)
library(patchwork)

#Import vectors
#South America
sa <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/South_America.gpkg")
br <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/Brazil_States.gpkg")

####Richness####
#Get richness of lifeforms
lf_indices <- readRDS("Data/Richness_by_lifeform.rds")
# Get xy
xy <- lf_indices$xy
lf_indices$xy <- NULL

#Remove palms and bamboos
names(lf_indices)
lf_indices <- lf_indices[!grepl("Palm_tree|Bamboo", names(lf_indices))]

#Get lifeforms
lf_names <- names(lf_indices)

#Reorder lifeforms
lf_names <- c("Tree", "Liana", "Shrub", "Subshrub", "Terrestrial_herb",
              "Epiphytic_herb")

#Reorder lf_indices
lf_indices <- lf_indices[lf_names]

#Rasterize indices
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_limite_integrador.gpkg")
r_base <- rast(ext = ext(af), res = 0.08333333)
#Test
#x <- lf_indices[[1]]
r <- pblapply(lf_indices, function(x){
  rx <- rasterize(as.matrix(xy), r_base,
                  values = x)
  rx[rx == 0] <- NA
  # Remove sites with less than 20 species
  rx[rx < 20] <- NA
  names(rx) <- "Richness"
  return(rx)
}) %>% rast()
plot(r[[1]])
plot(r)

#Get box limits to plot
bb_af <- ext(r[[1]])

#Replace _ by space in names r
names(r) <- gsub("_", " ", names(r))

#Test
#i <- 2
p <- pblapply(1:length(names(r)), function(i) {
  r_i <- r[[i]]
  #Convert raster to dataframe
  df_r <- as.data.frame(r_i, xy = TRUE) %>% na.omit()
  
  #Lifeform
  lf_i <- colnames(df_r)[3]
  
  #Change columns name
  colnames(df_r)[3] <- "Richness"
  
  #set braks
  b <- round(seq(min(df_r$Richness), max(df_r$Richness), length.out = 5), 0)
  
  #Get label
  my_label <- paste0("(", letters[i], ") ", lf_i)
  
  g <- ggplot() +
    geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
    geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
    geom_raster(data = df_r, aes(x, y, fill = Richness), alpha = 1) +
    scale_fill_gradientn(colors = rev(pals::brewer.rdylbu(10)),
                         name = "Richness", breaks = b) +
    geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
    
    # scale_fill_gradientn(
    #   colors = rev(c("#9DBF9E", "#FCB97D", "#A84268")),
    #   na.value = "grey80",
    #   #limits = c(0, 1),
    #   #oob = scales::squish,
    #   name = "Richness") + 
    geom_hline(yintercept = - 19, linetype = "dashed") +
    coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[2] + 0.5),
             ylim = c(bb_af[3] - 0.5, ymax=bb_af[4] + 0.5),
             expand = T) +
    theme(legend.position = c(0.75, 0.21),
          legend.direction = "vertical",
          text = element_text(family = "Arial"),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.75, 'cm'), #change legend key size
          legend.key.height = unit(0.87, 'cm'), #change legend key height
          legend.key.width = unit(1.3, 'cm'),
          legend.box = "horizontal",
          legend.background = element_rect(fill = "white", size = 0.5, colour = "black"),
          panel.background = element_rect(fill = 'aliceblue', colour = NA),
          panel.border = element_rect(colour = "black", size = 2, fill = NA),
          #plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    annotation_scale(pad_x = unit(2.5, "cm")) +
    xlab("Longitude") + ylab("Latitude") +
    metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6) +
    annotate("text", x = -57.8, y=-4, label = my_label, size = unit(5, "pt"),
             hjust = 0)
  g
})
#Get lifeforms
# lf_names <- pbsapply(lf_indices, function(x){ x$lifeform})
names(p) <- names(r)
# #Get letters to legend 
# my_legend <- paste0("(", letters[1:6], ") ", lf_names)

#Arrange plot
pp <- wrap_plots(p) + 
  plot_layout(ncol = 3, axes = "collect")

#Save
ggsave("Figures/Richness.png",
       pp, dpi = 600, units = "px", width = 1450,
       height = 980, scale = 7)

####Richness - With Occurrences ####
#Get richness of lifeforms
lf_indices_wo <- readRDS("Data/Richness_by_lifeform_withOccurrences.rds")
# Get xy
xy_wo <- lf_indices_wo$xy
lf_indices_wo$xy <- NULL

#Remove All, palms and bamboos
names(lf_indices_wo)
lf_indices_wo <- lf_indices_wo[!grepl("Palm_tree|Bamboo", names(lf_indices_wo))]

#Get lifeforms
lf_names_wo <- names(lf_indices_wo)

#Reorder lifeforms
lf_names_wo <- c("Tree", "Liana", "Shrub", "Subshrub", "Terrestrial_herb",
                 "Epiphytic_herb")

#Reorder lf_indices
lf_indices_wo <- lf_indices_wo[lf_names_wo]

#Rasterize indices
#Test
#x <- lf_indices[[1]]
r_wo <- pblapply(lf_indices_wo, function(x){
  rx <- rasterize(as.matrix(xy_wo), r_base,
                  values = x)
  rx[rx == 0] <- NA
  # Remove sites with less than 20 species
  rx[rx < 20] <- NA
  names(rx) <- "Richness"
  return(rx)
}) %>% rast()
plot(r_wo[[1]])
plot(r_wo)

#Get box limits to plot
bb_af <- ext(r_wo[[1]])

#Replace _ by space in names r
names(r_wo) <- gsub("_", " ", names(r_wo))

#Test
#i <- 2
p_wo <- pblapply(1:length(names(r_wo)), function(i) {
  r_i <- r_wo[[i]]
  #Convert raster to dataframe
  df_r <- as.data.frame(r_i, xy = TRUE) %>% na.omit()
  
  #Lifeform
  lf_i <- colnames(df_r)[3]
  
  #Change columns name
  colnames(df_r)[3] <- "Richness"
  
  #set braks
  b <- round(seq(min(df_r$Richness), max(df_r$Richness), length.out = 5), 0)
  
  #Get label
  my_label <- paste0("(", letters[i], ") ", lf_i)
  
  g <- ggplot() +
    geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
    geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
    geom_raster(data = df_r, aes(x, y, fill = Richness), alpha = 1) +
    scale_fill_gradientn(colors = rev(pals::brewer.rdylbu(10)),
                         name = "Richness", breaks = b) +
    geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
    
    # scale_fill_gradientn(
    #   colors = rev(c("#9DBF9E", "#FCB97D", "#A84268")),
    #   na.value = "grey80",
    #   #limits = c(0, 1),
    #   #oob = scales::squish,
    #   name = "Richness") + 
    geom_hline(yintercept = - 19, linetype = "dashed") +
    coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[2] + 0.5),
             ylim = c(bb_af[3] - 0.5, ymax=bb_af[4] + 0.5),
             expand = T) +
    theme(legend.position = c(0.75, 0.21),
          legend.direction = "vertical",
          text = element_text(family = "Arial"),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.75, 'cm'), #change legend key size
          legend.key.height = unit(0.87, 'cm'), #change legend key height
          legend.key.width = unit(1.3, 'cm'),
          legend.box = "horizontal",
          legend.background = element_rect(fill = "white", size = 0.5, colour = "black"),
          panel.background = element_rect(fill = 'aliceblue', colour = NA),
          panel.border = element_rect(colour = "black", size = 2, fill = NA),
          #plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    annotation_scale(pad_x = unit(2.5, "cm")) +
    xlab("Longitude") + ylab("Latitude") +
    metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6) +
    annotate("text", x = -57.8, y=-4, label = my_label, size = unit(5, "pt"),
             hjust = 0)
  g
})
#Get lifeforms
# lf_names <- pbsapply(lf_indices, function(x){ x$lifeform})
names(p_wo) <- names(r_wo)
# #Get letters to legend 
# my_legend <- paste0("(", letters[1:6], ") ", lf_names)

#Arrange plot
pp_wo <- wrap_plots(p_wo) + 
  plot_layout(ncol = 3, axes = "collect")

#Save
ggsave("Figures/Richness_with_Undersampled.png",
       pp_wo, dpi = 600, units = "px", width = 1450,
       height = 980, scale = 7)


####Predictors####
r_pred <- rast("Data/Variables/Explanatory_Variables.tiff")
plot(r_pred)


#Change names
names(r_pred) <- c("Aridity", "PET", "Bio06", "Bio02", "Bio15",
                   "Precipitation stability", "Temperature stability",
                   "Topographic heterogeneity",
                   "Spatia filter (Lat)", "Spatial filter (Long)")
plot(r_pred)

# Remove spatial rasters (not used)
r_pred <- r_pred[[c("Aridity", "PET", "Bio06", "Bio02", "Bio15",
                  "Precipitation stability", "Temperature stability",
                  "Topographic heterogeneity")]]

#Test
#i <- 1
p_pred <- pblapply(1:length(names(r_pred)), function(i) {
  r_i <- r_pred[[i]]
  #Convert raster to dataframe
  df_r <- as.data.frame(r_i, xy = TRUE) %>% na.omit()
  
  #Variable
  var_i <- colnames(df_r)[3]
  
  #Change columns name
  colnames(df_r)[3] <- "Value"
  
  # #set braks
  # b <- round(seq(min(df_r$Value, na.omit = T), max(df_r$Value, na.omit = T),
  #                length.out = 3), 1)
  # b[3] <- b[3] - 0.01
  
  #Get label
  my_label <- paste0("(", letters[i], ") ", var_i)
  
  g <- ggplot() +
    geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
    geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
    geom_raster(data = df_r, aes(x, y, fill = Value), alpha = 1) +
    scale_fill_gradientn(colors = rev(pals::brewer.spectral(10)),
                         name = var_i)+
    geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
    
    # scale_fill_gradientn(
    #   colors = rev(c("#9DBF9E", "#FCB97D", "#A84268")),
    #   na.value = "grey80",
    #   #limits = c(0, 1),
    #   #oob = scales::squish,
    #   name = "Richness") + 
    geom_hline(yintercept = - 19, linetype = "dashed") +
    coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[2] + 0.5),
             ylim = c(bb_af[3] - 0.5, ymax=bb_af[4] + 0.5),
             expand = T) +
    theme(legend.position = c(0.79, 0.22),
          legend.direction = "vertical",
          text = element_text(family = "Arial"),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.75, 'cm'), #change legend key size
          legend.key.height = unit(0.87, 'cm'), #change legend key height
          legend.key.width = unit(1.3, 'cm'),
          legend.box = "horizontal",
          legend.background = element_rect(fill = "white", size = 0.5, colour = "black"),
          panel.background = element_rect(fill = 'aliceblue', colour = NA),
          panel.border = element_rect(colour = "black", size = 2, fill = NA),
          #plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    annotation_scale(pad_x = unit(2.5, "cm")) +
    xlab("Longitude") + ylab("Latitude") +
    metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6) +
    annotate("text", x = -57.8, y=-4, label = my_label, size = unit(5, "pt"),
             hjust = 0)
  g
})
#Get variables names
names(p_pred) <- names(r_pred)

pp_pred2 <- wrap_plots(p_pred) + plot_layout(ncol = 4)
ggsave("Figures/Predictors.png",
       pp_pred2, dpi = 600, units = "px", width = 1350,
       height = 800, scale = 10.25)



# #Get letters to legend 
# my_legend <- paste0("(", letters[1:6], ") ", lf_names)


####Plot of correlation####
library(ggcorrplot)
library(rstatix)
pred_cor <- r_pred %>% as.data.frame() %>% na.omit() %>% cor() %>% round(digits = 2)
#Vector of correlations
v_cor <- pred_cor %>% pull_lower_triangle(diagonal = FALSE) %>% as.matrix() %>% #Get the lower part
  as.numeric() %>% abs() %>% na.omit()

#Convert to color
v_col <- ifelse(v_cor >= 0.7, "black", "gray45")
outline_color <- ifelse(v_cor >= 0.7, "black", "gray")
lab_size <- ifelse(v_cor >= 0.7, 3.5, 2.75)

g_cor <- ggcorrplot(pred_cor, method = "square", type = "lower",
                    lab = TRUE, lab_col = v_col, outline.color = outline_color,
                    lab_size = lab_size, ggtheme = ggpubr::theme_pubclean(),
                    legend.title = "Correlation") +
  theme(legend.position = "right")

g_cor
#Save
ggsave("Figures/Correlations.png",
       g_cor, dpi = 600, units = "px", width = 2500,
       height = 2500, scale = 2.25)


####Plot map with ecoregions ####
library(pals)
library(scales)

#Import ecoregions
eco <- vect("Data/Ecoregions_af.gpkg")
plot(eco)
#Import vectors
#South America
sa <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/South_America.gpkg")
br <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/Brazil_States.gpkg")
#Get names
names.eco <- eco$ECO_NAME
#Get colors
show_col(pals::okabe())

eco.col <- c("#D55E00", #Alto Paraná Atlantic forests
             "#0072B2", #Araucaria moist forests
             "green", #Bahia coastal forests"
             "#782AB6", #Bahia interior forests"
             "#B8860B",#"Cerrado" 
             "firebrick", #Brazilian Atlantic dry forests
             "#F0E442", #Caatinga Enclaves moist forests
             "#CC79A7", #Campos
             "#90AD1C", #Maranhão Babaçu forests
             "#FBE426",  #"Restingas/Mangrooves"
             "#C20088", #Pernambuco coastal forests"
             "#993F00", #Pernambuco interior forests"
             "#009E73", #Serra do Mar coastal forests
             "#56B4E9") #Uruguayan savanna
eco.col2 <- alpha(eco.col,
                  alpha = c(1, #Alto Paraná Atlantic forests
                            1, #Araucaria moist forests
                            1, #Bahia coastal forests"
                            1, #Bahia interior forests"
                            1, #"Cerrado" 
                            0, #Brazilian Atlantic dry forests
                            0, #Caatinga Enclaves moist forests
                            1, #Campos Rupestres montane savanna
                            0, #Maranhão Babaçu forests
                            0,  #Northeast Brazil restingas
                            0, #Pernambuco coastal forests"
                            0, #Pernambuco interior forests"
                            1, #Serra do Mar coastal forests
                            1)) #Uruguayan savanna
show_col(eco.col)

#Get box limits to plot
bb_sa <- ext(sa)
bb_af <- ext(eco)

#Main map South America
main_map <- ggplot() +
  geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
  geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
  geom_sf(data = eco, aes(fill = ECO_NAME), colour = NA) +
  scale_fill_manual(values = eco.col) +
  geom_sf(data = br, colour = "black", alpha = 0) +
  geom_hline(yintercept = - 19, linetype = "dashed") +
  coord_sf(xlim = c(-82, -34),
           ylim = c(-55, bb_sa[4]),
           expand = T) +
  #theme_void() +
  annotate("rect", xmin=bb_af[1] - 0.5, xmax=bb_af[3] + 0.5,
           ymin=bb_af[2] - 0.5, ymax=bb_af[4] + 0.5,
           alpha=0, fill="black", color = "black", linewidth = 2) +
  theme(legend.position = "none",
        legend.direction = "vertical",
        text = element_text(family = "Arial"),
        legend.title = element_text(family = "Roboto Black", size = 11),
        legend.text = element_text(family = "Roboto Medium", size = 10),
        legend.key.size = unit(0.75, 'cm'), #change legend key size
        legend.key.height = unit(0.87, 'cm'), #change legend key height
        legend.key.width = unit(0.87, 'cm'),
        legend.box = "horizontal",
        legend.background = element_rect(fill = NA, size = 0.5, colour = "black"),
        panel.background = element_rect(fill = 'aliceblue', colour = NA),
        panel.border = element_rect(colour = "black", size = 2, fill = NA)) +
  annotation_scale(location = "br") + 
  xlab("Longitude") + ylab("Latitude") +
  metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6)
main_map

#Zoomed map
g_eco <- ggplot() +
  geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
  geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
  geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
  geom_sf(data = eco, aes(fill = ECO_NAME), colour = NA) +
  scale_fill_manual(values = eco.col, name = "Ecoregions") +
  geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
  geom_hline(yintercept = - 19, linetype = "dashed") +
  coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[3] + 0.5),
           ylim = c(bb_af[2] - 0.5, ymax=bb_af[4] + 0.5),
           expand = T) +
  theme(text = element_text(family = "Arial"),
        legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.75, 'cm'), #change legend key size
        legend.key.height = unit(0.87, 'cm'), #change legend key height
        legend.key.width = unit(1.3, 'cm'),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", size = 0.5, colour = "black"),
        panel.background = element_rect(fill = 'aliceblue', colour = NA),
        panel.border = element_rect(colour = "black", size = 2, fill = NA),
        #plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  annotation_scale(pad_x = unit(2.5, "cm")) +
  xlab("Longitude") + ylab("Latitude") +
  metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6)
g_eco
# #Save map to edit in CANVA
# ggsave("Figures/Ecoregions_to_edit.png",
#        g_eco, dpi = 600, units = "px", width = 2500,
#        height = 2500, scale = 1.5)

#Join map with patchwork
f1 <- main_map + g_eco + 
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag.position = "topleft")
f1
ggsave("Figures/Ecoregions.png",
       f1, dpi = 600, width = 10, height = 4.4, scale = 1.25)

# #### Simplify ecoregions ####
# library(mapview)
# 
# #Import shapefiles
# af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_limite_integrador.gpkg")
# eco <- vect("c:/Users/wever/Downloads/Ecoregions2017/Ecoregions2017.shp")
# plot(af)
# #Crop
# eco_af <- crop(eco, af)
# plot(eco_af)
# eco_af$ECO_NAME %>% unique()
# mapview(eco_af, zcol = "ECO_NAME", burst = TRUE)
# 
# #Change names
# # eco_af$ECO_NAME[eco_af$ECO_NAME == "Atlantic Coast restingas"] <- "Atlantic Coast Restingas"
# # eco_af$ECO_NAME[eco_af$ECO_NAME == "Northeast Brazil restingas"] <- "Northeast Restingas"
# eco_af$ECO_NAME[eco_af$ECO_NAME == "Caatinga Enclaves moist forests"] <- "Caatinga"
# 
# #Remove
# to_remove <- c("Chiquitano dry forests", "Pantanal", "Atlantic Coast restingas",
#                "Northeast Brazil restingas", "Mangrooves",
#                "Fernando de Noronha-Atol das Rocas moist forests", 
#                "Humid Chaco",
#                "Southern Cone Mesopotamian savanna", 
#                "Southern Atlantic Brazilian mangroves")
# eco_af <- subset(eco_af, !(eco_af$ECO_NAME %in% to_remove))
# 
# #Plot
# eco_af
# 
# #Convert to raster
# r_base <- rast(ext = ext(af), res = 0.008333333, values = 0)
# eco_r <- rasterize(eco_af, r_base, field = "ECO_NAME")
# plot(eco_r)
# #Get levels
# l <- levels(eco_r)
# 
# #Fill NAs
# eco_r2 <- focal(eco_r, 29, "modal", na.policy="only")
# levels(eco_r2) <- l
# eco_r2 <- crop(eco_r2, af, mask = TRUE)
# plot(eco_r2)
# #Convert to vector again
# eco_v <- as.polygons(eco_r2)
# mapview(eco_v) + mapview(eco_af, zcol = "ECO_NAME") + mapview(af)
# 
# #Save
# writeVector(eco_v, "Data/Ecoregions_af.gpkg",  overwrite = T)
# 

#### Plot residuals ####
# Get life forms
lf_names <- c("Tree", "Liana", "Shrub", "Subshrub", "Terrestrial_herb",
              "Epiphytic_herb")

#### With modeled species ####
# Residuals
residuos <- pblapply(lf_names, function(x){
  # For test
  # x <- res[[1]]
  
  # Get data
  d_i <- res[[x]]$data
  
  # Remove pixels with less than 1 species
  d_i <- d_i %>% filter(Richness >= 1)
  
  # Check residuals
  index <- inla.stack.index(res[[x]]$stk, tag="est")$data
  d_i$predicted <- res[[x]]$model$summary.fitted.values[index, "mean"]
  d_i$residuals <- d_i$predicted - d_i$Richness
  r_residuals <- rasterize(as.matrix(d_i[,c("x", "y")]),
                           r_base, values = d_i$residuals)
  # plot(r_residuals, main = x)
  # Check moran
  #moran <- moranfast::moranfast(d_i$residuals, d_i$x, d_i$y)
  # moran
  return(r_residuals)
}) %>% rast()
names(residuos) <- lf_names
plot(residuos)

# Plot
p_residuals <- pblapply(lf_names, function(i) {
  # Get residuals
  r_i <- residuos[[i]]
  
  g <- ggplot() +
    geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
    geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
    geom_spatraster(data = r_i) +
    scale_fill_whitebox_c(name = "Residuals") +
    geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
    geom_hline(yintercept = - 19, linetype = "dashed") +
    coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[2] + 0.5),
             ylim = c(bb_af[3] - 0.5, ymax=bb_af[4] + 0.5),
             expand = T) +
    theme(legend.position = c(0.86, 0.21),
          legend.direction = "vertical",
          text = element_text(family = "Arial"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.72, 'cm'), #change legend key size
          legend.key.height = unit(0.82, 'cm'), #change legend key height
          legend.key.width = unit(1.1, 'cm'),
          legend.box = "horizontal",
          legend.background = element_rect(fill = "white", size = 0.5, colour = "black"),
          panel.background = element_rect(fill = 'aliceblue', colour = NA),
          panel.border = element_rect(colour = "black", size = 2, fill = NA),
          #plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    annotation_scale(pad_x = unit(2.5, "cm")) +
    xlab("Longitude") + ylab("Latitude") +
    metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6) +
    facet_wrap(.~ paste(i))
  g
})

names(p_residuals) <- lf_names

#Arrange plot
pp_residuals <- wrap_plots(p_residuals) + 
  plot_layout(ncol = 3, axes = "collect")

#Save
ggsave("Figures/Residuals.png",
       pp_residuals, dpi = 600, units = "px", width = 1000,
       height = 930, scale = 8)


#### With undersampled species ####
# Residuals
residuos_wo <- pblapply(lf_names, function(x){
  # For test
  # x <- res[[1]]
  
  # Get data
  d_i <- res_wo[[x]]$data
  
  # Remove pixels with less than 1 species
  d_i <- d_i %>% filter(Richness >= 1)
  
  # Check residuals
  index <- inla.stack.index(res_wo[[x]]$stk, tag="est")$data
  d_i$predicted <- res_wo[[x]]$model$summary.fitted.values[index, "mean"]
  d_i$residuals <- d_i$predicted - d_i$Richness
  r_residuals <- rasterize(as.matrix(d_i[,c("x", "y")]),
                           r_base, values = d_i$residuals)
  # plot(r_residuals, main = x)
  # Check moran
  #moran <- moranfast::moranfast(d_i$residuals, d_i$x, d_i$y)
  # moran
  return(r_residuals)
}) %>% rast()
names(residuos_wo) <- lf_names
plot(residuos_wo)

# Plot
p_residuals_wo <- pblapply(lf_names, function(i) {
  # Get residuals
  r_i <- residuos_wo[[i]]
  
  g <- ggplot() +
    geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
    geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
    geom_spatraster(data = r_i) +
    scale_fill_whitebox_c(name = "Residuals") +
    geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
    geom_hline(yintercept = - 19, linetype = "dashed") +
    coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[2] + 0.5),
             ylim = c(bb_af[3] - 0.5, ymax=bb_af[4] + 0.5),
             expand = T) +
    theme(legend.position = c(0.86, 0.21),
          legend.direction = "vertical",
          text = element_text(family = "Arial"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.72, 'cm'), #change legend key size
          legend.key.height = unit(0.82, 'cm'), #change legend key height
          legend.key.width = unit(1.1, 'cm'),
          legend.box = "horizontal",
          legend.background = element_rect(fill = "white", size = 0.5, colour = "black"),
          panel.background = element_rect(fill = 'aliceblue', colour = NA),
          panel.border = element_rect(colour = "black", size = 2, fill = NA),
          #plot.margin = unit(c(0,0,0,0), "cm")
    ) +
    annotation_scale(pad_x = unit(2.5, "cm")) +
    xlab("Longitude") + ylab("Latitude") +
    metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6) +
    facet_wrap(.~ paste(i))
  g
})

names(p_residuals_wo) <- lf_names

#Arrange plot
pp_residuals_wo <- wrap_plots(p_residuals_wo) + 
  plot_layout(ncol = 3, axes = "collect")

#Save
ggsave("Figures/Residuals_with_Undersampled.png",
       pp_residuals_wo, dpi = 600, units = "px", width = 1000,
       height = 930, scale = 8)


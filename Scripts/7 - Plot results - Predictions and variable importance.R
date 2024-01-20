#### Plot results of Richness models####

#Load packages
library(dplyr)
library(ggplot2)
library(pbapply)
library(ggthemes)
library(ggsci)

#Save plots
dir.create("Figures")

#### Effects/Predictions ####
lp <- list.files("Data/Models/Predictions/", full.names = TRUE)
lp <- pblapply(lp, readRDS)
dfp <- bind_rows(lp)

#Remove other lifeforms and 
d <- dfp %>%
  filter(lifeForm != "Other")
unique(d$lifeForm)

#Set colors
  #See colors
palette.colors(palette = "Okabe-Ito")
scales::show_col(ggsci::pal_npg("nrc")(10))
scales::show_col(palette.colors(palette = "Okabe-Ito"))
unique(d$group) %>% as.character()

myc <- c("#D55E00", #Aridity
         #"#0072B2", #Bio11
          "blue",  #Bio06
         "#DC0000FF", #Bio14
         "#E69F00", #Bio02
         "#F0E442", #Bio07
         "#3C5488FF", #Bio15
         "#F39B7FFF", #Temp_stab
         "#56B4E9", #Prec_stab
         "#7E6148FF", #Topoi_het
         "#000002", #Mid-domain
         #"#000002", #Spatial EV1
         "#999999") #Spatial (Ev2)

#Make factors
d$group <- factor(d$group,
                 levels = c("Aridity",
                            "Bio06", "Bio14",
                            "Bio02","Bio07", "Bio15",
                            "Temp_stab", "Prec_stab",
                            "Topoi_het", "Mid_domain", "EV2"),
                 labels = c("Aridity",
                            "Minimum\nTemperature of\nColdest Month",
                            "Precipitation\nof Driest\nMonth",
                            "Mean Diurnal\nRange",
                            "Temperature\nAnnual Range",
                            "Precipitation\nSeasonality",
                            "Temperature\nstability",
                            "Precipitation\nstability",
                            "Topographic\nheterogeneity",
                            "Mid-Domain", "Spatial\nfilter 2"))
d$lifeForm <- factor(d$lifeForm,
                       levels = c("All", "Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Herb"))
         
#Plot
gd <- ggplot(d, aes(x = x, y = predicted, colour = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
                  colour = NULL), alpha = 0.5) +
  geom_line() +
  scale_colour_manual(values = myc) +
  scale_fill_manual(values = myc) +
  facet_grid(lifeForm ~ group, scales = "free") +
  xlab("Variables") + ylab("Predicted richness") + 
  #theme_bw() +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", 
                                        colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        axis.text.x = element_text(angle = 45, hjust=1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x =  element_text(size = 20),
        axis.title.y =  element_text(size = 20),
        strip.text = element_text(size = 11))
X11()
gd
#Save as PNG file
ggsave("Figures/Predictions.png", gd, units = "px",
       dpi = 600, width = 3000,
       height = 1700, scale = 2.75)



#### Variables importance ####
li <- list.files("Data/Models/Partitioning/", full.names = TRUE)
li <- pblapply(li, readRDS)
dfi <- bind_rows(li)

#Remover others
dfi <- dfi %>% filter(lifeForm != "Other")

#Keep only variables with some importance
sum_imp_by_var <- dfi %>% group_by(Variable) %>% 
  summarise(Total = sum(Importance))
var_with_imp <- sum_imp_by_var %>% filter(Total > 0) %>% pull(Variable)
dfi <- dfi %>% filter(Variable %in% var_with_imp)


#Get colors
unique(dfi$Variable) %>% as.character() %>% sort()

myc <- c("#D55E00", #Aridity
         #"#0072B2", #Bio11
         "blue",  #Bio06
         "#DC0000FF", #Bio14
         "#E69F00", #Bio02
         "#F0E442", #Bio07
         "#3C5488FF", #Bio15
         "#F39B7FFF", #Temp_stab
         "#56B4E9", #Prec_stab
         "#7E6148FF", #Topoi_het
         "#000002", #Mid-domain
         #"#000002", #Spatial EV1
         "#999999") #Spatial (Ev3)



#Make factors
dfi$Variable <-factor(dfi$Variable,
                      c("Aridity",
                        "Bio06", "Bio14",
                        "Bio02","Bio07", "Bio15",
                        "Temp_stab", "Prec_stab",
                        "Topoi_het", "Mid_domain", "Spatial"),
                      labels = c("Aridity",
                                 "Minimum Temperature\nof Coldest Month",
                                 "Precipitation of\nDriest Month",
                                 "Mean Diurnal\nRange",
                                 "Temperature Annual\nRange",
                                 "Precipitation\nSeasonality",
                                 "Temperature\nstability",
                                 "Precipitation\nstability",
                                 "Topographic\nheterogeneity",
                                 "Mid-Domain", "Spatial\nfilter 2"))
unique(dfi$lifeForm)
dfi$lifeForm <- factor(dfi$lifeForm,
                       levels = c("All", "Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Herb"))

g_imp <- ggplot(data = dfi, aes(x = Variable, y = Importance, fill = Variable)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = myc) +
  facet_wrap(.~lifeForm) +
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
       g_imp, dpi = 600, units = "px", width = 2500,
       height = 1700, scale = 2.75)

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
lf_indices <- list.files("Data/PAM_indices/", full.names = TRUE)
#Remove others
lf_indices <- lf_indices[!grepl("Other", lf_indices)]

#Read data
lf_indices <- pblapply(lf_indices, readRDS)

#Get lifeforms
lf_names <- sapply(lf_indices, function(x) x$lifeform)
names(lf_indices) <- lf_names
#Reorder lifeforms
lf_names <- c("All", "Tree", "Liana", "Shrub", "Subshrub", "Herb")

#Reorder lf_indices
lf_indices <- lapply(lf_names, function(x) lf_indices[[x]])
names(lf_indices) <- lf_names

#Rasterize indices
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_dissolved.gpkg")
r_base <- rast(ext = ext(af), res = 0.08333333)
#Test
#x <- lf_indices[[1]]
r <- pblapply(lf_indices, function(x){
  rx <- rasterize(x$xy, r_base,
            values = x$Richness)
  rx[rx == 0] <- NA
  names(rx) <- "Richness"
  return(rx)
}) %>% rast()
plot(r[[1]])

#Get box limits to plot
bb_af <- ext(r[[1]])



#Test
#i <- 6
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
    coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[3] + 0.5),
             ylim = c(bb_af[2] - 0.5, ymax=bb_af[4] + 0.5),
             expand = T) +
    theme(legend.position = c(0.75, 0.19),
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
lf_names <- pbsapply(lf_indices, function(x){ x$lifeform})
names(p) <- lf_names
# #Get letters to legend 
# my_legend <- paste0("(", letters[1:6], ") ", lf_names)

#Arrange plot
pp <- (p$All + p$Tree + p$Liana) / (p$Shrub + p$Subshrub + p$Herb)

#Save
ggsave("Figures/Richness.png",
       pp, dpi = 600, units = "px", width = 2500,
       height = 1700, scale = 5)

####Predictors####
r_pred <- rast("Data/Variables/Explanatory_Variables.tiff")
plot(r_pred)
#Remove Bio11
r_pred <- r_pred[[!grepl("Bio11", names(r_pred))]]
#Change names
names(r_pred) <- c("Aridity", "PET",
                   "Bio06", "Bio14", "Bio17", 
                    "Bio02", "Bio04", "Bio07", "Bio15",
                   "Mid-domain",
                   "Precipitation stability", "Temperature stability",
                   "Topographic heterogeneity",
                   "Spatia filter (1)", "Spatial filter (2)")
plot(r_pred)

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
    coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[3] + 0.5),
             ylim = c(bb_af[2] - 0.5, ymax=bb_af[4] + 0.5),
             expand = T) +
    theme(legend.position = c(0.75, 0.19),
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

pp_pred2 <- cowplot::plot_grid(plotlist = p_pred)
ggsave("Figures/Predictors.png",
       pp_pred2, dpi = 600, units = "px", width = 2500,
       height = 2500, scale = 6)

# #Get letters to legend 
# my_legend <- paste0("(", letters[1:6], ") ", lf_names)

# #Arrange plot
# pp_pred <- (p_pred[["Aridity"]] + p_pred[["Bio06"]] + p_pred[["Bio14"]] + p_pred[["Mid-domain"]]) / 
#   (p_pred[["Bio07"]] + p_pred[["Bio15"]] + p_pred[["Precipitation\nstability"]] +
#      p_pred[["Temperature\nstability"]]) /
#   (p_pred[["Topographic\nheterogeneity"]] + p_pred[["Spatial\nfilter (1)"]] +  p_pred[["Spatial\nfilter (2)"]])
#   
# #Save
# ggsave("Figures/Predictors.png",
#        pp_pred, dpi = 600, units = "px", width = 2500,
#        height = 1900, scale = 5)

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


####Table with best models####
#Import candidate models
cm_l <- list.files("Data/Models/Candidate_models/", pattern = ".RDS",
                   full.names = TRUE, recursive = FALSE)

#To test
i <- cm_l[7]

bm <- pblapply(cm_l, function(i){
  cm_i <- readRDS(i) #Read i data
  #Get lifeform and endemism
  lf_i <- unique(cm_i$lifeForm)
  
  #Select best model
  bm <- cm_i %>%
    #filter(Highest_VIF <= 10) %>% #VIF < 10
    filter(p_dispersion_ration > 0.05) %>% #Dispersion ratio < 0.05
    #calculate delta AIC
    mutate(dAIC = AIC - min(AIC, na.rm = T), .after = AIC) %>% #Recalculate AIC
    filter(dAIC == 0) #AIC < 2
  #Select columns
  bm <- bm %>% dplyr::select(c(lifeForm, Formula, AIC, dAIC,
                               twologlik = loglik, Moran_I, rmse))
  return(bm)
  }) %>% bind_rows()
#Write table
# library(flextable)
# flextable(bm)
write.csv(bm, "Data/Best_models_table.csv", row.names = FALSE)

#Edit table
bm_df <- bm %>% dplyr::select(-Formula) %>% 
  flextable() %>% 
  align(align = "center", part = "all") %>% 
  theme_apa()
bm_df
#Save as docx
save_as_docx(bm_df, path = "Data/Best_models_table2.docx")

#Reshape table to inclue variables as colnames
#Get all variables
vars <- lapply(bm$Formula, function(x) all.vars(as.formula(x))) %>% unlist() %>% 
  unique()
#Exclude Richness
vars <- vars[-1]
#Get all lifeforms
lfs <- unique(bm$lifeForm)
#Exclude others and reoders
lfs %>% dput()
lfs <- c("All", "Tree", "Liana", "Shrub", "Subshrub",  "Herb")

lf_l <- pblapply(lfs, function(i){
  lf_i <- i
  bm_i <- bm %>% filter(lifeForm == lf_i)
  f_i <-  bm_i$Formula
  var_l <- lapply(vars, function(x){
    var_i <- x
    var_i2 <- paste0("I(", var_i, "^2)")
    
    lf_var <- ifelse(grepl(var_i, f_i), "L", "X")
    lf_var2 <- ifelse(grepl(var_i2, f_i, fixed = T), "Q", lf_var)
    
    df_var <- data.frame(var_i = lf_var2)
    colnames(df_var) <- var_i
    return(df_var)
  })
  d <- bind_cols(var_l) %>% mutate(lifeform = lf_i, .before = 1) %>% 
    mutate(AIC = round(bm_i$AIC,0),
           twologlik = round(bm_i$twologlik,1),
           Moran_I = round(bm_i$Moran_I,2),
           rmse = round(bm_i$rmse, 1))
  dl <- bind_rows(d)
  return(dl)
})
df <- bind_rows(lf_l)
#Reorder and rename variables
colnames(df) %>% dput()
df <- df %>% dplyr::select(lifeform, Aridity, 
                             Bio06, Bio14, Bio07, Bio15,
                             Temp_stab,
                             Prec_stab,
                             Mid_domain,
                             Topoi_het,
                             Ev3,
                             AIC, twologlik, 
                             Moran_I, rmse)
#Save
write.csv(df, "Data/Best_models_table_2.csv", row.names = FALSE)
#Beuutiful table with flextable
bdf <- flextable(df) %>% theme_vanilla() %>% 
  set_header_labels(Temp_stab = "Temperature Stability",
                    Prec_stab = "Precipitation Stability",
                    Mid_domain = "Mid-domain",
                    Topoi_het = "Topographic heterogeneity",
                    Ev3 = "Spatial (Ev3)",
                    Moran_I = "Moran Index") %>% 
  add_header_row(values = c("", "Energy", "Tolerance", "Seasonality",
                            "Stability", "", "", "", "", "", "", ""),
                 colwidths = c(1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1)) %>% 
  align(align = "center", part = "all") 
bdf
#Save as docx
save_as_docx(bdf, path = "Data/Best_models_table.docx")
save_as_image(bdf, path = "Data/Best_models_table.png")
save_as_pptx(bdf, path = "Data/Best_models_table.pptx")

####Plot map with hottest ecoregions
library(pals)
library(scales)

#Import ecoregions
eco <- st_read("../spatial_files/Data/Ecoregions_Atlantic_Forest_simplified.gpkg")[1]
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
             "firebrick", #Brazilian Atlantic dry forests
             "#F0E442", #Caatinga Enclaves moist forests
             "#CC79A7", #Campos Rupestres montane savanna
             "#90AD1C", #Maranhão Babaçu forests
             "#FBE426",  #Northeast Brazil restingas
             "#C20088", #Pernambuco coastal forests"
             "#993F00", #Pernambuco interior forests"
             "#009E73", #Serra do Mar coastal forests
             "#56B4E9") #Uruguayan savanna
eco.col2 <- alpha(eco.col,
                 alpha = c(1, #Alto Paraná Atlantic forests
                           1, #Araucaria moist forests
                           1, #Bahia coastal forests"
                           1, #Bahia interior forests"
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
bb_af <- ext(eco)


#Plot
g_eco <- ggplot() +
  geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
  geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
  geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
  geom_sf(data = eco, aes(fill = ECO_NAME), col = eco.col) +
  scale_fill_manual(values = eco.col) +
  geom_sf(data = br, fill = NA, size = 0.1, colour = "grey40") +
  geom_hline(yintercept = - 19, linetype = "dashed") +
  coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[3] + 0.5),
           ylim = c(bb_af[2] - 0.5, ymax=bb_af[4] + 0.5),
           expand = T) +
  theme(text = element_text(family = "Arial"),
        legend.position = "none",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.75, 'cm'), #change legend key size
        legend.key.height = unit(0.87, 'cm'), #change legend key height
        legend.key.width = unit(1.3, 'cm'),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", size = 0.5, colour = "black"),
        panel.background = element_rect(fill = 'aliceblue', colour = NA),
        #panel.border = element_rect(colour = "black", size = 2, fill = NA),
        #plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  annotation_scale(pad_x = unit(2.5, "cm")) +
  xlab("Longitude") + ylab("Latitude") +
  metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6)
g_eco
#Save map to edit in CANVA
ggsave("Figures/Ecoregions_to_edit.png",
       g_eco, dpi = 600, units = "px", width = 2500,
       height = 2500, scale = 1.5)

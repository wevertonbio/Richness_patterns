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
#PET
lp <- list.files("Data/Models/Predictions/",
                 full.names = TRUE)
lp_pet <- lp[grepl("pet", lp)]
lp_pet <- pblapply(lp_pet, readRDS)
d_pet <- bind_rows(lp_pet) %>% mutate("set" = "pet")

#Bio15
lp_bio15 <- lp[grepl("bio15", lp)]
lp_bio15 <- pblapply(lp_bio15, readRDS)
d_bio15 <- bind_rows(lp_bio15) %>% mutate("set" = "bio15")

#Join data
d <- rbind(d_pet, d_bio15)

#Remove bamboos and palm trees
d <- d %>% filter(!(lifeForm %in% c("Bamboo", "Palm_tree")))

#Set colors
  #See colors
palette.colors(palette = "Okabe-Ito")
scales::show_col(ggsci::pal_npg("nrc")(10))
scales::show_col(palette.colors(palette = "Okabe-Ito"))
unique(d$group) %>% as.character()

myc <- c("#D55E00", #Aridity
         "#800020", #PET
          "blue",  #Bio06
         #"#DC0000FF", #Bio14
         "#E69F00", #Bio02
        #"#F0E442", #Bio07
         "#3C5488FF", #Bio15
         "#F39B7FFF", #Temp_stab
         "#56B4E9", #Prec_stab
         "#7E6148FF", #Topoi_het
         "#000002") #Mid-domain
         #"grey50", #Spatial EV1
         #"grey45") #Spatial (Ev2)

#Make factors
unique(d$group)
d$group <- factor(d$group,
                 levels = c("Aridity",
                            "PET",
                            "Bio06",
                            #"Bio14",
                            "Bio02",
                            #"Bio07",
                            "Bio15",
                            "Temp_stab", "Prec_stab",
                            "Topo_het", "Mid_domain"
                            #"Long_EV", "Lat_EV"
                            ),
                 labels = c("Aridity",
                            "PET",
                            "Minimum\nTemp. of\nColdest Month",
                            #"Precipitation\nof Driest\nMonth",
                            "Mean Diurnal\nRange",
                            #"Temperature\nAnnual Range",
                            "Precipitation\nSeasonality",
                            "Temperature\nstability",
                            "Precipitation\nstability",
                            "Topographic\nheterogeneity",
                            "Mid-Domain"
                            #"Spatial\nfilter X",
                            #"Spatial\nfilter Y"
                            ))
#Rename colors
names(myc) <- c("Aridity",
                "PET",
                "Minimum\nTemp. of\nColdest Month",
                "Mean Diurnal\nRange",
                "Precipitation\nSeasonality",
                "Temperature\nstability",
                "Precipitation\nstability",
                "Topographic\nheterogeneity",
                "Mid-Domain")


unique(d$lifeForm)
d$lifeForm <- factor(d$lifeForm,
                     levels = c("All", "Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Herb"
                                  #"Bamboo", "Palm_tree"
                                ),
                     labels = c("All", "Tree", "Liana",
                                "Shrub", "Subshrub",
                                "Herb" 
                                #"Bamboo", "Palm tree"
                                ))
         
#Plot
pblapply(unique(d$set), function(x){
  d_x <- d[d$set == x,] %>% filter(!is.na(group))
  #Subset colors
  myc_x <- myc[names(myc) %in% unique(d_x$group)]
  d_x$group <- droplevels(d_x$group)
  
  gd <- ggplot(d_x, aes(x = x, y = predicted, colour = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
                    colour = NULL), alpha = 0.5) +
    geom_line() +
    scale_colour_manual(values = myc_x, na.translate = FALSE) +
    scale_fill_manual(values = myc_x, na.translate = FALSE) +
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
  # X11()
  # gd
  #Save as PNG file
  filename <- paste0("Figures/Predictions_", x, ".png")
  ggsave(filename, gd, units = "px",
         dpi = 600, width = 1750,
         height = 1400, scale= 3.5)
})

#### Variables importance ####
li <- list.files("Data/Models/Partitioning/", full.names = TRUE)
li_pet <- li[grepl("pet", li)]
li_pet <- pblapply(li_pet, readRDS)
dfi_pet <- bind_rows(li_pet) %>% mutate("set" = "pet")
li_bio15 <- li[grepl("bio15", li)]
li_bio15 <- pblapply(li_bio15, readRDS)
dfi_bio15 <- bind_rows(li_bio15) %>% mutate("set" = "bio15")

#Merge data
dfi <- rbind(dfi_pet, dfi_bio15)

#Remove bamboos and palm trees
dfi <- dfi %>% filter(!(lifeForm %in% c("Bamboo", "Palm_tree")))
unique(dfi$lifeForm)
# 
# 
# #Keep only variables with some importance
# sum_imp_by_var <- dfi %>% group_by(Variable) %>% 
#   summarise(Total = sum(Importance))
# var_with_imp <- sum_imp_by_var %>% filter(Total > 0) %>% pull(Variable)
# dfi <- dfi %>% filter(Variable %in% var_with_imp)

#Get colors
unique(dfi$Variable) %>% as.character() %>% sort()
myc <- c("#D55E00", #Aridity
         "#800020", #PET
         "blue",  #Bio06
         #"#DC0000FF", #Bio14
         "#E69F00", #Bio02
         #"#F0E442", #Bio07
         "#3C5488FF", #Bio15
         "#F39B7FFF", #Temp_stab
         "#56B4E9", #Prec_stab
         "#7E6148FF", #Topoi_het
         "#000002", #Mid-domain
         "#999999") #Spatial (Ev3)
#Rename colors
names(myc) <- c("Aridity",
                "PET",
                "Minimum\nTemperature of\nColdest Month",
                #"Precipitation\nof Driest\nMonth",
                "Mean Diurnal\nRange",
                #"Temperature\nAnnual Range",
                "Precipitation\nSeasonality",
                "Temperature\nstability",
                "Precipitation\nstability",
                "Topographic\nheterogeneity",
                "Mid-Domain", 
                "Spatial")


#Make factors
dfi$Variable <-factor(dfi$Variable,
                      levels = c("Aridity",
                                 "PET",
                                 "Bio06",
                                 #"Bio14",
                                 "Bio02",
                                 #"Bio07",
                                 "Bio15",
                                 "Temp_stab", "Prec_stab",
                                 "Topo_het", "Mid_domain",
                                 "Spatial"),
                      labels = c("Aridity",
                                 "PET",
                                 "Minimum\nTemperature of\nColdest Month",
                                 #"Precipitation\nof Driest\nMonth",
                                 "Mean Diurnal\nRange",
                                 #"Temperature\nAnnual Range",
                                 "Precipitation\nSeasonality",
                                 "Temperature\nstability",
                                 "Precipitation\nstability",
                                 "Topographic\nheterogeneity",
                                 "Mid-Domain", 
                                 "Spatial"))
unique(dfi$lifeForm)
dfi$lifeForm <- factor(dfi$lifeForm,
                       levels = c("All", "Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Herb"),
                       labels = c("All", "Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Herb"))

#Plot
pblapply(unique(dfi$set), function(x){
  dfi_x <- dfi[dfi$set == x,]
  
  if(x == "pet")
    dfi_x <- dfi_x %>% filter(Variable != "Precipitation\nSeasonality")
  
  if(x == "bio15")
    dfi_x <- dfi_x %>% filter(Variable != "PET")
  
  #Subset colors
  myc_x <- myc[names(myc) %in% unique(dfi_x$Variable)]
  dfi_x$Variable <- droplevels(dfi_x$Variable)
  
  #Plot
  g_imp <- ggplot(data = dfi_x, aes(x = Variable, y = Importance, fill = Variable)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = myc_x) +
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
  
  #Save
  filename <- paste0("Figures/VariableImportance_", x, ".png")
  ggsave(filename,
         g_imp, dpi = 600, units = "px", width = 2500,
         height = 1450, scale = 3)
  })


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

#Read data
lf_indices <- pblapply(lf_indices, readRDS)

#Get lifeforms
lf_names <- sapply(lf_indices, function(x) x$lifeform)
names(lf_indices) <- lf_names
#Reorder lifeforms
lf_names <- c("All", "Tree", "Liana", "Shrub", "Subshrub", "Herb", "Bamboo", "Palm_tree")

#Reorder lf_indices
lf_indices <- lapply(lf_names, function(x) lf_indices[[x]])
names(lf_indices) <- lf_names

#Rasterize indices
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_limite_integrador.gpkg")
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
names(p) <- names(r)
# #Get letters to legend 
# my_legend <- paste0("(", letters[1:6], ") ", lf_names)

#Arrange plot
pp <- wrap_plots(p) + 
  plot_layout(ncol = 3, axes = "collect")

#Save
ggsave("Figures/Richness.png",
       pp, dpi = 600, units = "px", width = 1600,
       height = 1600, scale = 7)

####Predictors####
r_pred <- rast("Data/Variables/Explanatory_Variables.tiff")
plot(r_pred)
#Remove Bio11

#Change names
names(r_pred) <- c("Aridity", "Bio06", "Bio02", "Bio15",
                   "Mid-domain",
                   "Precipitation stability", "Temperature stability",
                   "Topographic heterogeneity",
                   "Spatia filter (Lat)", "Spatial filter (Long)")
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
    coord_sf(xlim = c(bb_af[1] - 0.5, xmax=bb_af[2] + 0.5),
             ylim = c(bb_af[3] - 0.5, ymax=bb_af[4] + 0.5),
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

pp_pred2 <- wrap_plots(p_pred) + plot_layout(ncol = 5)
ggsave("Figures/Predictors.png",
       pp_pred2, dpi = 600, units = "px", width = 2500,
       height = 1200, scale = 7)

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
i <- cm_l[4]

bm <- pblapply(cm_l, function(i){
  cm_i <- readRDS(i) #Read i data
  #Get lifeform and endemism
  lf_i <- unique(cm_i$lifeForm)
  
  #Select best model
  bm <- cm_i %>%
    #filter(Highest_VIF <= 10) %>% #VIF < 10
    #filter(p_dispersion_ration > 0.05) %>% #Dispersion ratio < 0.05
    #calculate delta AIC
    mutate(dAIC = AIC - min(AIC, na.rm = T), .after = AIC) %>% #Recalculate AIC
    filter(dAIC == 0) #AIC < 2
  #Select columns
  bm <- bm %>% dplyr::select(c(lifeForm, Formula, AIC, dAIC,
                               twologlik = loglik, Moran_I = Moran_I_obs, rmse))
  return(bm)
  }) %>% bind_rows()
#Write table
# library(flextable)
# flextable(bm)
write.csv(bm, "Data/Best_models_table.csv", row.names = FALSE)

#Edit table
library(flextable)
bm_df <- bm %>% dplyr::select(-Formula) %>% 
  flextable() %>% 
  align(align = "center", part = "all")
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
lfs <- c("All", "Tree", "Liana", "Shrub", "Subshrub",  "Herb", "Bamboo", "Palm_tree")

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

#### Simplify ecoregions ####
library(mapview)

#Import shapefiles
af <- vect("https://github.com/wevertonbio/spatial_files/raw/main/Data/AF_limite_integrador.gpkg")
eco <- vect("c:/Users/wever/Downloads/Ecoregions2017/Ecoregions2017.shp")
plot(af)
#Crop
eco_af <- crop(eco, af)
plot(eco_af)
eco_af$ECO_NAME %>% unique()
mapview(eco_af, zcol = "ECO_NAME", burst = TRUE)

#Change names
# eco_af$ECO_NAME[eco_af$ECO_NAME == "Atlantic Coast restingas"] <- "Atlantic Coast Restingas"
# eco_af$ECO_NAME[eco_af$ECO_NAME == "Northeast Brazil restingas"] <- "Northeast Restingas"
eco_af$ECO_NAME[eco_af$ECO_NAME == "Caatinga Enclaves moist forests"] <- "Caatinga"

#Remove
to_remove <- c("Chiquitano dry forests", "Pantanal", "Atlantic Coast restingas",
               "Northeast Brazil restingas", "Mangrooves",
               "Fernando de Noronha-Atol das Rocas moist forests", 
               "Humid Chaco",
               "Southern Cone Mesopotamian savanna", 
               "Southern Atlantic Brazilian mangroves")
eco_af <- subset(eco_af, !(eco_af$ECO_NAME %in% to_remove))

#Plot
eco_af

#Convert to raster
r_base <- rast(ext = ext(af), res = 0.008333333, values = 0)
eco_r <- rasterize(eco_af, r_base, field = "ECO_NAME")
plot(eco_r)
#Get levels
l <- levels(eco_r)

#Fill NAs
eco_r2 <- focal(eco_r, 29, "modal", na.policy="only")
levels(eco_r2) <- l
eco_r2 <- crop(eco_r2, af, mask = TRUE)
plot(eco_r2)
#Convert to vector again
eco_v <- as.polygons(eco_r2)
mapview(eco_v) + mapview(eco_af, zcol = "ECO_NAME") + mapview(af)

#Save
writeVector(eco_v, "Data/Ecoregions_af.gpkg",  overwrite = T)



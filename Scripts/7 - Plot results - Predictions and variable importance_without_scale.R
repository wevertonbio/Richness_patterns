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
lp <- list.files("Data/Models/Predictions_without_scale/", full.names = TRUE)
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
unique(d$group)
myc <- c("#D55E00", #Aridity -
         "#E69F00", #Bio06
         "#4DBBD5FF", #Bio14
         "#F0E442", #Bio07
         "#DC0000FF", #Bio07 -
         "#3C5488FF", #Bio15 -
         "#F39B7FFF", #Temp_stab -
         "#56B4E9", #Prec_stab -
         "#CC79A7", #Mid_domain -
         "#7E6148FF", #Topo_het -
         "#999999") #Spatial (EV2) -

#Make factors
d$group <- factor(d$group,
                 levels = c("Aridity",
                               "Bio06", "Bio14",
                               "Bio07", "Bio15",
                               "Prec_stab", "Temp_stab",
                               "Mid_domain",
                               "Topoi_het", "EV2"),
                 labels = c("Aridity",
                            "Minimum Temperature of\nthe Coldest Month",
                            "Precipitation of\nDriest Month",
                            "Temperature Annual\nRange", "Precipitation\nSeasonality",
                            "Precipitation\nstability", "Temperature\nstability",
                            "Mid-domain",
                            "Topographic\nheterogeneity", "Spatial\nfilter"))
d$lifeForm <- factor(d$lifeForm,
                       levels = c("All", "Tree",
                                  "Shrub", "Subshrub",
                                  "Herb", "Liana"))
         
#Plot
gd <- ggplot(d, aes(x = x, y = predicted, colour = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group,
                  colour = NULL), alpha = 0.5) +
  geom_line() +
  scale_colour_manual(values = myc) +
  scale_fill_manual(values = myc) +
  facet_grid(lifeForm ~ group, scales = "free") +
  xlab("Variables (scaled and centered)") + ylab("Predicted richness") + 
  #theme_bw() +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", 
                                        colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        axis.text.x = element_text(angle = 45, hjust=1))
X11()
gd
#Save as PNG file
ggsave("Figures/Predictions_without_scale.png", gd, units = "px",
       dpi = 600, width = 2500,
       height = 1700, scale = 4)



#### Variables importance ####
li <- list.files("Data/Models/Partitioning_without_scale/", full.names = TRUE)
li <- pblapply(li, readRDS)
dfi <- bind_rows(li)

#Remover others
dfi <- dfi %>% filter(lifeForm != "Other")

#Get colors
unique(dfi$Variable)
myc <- c("#D55E00", #Aridity -
         "#E69F00", #Bio06
         "#4DBBD5FF", #Bio14
         "#F0E442", #Bio07
         "#DC0000FF", #Bio07 -
         "#3C5488FF", #Bio15 -
         "#F39B7FFF", #Temp_stab -
         "#56B4E9", #Prec_stab -
         "#CC79A7", #Mid_domain -
         "#7E6148FF", #Topo_het -
         "#999999") #Spatial (EV2) -

#Make factors
dfi$Variable <-factor(dfi$Variable,
                      levels = c("Aridity",
                                 "Bio06", "Bio14",
                                 "Bio07", "Bio15",
                                 "Prec_stab", "Temp_stab",
                                 "Mid_domain",
                                 "Topoi_het", "Spatial"),
                      labels = c("Aridity",
                                 "Bio06", "Bio14",
                                 "Bio07", "Bio15",
                                 "Precipitation\nstability", "Temperature\nstability",
                                 "Mid-domain",
                                 "Topographic\nheterogeneity", "Spatial\nfilter"))
unique(dfi$lifeForm)
dfi$lifeForm <- factor(dfi$lifeForm,
                       levels = c("All", "Tree", "Liana",
                                  "Shrub", "Subshrub",
                                  "Herb"))

g_imp <- ggplot(data = dfi, aes(x = Variable, y = Importance, fill = Variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = myc) +
  facet_wrap(.~lifeForm) +
  ggpubr::theme_pubclean() +
  #theme_stata() +
  #theme_wsj()+ 
  theme(axis.text.x = element_text(angle = 45,
                                   #vjust = 0.5
                                   hjust=1),
        legend.position = "bottom",
        panel.border=element_rect(colour="black",size=1, fill = NA))
g_imp
#Save
ggsave("Figures/VariableImportance_without_scale.png",
       g_imp, dpi = 600, units = "px", width = 2500,
       height = 1700, scale = 4)


####Plot maps####
#Load packages
library(rasterVis)
library(terra)
library(geobr)
library(tidyterra)
library(ggspatial)
library(scales)

#Import vectors
  #South America
sa <- vect("https://github.com/wevertonbio/ENM_Rscripts/raw/main/Vectors/South_America.gpkg")
br <- read_state() %>% vect()

####Richness
#Get richness of lifeforms
lf_indices <- list.files("Data/PAM_indices/", full.names = TRUE)
#Remove others
lf_indices <- lf_indices[!grepl("Other", lf_indices)]
#Read data
lf_indices <- pblapply(lf_indices, readRDS)

#Rasterize indices
af <- vect("https://github.com/wevertonbio/Get_and_Filter_Points/raw/main/Vectors/AF_dissolved..gpkg")
r_base <- rast(ext = ext(af), res = 0.08333333)
r <- pblapply(lf_indices, function(x){
  rx <- rasterize(x$xy, r_base,
            values = x$Richness)
  rx[rx == 0] <- NA
  return(rx)
}) %>% rast()
plot(r[[1]])

#Get box limits to plot
bb_af <- ext(r[[1]])

#Test
#i <- r[[1]]
p <- pblapply(r, function(i) {
  #Convert raster to dataframe
  df_r <- as.data.frame(i, xy = TRUE) %>% na.omit()
  
  g <- ggplot() +
    geom_sf(data = sa, fill = "grey77", size = 0.1, colour = "white") +
    geom_sf(data = br, fill = "grey80", size = 0.1, colour = "grey40") +
    geom_raster(data = df_r, aes(x, y, fill = last), alpha = 1) +
    scale_fill_gradientn(colors = rev(pals::brewer.rdylbu(10)),
                         name = "Richness") +
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
    metR::scale_x_longitude(ticks = 5) + metR::scale_y_latitude(ticks = 6)
  g
})
#Get lifeforms
lf_names <- pbsapply(lf_indices, function(x){ x$lifeform})
#Get letters to legend 
my_legend <- paste0("(", letters[1:6], ") ", lf_names)

#Arrange plot
pp <- (p[[1]] + p[[2]] + p[[3]]) / (p[[4]] + p[[5]] + p[[6]]) +
  plot_layout(tag_level = 'new') +
  plot_annotation(tag_levels = list(my_legend))

#Save
ggsave("Figures/Richness.png",
       pp, dpi = 600, units = "px", width = 2500,
       height = 1700, scale = 5)

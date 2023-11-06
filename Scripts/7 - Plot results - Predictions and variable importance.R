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
        panel.border = element_rect(fill = NA, colour = "grey20"))
X11()
gd
#Save as PNG file
ggsave("Figures/Predictions.png", gd, units = "px",
       dpi = 600, width = 2500,
       height = 1700, scale = 4)

#### Variables importance ####
li <- list.files("Data/Models/Partitioning", full.names = TRUE)
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
ggsave("Figures/VariableImportance.png",
       g_imp, dpi = 600, units = "px", width = 2500,
       height = 1700, scale = 4)

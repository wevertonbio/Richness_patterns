#### Plot results ####
library(dplyr)
library(INLA)
library(terra)
library(tidyterra)
library(ggplot2)
library(data.table)

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
                       "Precipitation\nstability", "Mean diurnal\nrange", 
                       "Min. Temp. of\n Coldest Month", "PET", 
                       "Aridity"))
# Lifeforms
fixed$lifeform <- factor(x = fixed$lifeform,
             levels = c("Tree", "Liana", "Shrub", "Subshrub", 
                        "Terrestrial_herb", "Epiphytic_herb"),
             labels = c("Tree", "Liana", "Shrub", "Subshrub", 
                        "Terrestrial herb", "Epiphytic herb"))

g_fixed <- ggplot(fixed, aes(x = v, y = mean, colour = Effect)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("firebrick", "gray", "forestgreen")) +
  coord_flip() +
  theme_minimal() +
  labs(y = "Estimates", x = "Variables") +
  facet_wrap(lifeform ~ ., scales = "fixed") + 
  ggpubr::theme_pubclean() +
  theme(legend.position = "bottom")
g_fixed

#### Variable importance ####
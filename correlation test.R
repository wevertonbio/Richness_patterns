#### Modified version of the t test to assess the correlation between richness lifeforms ####

# Load packages
library(dplyr)
library(SpatialPack)
library(parallel)
library(data.table)
library(ggcorrplot)
library(tidyr)
library(rstatix)

# Import richness
lf_indices <- readRDS("Data/Richness_by_lifeform.rds")

# Get coordinates
xy <- lf_indices$xy %>% as.matrix()

# Get lifeforms
lf <- c("Tree", "Liana", "Shrub", "Subshrub", "Terrestrial_herb", "Epiphytic_herb")

# Get dataframe with combinations
g <- expand.grid(lf, lf)

# Remove reverse dupicates
g <- g[!duplicated(t(apply(g, 1, sort))), ]


# Looping to calculate correlations
# Make cluster to run in paralell
cl <- parallel::makeCluster(7)
parallel::clusterExport(cl, varlist= c("g", "xy", "lf_indices"),
                        envir = environment()) #Send objects to nodes
parallel::clusterEvalQ(cl, {  #Send packages to nodes
  library(dplyr) #To fit negative binomial
  library(SpatialPack) #To check performance
})

t_cor <- pblapply(1:nrow(g), function(i){
  v1 <- g$Var1[i] %>% as.character()
  v2 <- g$Var2[i] %>% as.character()
  if(v1 == v2){
    data.frame(v1 = v1, v2 = v2, cor = 1, p = 0)
  } else {
    res <- modified.ttest(x = lf_indices[[v1]],
                   y = lf_indices[[v2]],
                   coords = xy)
    data.frame(v1 = v1, v2 = v2, 
               cor = res$corr, 
               p = res$p.value)
  }
}, cl = cl)
stopCluster(cl)

# Merge information
res_cor <- rbindlist(t_cor)

# Convert to correlation matrix
mcor <- res_cor %>%
  select(v1, v2, cor) %>%
  pivot_wider(names_from = v2, values_from = cor) %>%
  tibble::column_to_rownames("v1") %>%
  as.matrix()
# Fill matrix
mcor[upper.tri(mcor)] <- t(mcor)[upper.tri(mcor)]
# Round
mcor <- mcor %>% round(digits = 2)

# Save correlation matrix
saveRDS(mcor, "Data/Richness_correlations.rds")

# Import correlation matrix, if necessary
mcor <- readRDS("Data/Richness_correlations.rds")

# Set colors to highlight cor > 0.9
#Vector of correlations
v_cor <- mcor %>% pull_lower_triangle(diagonal = FALSE) %>% as.matrix() %>%
  #Get the lower part
  as.numeric() %>% abs() %>% na.omit()

#Convert to color
v_col <- ifelse(v_cor >= 0.9, "black", 
                ifelse(v_cor < 0.9 & v_cor >= 0.75, "gray20", "gray30"))
outline_color <- ifelse(v_cor >= 0.9, "black", "gray")
lab_size <- ifelse(v_cor >= 0.9, 3.5, 2.75)

# Rename life forms
row.names(mcor) <- gsub("_", " ", row.names(mcor))
colnames(mcor) <- gsub("_", " ", colnames(mcor))

# Plotar
g <- ggcorrplot(mcor, method = "square", type = "lower",
                lab = TRUE, lab_col = v_col, 
                #outline.color = outline_color,
                lab_size = lab_size, hc.order = TRUE,
                ggtheme = ggpubr::theme_pubclean()) +
  scale_fill_gradient2(low = "#6D9EC1", high = "#E46726", mid = "white",
                       midpoint = 0.75, limit = c(0.5, 1), space = "Lab",
                       name="Correlation") +
  theme(legend.position = "right", 
        legend.key.height = unit(1.75, 'cm'))
g
# Save
dir.create("Figures")
ggsave("Figures/Richness_Correlations.png",
       g, dpi = 600, units = "px", width = 2500,
       height = 2000, scale = 1.5)

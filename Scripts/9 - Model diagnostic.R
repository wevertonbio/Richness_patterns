#Model diagnostic#
library(dplyr)
library(DHARMa)
library(performance)

#Import best models
bml_list <- list.files("Data/Models/Best_models/", full.names = T)
bml <- pblapply(bml_list, readRDS)
names(bml) <- bml_list %>% basename() %>% fs::path_ext_remove()

#Check overdispersion and outliers
od <- pblapply(names(bml), function(x){
  over_x <- performance::check_overdispersion(bml[[x]], alternative = "greater")
  under_x <- performance::check_overdispersion(bml[[x]], alternative = "less")
  out_x <- performance::check_outliers(bml[[x]])
  
  data.frame(lifeform = x,
             overdispersion = over_x$dispersion_ratio,
             overdispersion_p = over_x$p_value,
             underdispersion = under_x$dispersion_ratio,
             underdispersion_p = under_x$p_value,
             outliers = length(which(out_x)))
})
od <- bind_rows(od)


#QQPlot
png(filename = "Figures/QQplot.png", units = "in", width = 8, height = 6, res = 300)
par(mfrow = c(2, 3)) # Create a 2 x 2 plotting matrix
pblapply(names(bml), function(x){
  DHARMa::plotQQunif(bml[[x]], testOutliers = FALSE, testDispersion = FALSE, main = x)
  })
dev.off()
dev.off()


# DHARMa::plotQQunif(sr, testOutliers = FALSE, testDispersion = FALSE, main = "Tree")
# DHARMa::plotQQunif(sr, testOutliers = FALSE, testDispersion = FALSE, main = "Tree")
# DHARMa::plotQQunif(sr, testOutliers = FALSE, testDispersion = FALSE, main = "Tree")
# DHARMa::plotQQunif(sr, testOutliers = FALSE, testDispersion = FALSE, main = "Tree")




####Generate hypothesis graph####
dir.create("Data/Hypothesis")

library(ggplot2)


# Create simulated data for intersecting lines
set.seed(123)  # Set a seed for reproducibility

# Coefficients of the models
beta_positive <- 1
beta_negative <- -1
intercept_positive <- 50 - beta_positive * 50  # Adjusted intercept to cross at 50
intercept_negative <- 50 - beta_negative * 50  # Adjusted intercept to cross at 50

data <- data.frame(
  X = 1:100,
  Y_positive = beta_positive * (1:100) + intercept_positive,
  Y_negative = beta_negative * (1:100) + intercept_negative)


#Negative relations
  #Energy
g.energy<- ggplot(data, aes(x = X, y = Y_negative)) +
  geom_line(size = 3, colour = "firebrick") +
  xlab("Aridity") + ylab("Richness") +
  theme_classic() +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50)) +
  annotate("text", x = min(data$X), y = -2, label = "More energy", size = 13, hjust = 0) +
  annotate("text", x = max(data$X), y = -2, label = "Less energy", size = 13, hjust = 1)

g.energy
ggsave("Data/Hypothesis/Energy.png", g.energy, dpi = 300, height = 8)

#Stability
g.stab <- ggplot(data, aes(x = X, y = Y_positive)) +
  geom_line(size = 3, colour = "forestgreen") +
  xlab("Stability in the last 21ky") + ylab("Richness") +
  theme_classic() +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50)) +
  annotate("text", x = min(data$X), y = -2, label = "Less stable", size = 13, hjust = 0) +
  annotate("text", x = max(data$X), y = -2, label = "More stable", size = 13, hjust = 1)

g.stab
ggsave("Data/Hypothesis/Stability.png", g.stab, dpi = 300, height = 8.5)

  #Mid-Domain
g.mid <- ggplot(data, aes(x = X, y = Y_negative)) +
  geom_line(size = 3, colour = "firebrick") +
  xlab("Distance to centroid (Latitude)") + ylab("Richness") +
  theme_classic() +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50))

g.mid
ggsave("Data/Hypothesis/MID.png", g.mid, dpi = 300, height = 8.5)

  #Seasonality (positive and negative)
g.seas <- ggplot(data, aes(x = X, y = Y_negative)) +
  geom_line(size = 3, colour = "firebrick") +
  geom_line(data =data, aes(x = X, y = Y_positive), size = 3, colour = "forestgreen") +
  xlab("Seasonality") + ylab("Richness") +
  theme_classic() +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50))

g.seas
ggsave("Data/Hypothesis/Seasonality.png", g.seas, dpi = 300, height = 8.5)

#Seasonality (only negative)
g.seas2 <- ggplot(data, aes(x = X, y = Y_negative)) +
  geom_line(size = 3, colour = "firebrick") +
  xlab("Seasonality") + ylab("Richness") +
  theme_classic() +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50))

g.seas2
ggsave("Data/Hypothesis/Seasonality2.png", g.seas2, dpi = 300, height = 8.5)

#Tolerance (positive and negative)
g.tol <- ggplot(data, aes(x = X, y = Y_negative)) +
  geom_line(size = 3, colour = "firebrick") +
  geom_line(data =data, aes(x = X, y = Y_positive), size = 3, colour = "forestgreen") +
  xlab("Extreme of temperature\nor precipitation") + ylab("Richness") +
  theme_classic() +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50))

g.tol
ggsave("Data/Hypothesis/Tolerance.png", g.tol, dpi = 300, height = 8.5)

#Tolerance (Positive)
g.tol2 <- ggplot(data, aes(x = X, y = Y_positive)) +
  geom_line(size = 3, colour = "forestgreen") +
  xlab("Extremes of temperature\nor precipitation") + ylab("Richness") +
  theme_classic() +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50)) +
  annotate("text", x = min(data$X), y = -2, label = "Colder/Dryer", size = 13, hjust = 0) +
  annotate("text", x = max(data$X), y = -2, label = "Warmer/Wetter", size = 13, hjust = 1)

g.tol2
ggsave("Data/Hypothesis/Tolerance2.png", g.tol2, dpi = 300, height = 8.5)

#Topographic 
g.topo <- ggplot(data, aes(x = X, y = Y_positive)) +
  geom_line(size = 3, colour = "forestgreen") +
  xlab("Diversity of habitats") + ylab("Richness") +
  theme_classic() +
  theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.line.x = element_line(arrow = grid::arrow(length = unit(0.75, "cm"), 
                                                       ends = "last", type = "closed"),
                                   linewidth = 2),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 50),
        axis.title.x = element_text(size = 50))

g.topo
ggsave("Data/Hypothesis/topographic.png", g.topo, dpi = 300, height = 8.5)

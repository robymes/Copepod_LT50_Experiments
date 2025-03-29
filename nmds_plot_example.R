# SIMULATE DATA

# Load packages 
library(vegan)      # for NMDS functions
library(tidyverse)  # for data wrangling and ggplot2
library(goeveg)     # for NMDS scree plot

# Reproducibility  
set.seed(1)  

# Create "site1" to "site20"
site_names <- paste0("site", 1:20) 

# Create "species1" to "species8" and "Pine Forest" and "Oak Forest" 
species_names <- paste0("species", 1:8) 
habitat_names <- rep(c("PineForest", "OakForest"), each = 80)  

# Create data frame 
df <- 
  data.frame(
    Site = rep(site_names, each = 8),
    Species = rep(species_names, times = 1),
    Habitat = rep(c("PineForest","OakForest"), each = 80))

# Create column of species abundances. Sample counts of 0:51 for the Oak Forest
# and 0:10 for the Pine Forest and get a unique count for each species at each
# site using `n()`.
df <- df %>% 
  mutate(Abundance = ifelse(
    Habitat == "PineForest", 
    ifelse(Species == "species1" | Species == "species5", 0, 
           sample(0:51, n(), replace = TRUE)),
    sample(0:10, n(), replace = TRUE)))
head(df)

#-------------------------------------------------------------------------------

# STACKED BAR PLOT

ggplot(df, aes(x = Habitat, y = Abundance, fill = Species)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic()

#-------------------------------------------------------------------------------

# PIVOT DATA FRAME

com_matrix <- df %>% 
  
  # Turn our Site, Species, and Forest columns into factors 
  mutate_at(vars(Site), as.factor) %>% 
  
  # pivot to wide format
  pivot_wider(names_from = Species,
              values_from = Abundance) %>% 
  
  # change our column "site" to our rownames
  column_to_rownames(var = "Site")

head(com_matrix)

#-------------------------------------------------------------------------------

# RUN nMDS ANALYSIS USING BRAY-CURTIS DISTANCE

nmds <-
  metaMDS(com_matrix[,2:9],
          distance = "bray",
          k = 2)


print(nmds)

# How to interpret the "Stress" measure:
# < 0.05 excellent fit of nMDS
# < 0.01 good fit with little risk of misinterpretation
# < 0.02 fair to use but values close to 0.02 are at risk of misinterpretation
# > 0.02 poorly represents the data

#-------------------------------------------------------------------------------

# CREATE STRESS PLOT

# Shows how Stress is calculated
stressplot(nmds)

#-------------------------------------------------------------------------------

# CREATE SCREE PLOT

# Usually Stress decreases with increasing number of dimensions
# However, the fewer dimensions we use, the easier our data is to interpret
dimcheck_out <- 
  dimcheckMDS(com_matrix[,2:9],
              distance = "bray",
              k = 6)

# Look out for number of dimensions after which you get "diminishing returns" in
# the decrease of Stress
print(dimcheck_out)

#-------------------------------------------------------------------------------

# CREATE nMDS PLOT

# red crosses are species, black points are sites
plot(nmds)

# get habitat type: which points come from which habitat?
habitat_type <- df %>% 
  distinct(Site, Habitat)

# Extract NMDS scores for sites 
nmds_SiteScores <-
  # get nmds scores 
  as.data.frame(scores(nmds)$sites) %>%
  # change rownames (site) to a column 
  rownames_to_column(var = "Site") %>% 
  # join our habitat type (grouping variable) to each site 
  left_join(habitat_type)

# Extract NMDS scores for species  
nmds_SpeciesScores <- 
  as.data.frame(scores(nmds, "species"))

# create a column of species, from the rownames of species.scores
nmds_SpeciesScores$species <- rownames(nmds_SpeciesScores)

# create updated plot
ggplot() + 
  
  # add site scores
  geom_point(data = nmds_SiteScores, 
             aes(x=NMDS1, y=NMDS2, colour = Habitat), 
             size = 2) + 
  
  # add species scores 
  geom_text(data = nmds_SpeciesScores, 
            aes(x=NMDS1, y=NMDS2, label = species)) +
  
  theme_classic()

#-------------------------------------------------------------------------------

# IMPROVE PREVIOUS PLOT

# get centroid (central value for each habitat)
Habitat_Centroid <- 
  nmds_SiteScores %>% 
  group_by(Habitat) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

# extract convex hull (border around points belonging to the same habitat)
habitat.hull <- 
  nmds_SiteScores %>% 
  group_by(Habitat) %>%
  slice(chull(NMDS1, NMDS2))

# extract stress value to print it on the plot (VERY IMPORTANT TO SHOW)
nmds_stress <- nmds$stress

#-------------------------------------------------------------------------------

# CREATE IMPROVED PLOT

# use ggplot to plot 
ggplot() + 
  
  # add site scores
  geom_point(data = nmds_SiteScores, 
             aes(x=NMDS1, y=NMDS2, colour = Habitat), size = 2) + 
  
  # add species scores 
  geom_text(data = nmds_SpeciesScores, 
            aes(x=NMDS1, y=NMDS2, label = species)) +
  
  # add centroid 
  geom_point(data = Habitat_Centroid, 
             aes(x = axis1, y = axis2, color = Habitat), 
             size = 5, shape = 17) +
  
  # add convex hull
  geom_polygon(data = habitat.hull, 
               aes(x = NMDS1, y = NMDS2, fill = Habitat, group = Habitat), 
               alpha = 0.30) +
  
  # add stress value
  annotate("text", x = 0.75, y = 0.65, 
           label = paste("2d stress =", round(nmds_stress, 3))) +
  
  # edit theme (axis labels are meaningless here, so they are removed)
  labs(x = "NMDS1", y = "NMDS2") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, linewidth = .5),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(.25, "cm"))

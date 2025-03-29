# Loading required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Creating data frames for the two habitats
habitat1 <- data.frame(
  Habitat = "Tubeworms & Mussels",
  Species = c("Other", "Rhogobius sp.", "Nilva torifera", "Stygiopontius hispidulus", 
              "Aphotopontius mammillatus", "Ceuthoecetes acanthothrix"),
  Percentage = c(9.33, 5.60, 6.15, 10.98, 31.94, 36.00)
)

habitat2 <- data.frame(
  Habitat = "Pompeii Worms",
  Species = c("Other", "Stygiopontius hispidulus"),
  Percentage = c(8.49, 91.51)
)

# Merging the two data frames
complete_data <- rbind(habitat1, habitat2)

# Viewing the created data frame
print(complete_data)

# Custom color palette for species
species_colors <- c("Other" = "#999999", 
                    "Rhogobius sp." = "#E69F00", 
                    "Nilva torifera" = "#56B4E9", 
                    "Stygiopontius hispidulus" = "#009E73", 
                    "Aphotopontius mammillatus" = "#F0E442", 
                    "Ceuthoecetes acanthothrix" = "#0072B2")

# Creating the stacked column chart
# First, create a named vector for italicized labels
species_labels <- c(
  "Other" = "Other",
  "Rhogobius sp." = expression(italic("Rhogobius")~"sp."),
  "Nilva torifera" = expression(italic("Nilva torifera")),
  "Stygiopontius hispidulus" = expression(italic("Stygiopontius hispidulus")),
  "Aphotopontius mammillatus" = expression(italic("Aphotopontius mammillatus")),
  "Ceuthoecetes acanthothrix" = expression(italic("Ceuthoecetes acanthothrix"))
)

chart <- ggplot(complete_data, aes(x = Habitat, y = Percentage, fill = Species)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = species_colors, 
                    labels = species_labels) +
  labs(title = "Species Distribution by Habitat",
       x = "Habitat",
       y = "Population Percentage (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text.align = 0  # 0 = left alignment, 1 = right alignment
  )

# Displaying the chart
print(chart)

# Saving the chart (optional)
# ggsave("habitat_species_distribution.png", chart, width = 10, height = 6, dpi = 300)
library(grid)
library(tidyverse)
library(shadowtext)

#------------------------------------------------------------------------------

#Creation of Data Frame

names <- c(
  "Other", "Stygiopontius hispidulus"
)

# Name is an ordered factor. We do this to ensure the bars are sorted.

data <- data.frame(
  count = c(8.49, 91.51), 
  name = factor(names, levels = names),
  y = seq(length(names)) * 0.9,
  percentages = c("8.49%", "91.51%")
)

#------------------------------------------------------------------------------

# The colors

BLUE <- "#076fa2"
RED <- "#E3120B"
BLACK <- "#202020"
GREY <- "grey50"

#------------------------------------------------------------------------------

#Basic Bar Chart

plt <- ggplot(data) +
  geom_col(aes(count, name), fill = RED, width = 0.6) 

plt

#------------------------------------------------------------------------------

#Customized Layout

plt <- plt + 
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 90, by = 5), 
    expand = c(0, 0), # The horizontal axis does not extend to either side
    position = "top"  # Labels are located on the top
  ) +
  # The vertical axis only extends upwards 
  scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    axis.title = element_blank(),
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank(),
    # But customize labels for the horizontal axis
    axis.text.x = element_text(family = "Econ Sans Cnd", size = 16)
  )

plt

#------------------------------------------------------------------------------

#Add Labels

plt <- plt + 
  geom_shadowtext(
    data = subset(data, count >= 5),
    aes(count, y = name, label = percentages),
    hjust = 0,
    nudge_x = 0.3,
    colour = RED,
    bg.colour = "white",
    bg.r = 0.2,
    family = "Aptos",
    size = 6.5
  ) + 
  geom_text(
    data = subset(data, count >= 5),
    aes(0, y = name, label = name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Aptos",
    fontface = "italic",
    size = 6.5
  )

plt

#------------------------------------------------------------------------------

#Add Annotations

plt <- plt +
  labs(
    title = "Pompeii worm Habitat",
    subtitle = "Dirivultid Copepod Community Diversity"
  ) + 
  theme(
    plot.title = element_text(
      family = "Econ Sans Cnd", 
      face = "bold",
      size = 22
    ),
    plot.subtitle = element_text(
      family = "Econ Sans Cnd",
      size = 20
    )
  )
plt

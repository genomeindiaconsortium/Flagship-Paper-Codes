## -----------------------------
## Libraries
## -----------------------------
library(plotly)
library(dplyr)
library(tidyr)

## -----------------------------
## Input data
## -----------------------------
df <- data.frame(
  Group = c("Tribe", "Tribe", "Tribe", "Tribe",
            "NonTribe", "NonTribe", "NonTribe", "CAO"),
  Subgroup = c("IE", "DR", "TB", "AA",
               "IE", "DR", "TB", "CAO"),
  PopSpecific = c(10, 26, 9, 10, 54, 24, 1, 30),
  Novel = c(1, 10, 3, 2, 5, 3, 1, 1)
)

## -----------------------------
## Data transformation
## -----------------------------
df_long <- df %>%
  pivot_longer(
    cols = c(PopSpecific, Novel),
    names_to = "VariantType",
    values_to = "Count"
  ) %>%
  filter(Count > 0) %>%
  mutate(
    group_id    = Group,
    subgroup_id = paste(Group, Subgroup, sep = "_"),
    variant_id  = paste(subgroup_id, VariantType, sep = "_")
  )

## -----------------------------
## Build hierarchy (IDs stay intact)
## -----------------------------
level1 <- df_long %>%
  distinct(id = group_id) %>%
  mutate(parent = "", value = 0, hover = "", label = id)

level2 <- df_long %>%
  distinct(id = subgroup_id, parent = group_id) %>%
  mutate(value = 0, hover = "", label = sub(".*_", "", id))

level3 <- df_long %>%
  transmute(
    id     = variant_id,
    parent = subgroup_id,
    value  = Count,
    hover  = paste0(VariantType, ": ", Count),
    label  = as.character(Count)  # numbers only in outer ring
  )

sunburst_data <- bind_rows(level1, level2, level3)

## -----------------------------
## Color mapping
## -----------------------------
subgroup_colors <- c(
  "Tribe_DR"     = "#000080",
  "Tribe_TB"     = "#7ccd7c",
  "Tribe_AA"     = "#4682b4",
  "Tribe_IE"     = "#cd3333",
  "NonTribe_DR"  = "#8A8AF4",
  "NonTribe_TB"  = "#C7F2C7",
  "NonTribe_IE"  = "#FAC4C4",
  "CAO_CAO"      = "#4c4c4c"
)

sunburst_data$color <- NA

# Leaf colors
sunburst_data$color[sunburst_data$value > 0 &
                      grepl("PopSpecific$", sunburst_data$id)] <- "#999999"
sunburst_data$color[sunburst_data$value > 0 &
                      grepl("Novel$", sunburst_data$id)] <- "gold2"

# Subgroup colors
sunburst_data$color[is.na(sunburst_data$color)] <-
  subgroup_colors[sunburst_data$id]

# Group colors
sunburst_data$color[sunburst_data$id == "Tribe"]    <- "#f7efd2"
sunburst_data$color[sunburst_data$id == "NonTribe"] <- "#dcf0fa"
sunburst_data$color[sunburst_data$id == "CAO"]      <- "#404040"

## -----------------------------
## Plot (CORRECT WAY)
## -----------------------------
fig <- plot_ly(
  type    = "sunburst",
  ids     = sunburst_data$id,      # stable internal IDs
  labels  = sunburst_data$label,   # what is displayed
  parents = sunburst_data$parent,
  values  = sunburst_data$value,
  textinfo = "label",
  hovertext = sunburst_data$hover,
  hoverinfo = "text",
  marker = list(
    colors = sunburst_data$color,
    line = list(color = "white", width = 1)
  ),
  insidetextorientation = "radial"
)

fig$x$layout$sunburst <- list(
  domain = list(x = c(0.25, 0.75), y = c(0.20, 0.80))
)

fig


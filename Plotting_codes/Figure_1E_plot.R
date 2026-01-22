library(readxl)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(scales)

# ---- 1. Read input ----
data_raw1 <- read_excel(
  "PlotInputs.xlsx",
  sheet = "Figure1D"
)
data_raw <- data_raw1 %>% 
  filter(Clusters != "TB_NT")

# ---- 2. Cluster colors ----
cluster_colors <- c(
  "IE_T"  = "#cd3333",
  "DR_T"  = "#000080",
  "AA_T"  = "#4682b4",
  "TB_T"  = "#7ccd7c",
  "IE_NT" = "#FAC4C4",
  "DR_NT" = "#8A8AF4",
  "TB_NT" = "#C7F2C7"
)

# ---- 3. Clean & sort ----
data_clean <- data_raw %>%
  filter(!is.na(Clusters)) %>%
  mutate(
    sample_rank            = as.integer(sample_rank),
    mean_cumulative_genome = as.numeric(mean_cumulative_genome),
    sd_cumulative_genome   = as.numeric(sd_cumulative_genome)
  ) %>%
  arrange(Clusters, sample_rank)

# ---- 4. Parameters ----
MAX_SAMPLES <- 500      
N_BOOT      <- 50       
set.seed(123)       

# ---- 5. Generate bootstrap curves ----
boot_list <- vector("list", N_BOOT)

for (b in seq_len(N_BOOT)) {
  boot_list[[b]] <- data_clean %>%
    group_by(Clusters) %>%
    arrange(sample_rank, .by_group = TRUE) %>%
    
    slice_head(n = MAX_SAMPLES) %>%
    
    mutate(
      Cum_boot = rnorm(
        n(),
        mean = mean_cumulative_genome,
        sd   = sd_cumulative_genome
      ),
      Cum_boot = pmax(Cum_boot, 0)         
    ) %>%
    ungroup() %>%
    transmute(
      Cluster    = Clusters,
      Boot       = b,
      Step       = sample_rank,
      Cumulative = Cum_boot
    )
}

curves_all <- bind_rows(boot_list)

curve_summary <- curves_all %>%
  group_by(Cluster, Step) %>%
  summarise(
    mean_unique = mean(Cumulative, na.rm = TRUE),
    lower       = quantile(Cumulative, 0.025, na.rm = TRUE),
    upper       = quantile(Cumulative, 0.975, na.rm = TRUE),
    .groups     = "drop"
  )

p <- ggplot(curve_summary,
            aes(x = Step, y = mean_unique,
                color = Cluster, fill = Cluster)) +
  geom_line(linewidth = 3.5, alpha = 0.95) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  labs(
    x     = "Number of samples",
    y     = "Cumulative novel variants discovered",
    color = "Cluster",
    fill  = "Cluster"
  ) +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(size = 20, colour = "black", hjust = 1),
    axis.text.y   = element_text(size = 24, colour = "black"),
    axis.title.x  = element_text(size = 24, face = "bold"),
    axis.title.y  = element_text(size = 24, face = "bold"),
    legend.title  = element_text(size = 24, face = "bold"),
    legend.text   = element_text(size = 20, colour = "black"),
    plot.title    = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.position = "", #bottom
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black"),
    legend.justification = "center"
  )

p

ggsave(
  "Figure1E.png",
  p,
  width = 15,
  height = 8,
  dpi = 350
)

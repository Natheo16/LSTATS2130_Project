# Q1

library(ggplot2)
library(tidyr)
library(dplyr)

data <- read.table("Weight.txt", header = TRUE)

# Measurement times in days
t.obs <- c(0, 7, 14, 30*(1:6), 30*seq(12, 30, by=6))

# Rename columns to actual time values
colnames(data) <- c("Boy", t.obs)

# Add subject ID
data$id <- 1:nrow(data)

# Reshape to usable format
data_long <- data %>%
  pivot_longer(cols = -c(id, Boy), names_to = "time", values_to = "weight") %>%
  mutate(time = as.numeric(time),
         sex  = ifelse(Boy == 1, "Boys", "Girls"))

# --- Plot 1: Individual trajectories ---
p1 <- ggplot(data_long, aes(x = time, y = weight, group = id)) +
  geom_line(alpha = 0.3, linewidth = 0.4) +
  facet_wrap(~sex) +
  labs(x = "Age (days)", y = "Weight (kg)",
       title = "Individual weight trajectories by sex") +
  theme_bw()

ggsave("Q1_trajectories.png", p1, width = 10, height = 5, dpi = 150)

# --- Plot 2: Mean weight over time ---
mean_data <- data_long %>%
  group_by(sex, time) %>%
  summarise(mean_weight = mean(weight), .groups = "drop")

p2 <- ggplot(mean_data, aes(x = time, y = mean_weight, color = sex)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Boys" = "steelblue", "Girls" = "tomato")) +
  labs(x = "Age (days)", y = "Mean weight (kg)",
       title = "Mean weight over time by sex", color = "Sex") +
  theme_bw()

ggsave("Q1_mean_growth.png", p2, width = 7, height = 5, dpi = 150)

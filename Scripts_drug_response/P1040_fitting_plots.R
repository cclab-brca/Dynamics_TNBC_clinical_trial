# Load necessary libraries
library(dplyr)
library(ggplot2)
library(lmerTest)
library(tidyr)
library(emmeans)
library(purrr)

# Prepare the data
P1040 <- P1040 %>%
  dplyr::select(Model:Day, Volume, Perc_vol_change) %>%
  rename(Time = Day, Size = Volume) %>%
  filter(Group != "") %>%
  mutate(Size = as.numeric(Size),
         Group = factor(Group),
         lSize = log(Size)) %>%
  drop_na(Size) %>%
  filter(Size > 0, is.finite(lSize)) %>%
  mutate(Groups = paste(Group, Cohort, sep = " -> "),
         Cohort = factor(Cohort, levels = c("Untreated", "CT", "Olaparib")))

# Define plot function
create_plot <- function(data, title) {
  ggplot(data, aes(x = Time, y = Size, group = Mouse, color = Cohort)) +
    geom_line() +
    facet_wrap(~ Cohort, ncol = 1, strip.position = "right") +
    labs(x = "Days since start of treatment", y = "Tumour volume (mm³)", title = title) +
    theme_minimal() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", angle = -90, hjust = 0),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none") +
    scale_color_manual(values = c("Untreated" = "black", "CT" = "#C43826", "Olaparib" = "#858FB0"))
}

# Subset the data by previous treatment
Untreated_PreviousTreatment <- P1040 %>% filter(Previous_Treatment == "Untreated")
CT_Previous_Treatment <- P1040 %>% filter(Previous_Treatment == "CT")
Olaparib_Previous_Treatment <- P1040 %>% filter(Previous_Treatment == "Olaparib")

# Generate plots for each treatment group
plot_Untreated <- create_plot(Untreated_PreviousTreatment, "Untreated Previous Treatment")
plot_CT <- create_plot(CT_Previous_Treatment, "CT Previous Treatment")
plot_Olaparib <- create_plot(Olaparib_Previous_Treatment, "Olaparib Previous Treatment")

# Print the plots
print(plot_Untreated)
print(plot_CT)
print(plot_Olaparib)

# Fit linear mixed-effects models
lmer_Untreated <- lmerTest::lmer(lSize ~ Time * Cohort + (Time | Mouse), data = Untreated_PreviousTreatment, control = lmerControl(optimizer = "bobyqa"))
lmer_CT <- lmerTest::lmer(lSize ~ Time * Cohort + (Time | Mouse), data = CT_Previous_Treatment, control = lmerControl(optimizer = "bobyqa"))
lmer_Olaparib <- lmerTest::lmer(lSize ~ Time * Cohort + (Time | Mouse), data = Olaparib_Previous_Treatment, control = lmerControl(optimizer = "bobyqa"))

# Function to create prediction data
create_pred_data <- function(model, data) {
  new_data <- expand.grid(Time = seq(min(data$Time), max(data$Time), length.out = 100), Cohort = unique(data$Cohort))
  new_data$Predicted <- predict(model, newdata = new_data, re.form = NA)
  new_data
}

# Generate predictions for each treatment group
pred_untreated <- create_pred_data(lmer_Untreated, Untreated_PreviousTreatment)
pred_ct <- create_pred_data(lmer_CT, CT_Previous_Treatment)
pred_olaparib <- create_pred_data(lmer_Olaparib, Olaparib_Previous_Treatment)

# Combine predictions
all_pred_data <- bind_rows(pred_untreated, pred_ct, pred_olaparib, .id = "Previous_Treatment")

# Plot model predictions
ggplot(all_pred_data, aes(x = Time, y = exp(Predicted), color = Cohort)) +
  geom_line() +
  facet_wrap(~ Previous_Treatment, scales = "free_y", ncol = 3) +
  labs(x = "Days", y = "Predicted tumour size", title = "Tumour growth by previous and current treatments") +
  theme_bw() +
  scale_y_log10() +
  scale_color_manual(values = c("Untreated" = "black", "CT" = "#C43826", "Olaparib" = "#858FB0")) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
  )

# Create mouse-level prediction data
create_mouse_pred_data <- function(model, data) {
  mice <- unique(data$Mouse)
  map_dfr(mice, function(mouse) {
    new_data <- data.frame(Time = seq(min(data$Time), max(data$Time), length.out = 100), Mouse = mouse, Cohort = unique(data$Cohort[data$Mouse == mouse]))
    new_data$Predicted <- predict(model, newdata = new_data, re.form = NULL)
    new_data$Previous_Treatment <- unique(data$Previous_Treatment)
    new_data
  })
}

# Generate mouse-level predictions for each treatment group
pred_untreated <- create_mouse_pred_data(lmer_Untreated, Untreated_PreviousTreatment)
pred_ct <- create_mouse_pred_data(lmer_CT, CT_Previous_Treatment)
pred_olaparib <- create_mouse_pred_data(lmer_Olaparib, Olaparib_Previous_Treatment)

# Combine all prediction data
all_pred_data <- bind_rows(pred_untreated, pred_ct, pred_olaparib)
all_actual_data <- bind_rows(Untreated_PreviousTreatment, CT_Previous_Treatment, Olaparib_Previous_Treatment)

# Remove predicted values greater than 1700
all_pred_data <- all_pred_data %>%
  mutate(Predicted_exp = exp(Predicted)) %>%
  filter(Predicted_exp <= 1700)

# Plot individual mouse predictions
ggplot() +
  geom_point(data = all_actual_data, aes(x = Time, y = Size, group = Mouse, color = Cohort), alpha = 0.7) +
  geom_line(data = all_pred_data, aes(x = Time, y = Predicted_exp, group = Mouse, color = Cohort), linetype = "dashed", alpha = 0.5) +
  facet_grid(Cohort ~ Previous_Treatment, scales = "free_y") +
  labs(x = "Days", y = "Tumour volume", title = "Individual mouse tumour growth by previous and current treatments") +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 1700), labels = scales::comma, limits = c(0, 1700)) +
  scale_color_manual(values = c("Untreated" = "black", "CT" = "#C43826", "Olaparib" = "#858FB0")) +
  coord_cartesian(xlim = c(0, 90)) +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )
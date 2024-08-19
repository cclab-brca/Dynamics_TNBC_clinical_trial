# Load necessary libraries
library(dplyr)
library(ggplot2)
library(lmerTest)
library(tidyr)
library(emmeans)
library(lme4)

# Interaction between Untreated post-treatment and other groups
fixed_effects_summary <- data.frame(
  Previous_Treatment = character(),
  Parameter = character(),
  Estimate = numeric(),
  StdError = numeric(),
  Pr = numeric(),
  CI_low = numeric(),
  CI_high = numeric(),
  stringsAsFactors = FALSE
)

subset_models <- list(
  Untreated = list(data = Untreated_PreviousTreatment, model = lmer_Untreated),
  CT = list(data = CT_Previous_Treatment, model = lmer_CT),
  Olaparib = list(data = Olaparib_Previous_Treatment, model = lmer_Olaparib)
)

for (prev_treatment in names(subset_models)) {
  m_t <- subset_models[[prev_treatment]]$model
  summary_m_t <- summary(m_t)
  
  # Extract fixed effects and confidence intervals
  fixed_eff_df <- as.data.frame(summary_m_t$coefficients)
  fixed_eff_df$Parameter <- rownames(summary_m_t$coefficients)
  fixed_eff_df$Previous_Treatment <- prev_treatment
  
  ci <- confint(m_t, level = 0.95)
  
  fixed_eff_df$CI_low <- sapply(fixed_eff_df$Parameter, function(param) ci[param, 1])
  fixed_eff_df$CI_high <- sapply(fixed_eff_df$Parameter, function(param) ci[param, 2])
  
  fixed_effects_summary <- rbind(
    fixed_effects_summary, 
    fixed_eff_df[, c("Previous_Treatment", "Parameter", "Estimate", "Std. Error", "Pr(>|t|)", "CI_low", "CI_high")]
  )
}

fixed_effects_summary <- fixed_effects_summary %>%
  mutate(Model = case_when(
    Previous_Treatment == "Untreated" ~ "lmer_Untreated",
    Previous_Treatment == "CT" ~ "lmer_CT",
    Previous_Treatment == "Olaparib" ~ "lmer_Olaparib"
  ))

# Save the fixed effects summary
write.csv2(fixed_effects_summary, file = "P1040_revisions_M1.csv", row.names = FALSE)

# Visualization of interaction effects
interaction_data <- fixed_effects_summary %>%
  filter(grepl("Time:Cohort", Parameter)) %>%
  separate(Parameter, into = c("Time", "Cohort"), sep = ":") %>%
  mutate(Cohort = gsub("Cohort", "", Cohort),
         Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")))

# Plot option A
ggplot(interaction_data, aes(x = Cohort, y = Estimate, color = Cohort)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  facet_wrap(~ Previous_Treatment, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("CT" = "#C43826", "Olaparib" = "#858FB0")) +
  labs(x = "Current Treatment", y = "Change in growth rate compared to Untreated") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Plot option B
pdf("M1_P1040.pdf", width = 12, height = 9)
ggplot(interaction_data, aes(x = Cohort, y = Estimate, color = Cohort)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2, position = position_dodge(0.5)) +
  facet_grid(Cohort ~ Previous_Treatment, scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("CT" = "#C43826", "Olaparib" = "#858FB0")) +
  labs(x = "Current Treatment", y = "Change in growth rate",
       title = "Change in Growth Rate by Previous Treatment and Current Treatment") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
  )
dev.off()

# AUC between Untreated post-treatment and other groups
emm_pairs_summary <- data.frame()

for (prev_treatment in names(model_list)) {
  m_t <- model_list[[prev_treatment]]
  emm <- emmeans(m_t, specs = ~ Cohort, at = list(Time = c(0, 77)))
  emm_pairs <- pairs(emm, scale = 77)
  
  emm_pairs_df <- summary(emm_pairs, infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    mutate(Previous_Treatment = prev_treatment)
  
  emm_pairs_summary <- rbind(emm_pairs_summary, emm_pairs_df)
}

emm_pairs_summary <- emm_pairs_summary %>%
  mutate(Model = case_when(
    Previous_Treatment == "Untreated" ~ "lmer_Untreated",
    Previous_Treatment == "CT" ~ "lmer_CT",
    Previous_Treatment == "Olaparib" ~ "lmer_Olaparib"
  ))

# Save the estimated marginal means pairwise comparisons
write.csv2(emm_pairs_summary, file = "P1040_revisions_M2.csv", row.names = FALSE)

# Visualization of AUC differences
plot_data <- emm_pairs_summary %>%
  filter(grepl("Untreated -", contrast)) %>%
  mutate(Drug = case_when(
    grepl("- CT", contrast) ~ "CT",
    grepl("- Olaparib", contrast) ~ "Olaparib"
  )) %>%
  filter(!is.na(Drug)) %>%
  mutate(
    Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")),
    Drug = factor(Drug, levels = c("CT", "Olaparib"))
  )

group_colors <- c("CT" = "#C43826", "Olaparib" = "#858FB0")

pdf("M2_P1040.pdf", width = 12, height = 9)
ggplot(plot_data, aes(x = Drug, y = estimate, fill = Drug)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  facet_grid(Drug ~ Previous_Treatment, scales = "fixed", space = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(-10, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Previous Treatment", y = "Estimated Difference in AUC",
       title = "Estimated Difference in Tumour Size between Previous\nTreatment and Current Treatment") +
  theme_bw() +
  theme(
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none"
  )
dev.off()

##### Supplement: All contrast estimations for CT vs Olaparib ####
plot_data <- emm_pairs_summary %>%
  mutate(
    Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")),
    Contrast = factor(contrast, levels = c("Untreated - CT", "Untreated - Olaparib", "CT - Olaparib"))
  )

group_colors <- c("Untreated - CT" = "#C43826", "Untreated - Olaparib" = "#858FB0", "CT - Olaparib" = "#5D8AA8")

pdf("M2_P1040_suppl.pdf", width = 12, height = 9)
ggplot(plot_data, aes(x = Contrast, y = estimate, fill = Contrast)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  facet_grid(Contrast ~ Previous_Treatment, scales = "fixed", space = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(-10, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Previous Treatment", y = "Estimated Difference in log(Size)",
       title = "Estimated Difference in Tumour Size between Previous\nTreatment and Current Treatment") +
  theme_bw() +
  theme(
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none"
  )
dev.off()

# Predicted volumes between Untreated and other Cohorts
combined_data_list <- list()

for (prev_treatment in names(model_list)) {
  m_t <- model_list[[prev_treatment]]
  
  tumrs <- subset_models[[prev_treatment]]$data
  
  dc <- distinct(select(tumrs, Mouse, Cohort))
  pframe <- expand_grid(Mouse = unique(tumrs$Mouse), Time = 84) %>%
    full_join(dc, by = "Mouse") %>%
    mutate(lSize = predict(m_t, newdata = .))
  
  comb <- bind_rows(list(data = tumrs, model = pframe), .id = "Type")
  comb$ModelSize <- exp(comb$lSize)
  comb$Previous_Treatment <- prev_treatment
  
  combined_data_list[[prev_treatment]] <- comb
}

final_combined_data <- bind_rows(combined_data_list)

predicted_final_bxp <- final_combined_data %>%
  filter(Time == 84) %>%
  distinct(Previous_Treatment, Mouse, Cohort, .keep_all = TRUE) %>%
  rename(Current_Treatment = Cohort, Size = ModelSize)

# Save the final predictions
write.csv(predicted_final_bxp, file = "P1040_revisions_M4.csv", row.names = FALSE)

# Visualization of predicted final volumes
predicted_final_bxp <- predicted_final_bxp %>%
  mutate(Current_Treatment = factor(Current_Treatment, levels = c("Untreated", "CT", "Olaparib")),
         Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")))

group_colors <- c("Untreated" = "black", "CT" = "#C43826", "Olaparib" = "#858FB0")

pdf("M4_P1040.pdf", width = 12, height = 9)
ggplot(predicted_final_bxp, aes(x = Current_Treatment, y = lSize)) +
  geom_boxplot(aes(fill = Current_Treatment), width = 0.7, show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(aes(color = Current_Treatment), width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_grid(. ~ Previous_Treatment, scales = "free_x", space = "free_x") +
  labs(x = "Current Treatment", y = "Predicted log(volume) at day 84",
       title = "Log Tumour Volume Predictions\nby Previous and Current Treatments") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    panel.spacing = unit(0.5, "lines")
  ) +
  coord_cartesian(ylim = c(min(predicted_final_bxp$lSize) * 0.9, max(predicted_final_bxp$lSize) * 1.1)) +
  scale_y_continuous(breaks = seq(floor(min(predicted_final_bxp$lSize)), ceiling(max(predicted_final_bxp$lSize)), by = 1))
dev.off()


####Adding the Untreated mice predictions ####

# Adding the Untreated mice predictions
model_list <- list(
  Untreated = lmer_Untreated,
  CT = lmer_CT,
  Olaparib = lmer_Olaparib
)

combined_data_list <- list()

for (prev_treatment in names(model_list)) {
  m_t <- model_list[[prev_treatment]]
  
  # Select the appropriate dataset
  tumrs <- switch(prev_treatment,
                  Untreated = Untreated_PreviousTreatment,
                  CT = CT_Previous_Treatment,
                  Olaparib = Olaparib_Previous_Treatment)
  
  # Create a prediction frame for all cohorts, including Untreated
  dc <- distinct(select(tumrs, Mouse))
  pframe <- expand_grid(Mouse = unique(tumrs$Mouse), Time = 84, Cohort = c("Untreated", "CT", "Olaparib")) %>%
    full_join(dc, by = "Mouse") %>%
    mutate(lSize = predict(m_t, newdata = .))
  
  # Combine actual data with predicted data
  comb <- bind_rows(list(data = tumrs, model = pframe), .id = "Type")
  comb$ModelSize <- exp(comb$lSize)
  comb$Previous_Treatment <- prev_treatment
  
  combined_data_list[[prev_treatment]] <- comb
}

# Combine all data into one dataframe
final_combined_data <- bind_rows(combined_data_list)

# Prepare data for plotting
predicted_final_bxp <- final_combined_data %>%
  filter(Time == 84) %>%
  distinct(Previous_Treatment, Mouse, Cohort, .keep_all = TRUE) %>%
  rename(Current_Treatment = Cohort, Size = ModelSize)

# Save the predicted data to a CSV file
write.csv(predicted_final_bxp, file = "P1040_revisions_M4.csv", row.names = FALSE)

# Prepare the data for plotting with appropriate factor levels
predicted_final_bxp <- predicted_final_bxp %>%
  mutate(Current_Treatment = factor(Current_Treatment, levels = c("Untreated", "CT", "Olaparib")),
         Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")))

# Define color palette
group_colors <- c("Untreated" = "black", "CT" = "#C43826", "Olaparib" = "#858FB0")

# Create the plot and save it as a PDF
pdf("M4_P1040.pdf", width = 12, height = 9)

ggplot(predicted_final_bxp, aes(x = Current_Treatment, y = lSize)) +
  geom_boxplot(aes(fill = Current_Treatment), width = 0.7, show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(aes(color = Current_Treatment), width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_grid(. ~ Previous_Treatment, scales = "free_x", space = "free_x") +
  labs(x = "Current Treatment", 
       y = "Predicted log(volume) at day 84", 
       title = "Log Tumour Volume Predictions\nby Previous and Current Treatments") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 14),
    panel.spacing = unit(0.5, "lines")
  ) +
  coord_cartesian(ylim = c(min(predicted_final_bxp$lSize) * 0.9, max(predicted_final_bxp$lSize) * 1.1)) +
  scale_y_continuous(breaks = seq(floor(min(predicted_final_bxp$lSize)), ceiling(max(predicted_final_bxp$lSize)), by = 1))

dev.off()

P1040 <- read.csv("/Users/guerre01/Downloads/Sequential drug treatment in vivo trial/P1040.csv", 
                  header = TRUE, stringsAsFactors = FALSE)

P1040 <- dplyr:::select(P1040, Model:Day, Volume, Perc_vol_change)
P1040 <- rename(P1040, Time = Day, Size = Volume)
P1040 <- filter(P1040, Group != "")


P1040$Size <- as.numeric(P1040$Size)
P1040$Group <- factor(P1040$Group)
P1040 <- P1040 %>% drop_na(Size) %>% as.data.frame()
P1040 <- P1040 %>% filter(Size >0)
P1040$lSize <- log(P1040$Size)
P1040 <- P1040[is.finite(P1040$lSize), ]
P1040$Groups <- paste(P1040$Group, P1040$Cohort, sep = " -> ")

P1040 <- P1040 %>%
  mutate(Cohort = factor(Cohort, levels = c("Untreated", "CT", "Olaparib")))

ggplot(P1040, aes(x = Time, y = Size, group = Mouse, color = Cohort)) +
  geom_line() +
  facet_grid(Cohort ~ Previous_Treatment) +
  labs(x = "Days since start of treatment",
       y = "Tumour volume (mm³)",
       title = "Previous Treatment Arm Group") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_color_manual(values = c("Untreated" = "black", 
                                "CT" = "#C43826", 
                                "Olaparib" = "#858FB0"))

#  subset for Untreated Previous Treatment
Untreated_PreviousTreatment <- Untreated_PreviousTreatment %>%
  mutate(Cohort = factor(Cohort, levels = c("Untreated", "CT", "Olaparib")))

Untreated_PreviousTreatment <- P1040 %>%
  filter(Previous_Treatment == "Untreated")

#  subset for CT Previous Treatment
CT_Previous_Treatment <- CT_Previous_Treatment %>%
  mutate(Cohort = factor(Cohort, levels = c("Untreated", "CT", "Olaparib")))

CT_Previous_Treatment <- P1040 %>%
  filter(Previous_Treatment == "CT")

#  subset for Olaparib Previous Treatment
Olaparib_Previous_Treatment <- Olaparib_Previous_Treatment %>%
  mutate(Cohort = factor(Cohort, levels = c("Untreated", "CT", "Olaparib")))

Olaparib_Previous_Treatment <- P1040 %>%
  filter(Previous_Treatment == "Olaparib")




create_plot <- function(data, title) {
  ggplot(data, aes(x = Time, y = Size, group = Mouse, color = Cohort)) +
    geom_line() +
    facet_wrap(~ Cohort, ncol = 1, strip.position = "right") +
    labs(x = "Days since start of treatment",
         y = "Tumour volume (mm³)",
         title = title) +
    theme_minimal() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", angle = -90, hjust = 0),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none") +
    scale_color_manual(values = c("Untreated" = "black", 
                                  "CT" = "#C43826", 
                                  "Olaparib" = "#858FB0"))
}

plot_Untreated <- create_plot(Untreated_PreviousTreatment, "Untreated Previous Treatment")
plot_CT <- create_plot(CT_Previous_Treatment, "CT Previous Treatment")
plot_Olaparib <- create_plot(Olaparib_Previous_Treatment, "Olaparib Previous Treatment")


print(plot_Untreated)
print(plot_CT)
print(plot_Olaparib)


levels(Untreated_PreviousTreatment$Cohort)
levels(CT_Previous_Treatment$Cohort)
levels(Olaparib_Previous_Treatment$Cohort)


# For Untreated_PreviousTreatment subset
lmer_Untreated <- lmerTest::lmer(lSize ~ Time * Cohort + (Time | Mouse),
                                 data = Untreated_PreviousTreatment,
                                 control = lmerControl(optimizer = "bobyqa"))

# For CT_Previous_Treatment subset
lmer_CT <- lmerTest::lmer(lSize ~ Time * Cohort + (Time | Mouse),
                          data = CT_Previous_Treatment,
                          control = lmerControl(optimizer = "bobyqa"))

# For Olaparib_Previous_Treatment subset
lmer_Olaparib <- lmerTest::lmer(lSize ~ Time * Cohort + (Time | Mouse),
                                data = Olaparib_Previous_Treatment,
                                control = lmerControl(optimizer = "bobyqa"))

#### Plot unique fittting per unique Model per Condition ####
library(ggplot2)
library(dplyr)

create_pred_data <- function(model, data) {
  new_data <- expand.grid(
    Time = seq(min(data$Time), max(data$Time), length.out = 100),
    Cohort = unique(data$Cohort)
  )
  new_data$Predicted <- predict(model, newdata = new_data, re.form = NA)
  return(new_data)
}

# Create prediction data for each model
pred_untreated <- create_pred_data(lmer_Untreated, Untreated_PreviousTreatment)
pred_ct <- create_pred_data(lmer_CT, CT_Previous_Treatment)
pred_olaparib <- create_pred_data(lmer_Olaparib, Olaparib_Previous_Treatment)

# Add Previous_Treatment information
pred_untreated$Previous_Treatment <- "Untreated"
pred_ct$Previous_Treatment <- "CT"
pred_olaparib$Previous_Treatment <- "Olaparib"

all_pred_data <- rbind(pred_untreated, pred_ct, pred_olaparib)

curve_colors <- c("CT" = "#C43826",
                  "Olaparib" = "#858FB0",
                  "Untreated" = "black")

ggplot(all_pred_data, aes(x = Time, y = exp(Predicted), color = Cohort)) +
  geom_line() +
  facet_wrap(~ Previous_Treatment, scales = "free_y", ncol = 3) +
  labs(x = "Days", y = "Predicted tumour size", 
       title = "Tumour growth by previous and current treatments") +
  theme_bw() +
  scale_y_log10() +
  scale_color_manual(values = curve_colors) +
  theme(
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

##Fitting mouse
library(ggplot2)
library(dplyr)
library(purrr)

# Function to create prediction data for each mouse, including random effects
create_mouse_pred_data <- function(model, data) {
  mice <- unique(data$Mouse)
  
  pred_data <- map_dfr(mice, function(mouse) {
    new_data <- data.frame(
      Time = seq(min(data$Time), max(data$Time), length.out = 100),
      Mouse = mouse,
      Cohort = unique(data$Cohort[data$Mouse == mouse])
    )
    new_data$Predicted <- predict(model, newdata = new_data, re.form = NULL)
    new_data$Previous_Treatment <- unique(data$Previous_Treatment)
    return(new_data)
  })
  
  return(pred_data)
}

# Create prediction data for each model
pred_untreated <- create_mouse_pred_data(lmer_Untreated, Untreated_PreviousTreatment)
pred_ct <- create_mouse_pred_data(lmer_CT, CT_Previous_Treatment)
pred_olaparib <- create_mouse_pred_data(lmer_Olaparib, Olaparib_Previous_Treatment)

# Combine all prediction data
all_pred_data <- rbind(pred_untreated, pred_ct, pred_olaparib)

# Combine actual data
all_actual_data <- rbind(Untreated_PreviousTreatment, CT_Previous_Treatment, Olaparib_Previous_Treatment)

# Define color palette
curve_colors <- c("CT" = "#C43826", "Olaparib" = "#858FB0", "Untreated" = "black")

# Remove predicted values greater than 1700
all_pred_data <- all_pred_data %>%
  group_by(Mouse, Cohort, Previous_Treatment) %>%
  mutate(Predicted_exp = exp(Predicted)) %>%
  filter(Predicted_exp <= 1700) %>%
  ungroup()

ggplot() +
  geom_point(data = all_actual_data, aes(x = Time, y = Size, group = Mouse, color = Cohort),alpha = 0.7) +
  geom_line(data = all_pred_data, aes(x = Time, y = Predicted_exp, group = Mouse, color = Cohort), linetype = "dashed", alpha = 0.5) +
  facet_grid(Previous_Treatment ~ Cohort, scales = "free_y") +
  labs(x = "Days since Treatment", y = "Tumor Volume", 
       title = "Individual mouse tumour growth by previous and current treatments") +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 1700),
                     labels = scales::comma,
                     limits = c(0, 1700)) +
  scale_color_manual(values = curve_colors) +
  coord_cartesian(xlim = c(0, 90)) +
  scale_x_continuous(breaks = seq(0, 90, by = 20)) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")  # Increase space between facets
  )

#Interaction between Untreated post-treatment and other Groups #####
library(dplyr)
library(lmerTest)

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
  tumrs <- subset_models[[prev_treatment]]$data
  m_t <- subset_models[[prev_treatment]]$model
  
  summary_m_t <- summary(m_t)
  fixed_eff_df <- as.data.frame(summary_m_t$coefficients)
  fixed_eff_df$Parameter <- rownames(summary_m_t$coefficients)
  fixed_eff_df$Previous_Treatment <- prev_treatment
  
  ci <- confint(m_t,level = 0.95)
  
  fixed_eff_df$CI_low <- NA
  fixed_eff_df$CI_high <- NA
  
  for (param in fixed_eff_df$Parameter) {
    if (param %in% rownames(ci)) {
      ci_param <- ci[param, ]
      fixed_eff_df$CI_low[fixed_eff_df$Parameter == param] <- ci_param[1]
      fixed_eff_df$CI_high[fixed_eff_df$Parameter == param] <- ci_param[2]
    }
  }
  
  fixed_effects_summary <- rbind(
    fixed_effects_summary, 
    fixed_eff_df[, c("Previous_Treatment", "Parameter", "Estimate", "Std. Error", "Pr(>|t|)", "CI_low", "CI_high")]
  )
}

fixed_effects_summary <- fixed_effects_summary %>%
  mutate(Model = case_when(
    Previous_Treatment == "Untreated" ~ "lmer_Untreated",
    Previous_Treatment == "CT" ~ "lmer_CT",
    Previous_Treatment == "Olaparib" ~ "lmer_Olaparib",
    TRUE ~ NA_character_  # This handles any unexpected values
  ))

# Reorder columns to put the new 'Model' column right after 'Previous_Treatment'
fixed_effects_summary <- fixed_effects_summary %>%
  dplyr:::select(Previous_Treatment, Model, everything())

names(fixed_effects_summary) <- c("Previous_Treatment", "Model", "Parameter", "Estimate", "StdError", "Pr", "CI_low", "CI_high")
fixed_effects_summary$Pr <- floor(fixed_effects_summary$Pr * 100) / 100

write.csv2(fixed_effects_summary, file = "P1040_revisions_M1.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)

interaction_data <- fixed_effects_summary %>%
  filter(grepl("Time:Cohort", Parameter)) %>%
  separate(Parameter, into = c("Time", "Cohort"), sep = ":") %>%
  mutate(Cohort = gsub("Cohort", "", Cohort),
         Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")))

interaction_data <- interaction_data %>%
  filter(!is.na(Estimate))

#option A
ggplot(interaction_data, aes(x = Cohort, y = Estimate, color = Cohort)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                width = 0.2, position = position_dodge(0.5)) +
  facet_wrap(~ Previous_Treatment, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("CT" = "#C43826", 
                                "Olaparib" = "#858FB0")) +
  labs(x = "Current Treatment", 
       y = "Change in growth rate compared to Untreated") +
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

#option B
ggplot(interaction_data, aes(x = Cohort, y = Estimate, color = Cohort)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                width = 0.2, position = position_dodge(0.5)) +
  facet_grid(Cohort ~ Previous_Treatment, scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("CT" = "#C43826", 
                                "Olaparib" = "#858FB0")) +
  labs(x = "Current Treatment", 
       y = "Change in growth rate compared to Untreated") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",  # This removes the legend
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


#AUC between Untreated post-treatment and other Groups####
library(emmeans)
library(dplyr)

# Initialize list for storing models
model_list <- list(
  Untreated = lmer_Untreated,
  CT = lmer_CT,
  Olaparib = lmer_Olaparib
)

# Initialize data frame for storing emmeans pairs summary
emm_pairs_summary <- data.frame()

# Iterate over each unique model
for (prev_treatment in names(model_list)) {
  m_t <- model_list[[prev_treatment]]
  
  # Compute estimated marginal means for 'Cohort' at Time c(0, 77)
  emm <- emmeans(m_t, specs = ~ Cohort, at = list(Time = c(0, 77)))
  
  # Generate pairwise comparisons
  emm_pairs <- pairs(emm, scale = 77)
  
  # Convert summary of emm_pairs to a data frame and append model identifier
  emm_pairs_df <- summary(emm_pairs, infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    mutate(Previous_Treatment = prev_treatment)
  
  # Append to the summary data frame
  emm_pairs_summary <- rbind(emm_pairs_summary, emm_pairs_df)
}

# Add Model column based on Previous_Treatment
emm_pairs_summary <- emm_pairs_summary %>%
  mutate(Model = case_when(
    Previous_Treatment == "Untreated" ~ "lmer_Untreated",
    Previous_Treatment == "CT" ~ "lmer_CT",
    Previous_Treatment == "Olaparib" ~ "lmer_Olaparib",
    TRUE ~ NA_character_
  )) %>%
  dplyr:::select(Previous_Treatment, Model, everything())

write.csv2(emm_pairs_summary, file = "P1040_revisions_M2.csv", row.names = FALSE)


library(dplyr)
library(tidyr)
library(ggplot2)

plot_data <- emm_pairs_summary %>%
  filter(grepl("Untreated -", contrast)) %>%
  mutate(Drug = case_when(
    grepl("- CT", contrast) ~ "CT",
    grepl("- Olaparib", contrast) ~ "Olaparib",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Drug))

plot_data <- plot_data %>%
  mutate(
    Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")),
    Drug = factor(Drug, levels = c("CT", "Olaparib"))
  )

group_colors <- c("CT" = "#C43826", "Olaparib" = "#858FB0")

ggplot(plot_data, aes(x = Drug, y = estimate, fill = Drug)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  facet_grid(Drug ~ Previous_Treatment, scales = "fixed", space = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(-10, 100), 
                     breaks = seq(0, 100, by = 20)) +
  labs(x = "Previous Treatment",  # Remove x-axis label
       y = "Estimated Difference in log(Size)", 
       fill = "Drug") +
  theme_bw() +
  theme(
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )
  )

##### Supplement All contrast estimations CT vs Olaparib ####
library(dplyr)
library(ggplot2)

plot_data <- plot_data %>%
  mutate(
    Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")),
    Contrast = factor(contrast, levels = c("Untreated - CT", "Untreated - Olaparib", "CT - Olaparib"))
  )


group_colors <- c("Untreated - CT" = "#C43826", "Untreated - Olaparib" = "#858FB0", "CT - Olaparib" = "#5D8AA8")


ggplot(plot_data, aes(x = Contrast, y = estimate, fill = Contrast)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), color = "black") +
  facet_grid(Contrast ~ Previous_Treatment, scales = "fixed", space = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(-10, 100), 
                     breaks = seq(0, 100, by = 20)) +
  labs(x = "Previous Treatment",  # Remove x-axis label
       y = "Estimated Difference in log(Size)", 
       fill = "Contrast") +
  theme_bw() +
  theme(
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) 

#Time plus Time*GroupCT between Untreated post-treatment and other Groups####
K <- c("(Intercept)" = 0, 
       "Time" = 1,
       "CohortCT" = 0, 
       "CohortOlaparib" = 0, 
       "Time:CohortCT" = 0, 
       "Time:CohortOlaparib" = 1)

Va <- vcov(lmer_CT)
ddf <- get_Lb_ddf(lmer_CT, K)

#from this we can construct the usual t statistc and corresponding p value
b.hat <- fixef(lmer_CT)
Lb.hat <- sum(K * b.hat)
#confint
alpha = 0.05
SE <- as.numeric(sqrt(K %*% Va %*% K))
Lb.hat + SE * qt(c(0.5 * alpha, 1 - 0.5* alpha), ddf)
Lb.hat

library(lme4)
library(lmerTest)
library(dplyr)

model_list <- list(
  Untreated = lmer_Untreated,
  CT = lmer_CT,
  Olaparib = lmer_Olaparib
)

results_summary <- data.frame(
  Previous_Treatment = character(),
  Current_Treatment = character(),
  Lb_hat = numeric(),
  SE = numeric(),
  CI_low = numeric(),
  CI_high = numeric(),
  Pr = numeric(),
  stringsAsFactors = FALSE
)

alpha = 0.05

for (prev_treatment in names(model_list)) {
  m_t <- model_list[[prev_treatment]]
  
  current_treatments <- levels(m_t@frame$Cohort)[-1]  # Exclude the reference level (Untreated)
  
  for (current_treatment in current_treatments) {
    interaction_term <- paste0("Time:Cohort", current_treatment)
    
    if (interaction_term %in% names(fixef(m_t))) {
      K <- setNames(rep(0, length(fixef(m_t))), names(fixef(m_t)))
      K["Time"] <- 1
      K[interaction_term] <- 1
      
      Va <- vcov(m_t)
      b.hat <- fixef(m_t)
      Lb.hat <- sum(K * b.hat)
      SE <- sqrt(sum((K %*% Va) * K))
      ddf <- df.residual(m_t)
      
      confint <- Lb.hat + c(-1, 1) * SE * qt(1 - 0.5 * alpha, ddf)
      
      # Calculate p-value
      t_stat <- Lb.hat / SE
      p_value <- 2 * (1 - pt(abs(t_stat), ddf))
      
      results_summary <- rbind(results_summary, 
                               data.frame(Previous_Treatment = prev_treatment,
                                          Current_Treatment = current_treatment,
                                          Lb_hat = Lb.hat,
                                          SE = SE,
                                          CI_low = confint[1],
                                          CI_high = confint[2],
                                          Pr = p_value))
    }
  }
}

results_summary <- results_summary %>%
  mutate(Model = case_when(
    Previous_Treatment == "Untreated" ~ "lmer_Untreated",
    Previous_Treatment == "CT" ~ "lmer_CT",
    Previous_Treatment == "Olaparib" ~ "lmer_Olaparib",
    TRUE ~ NA_character_
  )) %>%
  dplyr:::select(Previous_Treatment, Model, Current_Treatment, Lb_hat, SE, Pr, CI_low, CI_high)

# Rename columns
names(results_summary) <- c("Previous_Treatment", "Model", "Current_Treatment", "Estimate", "StdError", "Pr", "CI_low", "CI_high")

# Round p-values
results_summary$Pr <- floor(results_summary$Pr * 100) / 100

# Write the results to a CSV file
write.csv2(results_summary, file = "P1040_revisions_M3.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)

group_colors <- c("CT" = "#C43826", "Olaparib" = "#858FB0")

ggplot( , aes(x = Previous_Treatment, y = Estimate, color = Current_Treatment)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                position = position_dodge(width = 0.5), width = 0.2) +
  facet_grid(Current_Treatment ~ Previous_Treatment, scales = "free_y", space = "free_y") +
  scale_color_manual(values = group_colors, name = "Current Treatment") +
  labs(x = "Previous Treatment", 
       y = "Change in growth rate compared to Untreated", 
       title = "Treatment Effects on Tumor Growth Rate") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing = unit(1, "lines")
  )


###Adding the Time =1 only for the previous code ####

model_list <- list(
  Untreated = lmer_Untreated,
  CT = lmer_CT,
  Olaparib = lmer_Olaparib
)

results_summary <- data.frame(
  Previous_Treatment = character(),
  Current_Treatment = character(),
  Parameter = character(),
  Lb_hat = numeric(),
  SE = numeric(),
  CI_low = numeric(),
  CI_high = numeric(),
  Pr = numeric(),
  stringsAsFactors = FALSE
)

alpha = 0.05

for (prev_treatment in names(model_list)) {
  m_t <- model_list[[prev_treatment]]
  
  # Estimate for Time parameter alone
  K_time <- setNames(rep(0, length(fixef(m_t))), names(fixef(m_t)))
  K_time["Time"] <- 1
  
  Va <- vcov(m_t)
  b.hat <- fixef(m_t)
  Lb.hat_time <- sum(K_time * b.hat)
  SE_time <- sqrt(sum((K_time %*% Va) * K_time))
  ddf_time <- df.residual(m_t)
  
  confint_time <- Lb.hat_time + c(-1, 1) * SE_time * qt(1 - 0.5 * alpha, ddf_time)
  
  # Calculate p-value for Time
  t_stat_time <- Lb.hat_time / SE_time
  p_value_time <- 2 * (1 - pt(abs(t_stat_time), ddf_time))
  
  results_summary <- rbind(results_summary, 
                           data.frame(Previous_Treatment = prev_treatment,
                                      Current_Treatment = "Untreated",
                                      Parameter = "Time",
                                      Lb_hat = Lb.hat_time,
                                      SE = SE_time,
                                      CI_low = confint_time[1],
                                      CI_high = confint_time[2],
                                      Pr = p_value_time))
  
  current_treatments <- levels(m_t@frame$Cohort)[-1]  # Exclude the reference level (Untreated)
  
  for (current_treatment in current_treatments) {
    interaction_term <- paste0("Time:Cohort", current_treatment)
    
    if (interaction_term %in% names(fixef(m_t))) {
      K <- setNames(rep(0, length(fixef(m_t))), names(fixef(m_t)))
      K["Time"] <- 1
      K[interaction_term] <- 1
      
      Va <- vcov(m_t)
      b.hat <- fixef(m_t)
      Lb.hat <- sum(K * b.hat)
      SE <- sqrt(sum((K %*% Va) * K))
      ddf <- df.residual(m_t)
      
      confint <- Lb.hat + c(-1, 1) * SE * qt(1 - 0.5 * alpha, ddf)
      
      # Calculate p-value
      t_stat <- Lb.hat / SE
      p_value <- 2 * (1 - pt(abs(t_stat), ddf))
      
      results_summary <- rbind(results_summary, 
                               data.frame(Previous_Treatment = prev_treatment,
                                          Current_Treatment = current_treatment,
                                          Parameter = paste("Time +", interaction_term),
                                          Lb_hat = Lb.hat,
                                          SE = SE,
                                          CI_low = confint[1],
                                          CI_high = confint[2],
                                          Pr = p_value))
    }
  }
}

results_summary <- results_summary %>%
  mutate(Model = case_when(
    Previous_Treatment == "Untreated" ~ "lmer_Untreated",
    Previous_Treatment == "CT" ~ "lmer_CT",
    Previous_Treatment == "Olaparib" ~ "lmer_Olaparib",
    TRUE ~ NA_character_
  )) %>%
  dplyr:::select(Previous_Treatment, Model, Current_Treatment, Parameter, Lb_hat, SE, Pr, CI_low, CI_high)

names(results_summary) <- c("Previous_Treatment", "Model", "Current_Treatment", "Parameter", "Estimate", "StdError", "Pr", "CI_low", "CI_high")

results_summary$Pr <- floor(results_summary$Pr * 100) / 100

write.csv2(results_summary, file = "P1040_revisions_M3.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)

plot_data <- results_summary %>%
  mutate(
    Parameter_Type = ifelse(Parameter == "Time", "Untreated", "Interaction"),
    Current_Treatment = factor(ifelse(Current_Treatment == "Untreated", "Untreated", Current_Treatment),
                               levels = c("Untreated", "CT", "Olaparib"))
  )


group_colors <- c("CT" = "#C43826", "Olaparib" = "#858FB0", "Untreated" = "black")

ggplot(plot_data, aes(x = Previous_Treatment, y = Estimate, color = Current_Treatment)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                position = position_dodge(width = 0.5), width = 0.2) +
  facet_grid(Current_Treatment ~ Previous_Treatment, scales = "free_y", space = "free_y") +
  scale_color_manual(values = group_colors, name = "Current Treatment") +
  labs(x = "Previous Treatment", 
       y = "Change in growth rate") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing = unit(1, "lines")
  ) 


plot_data <- plot_data %>%
  mutate(
    Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")),
    Current_Treatment = factor(Current_Treatment, levels = c("Untreated", "CT", "Olaparib"))
  )

group_colors <- c("CT" = "#C43826", "Olaparib" = "#858FB0", "Untreated" = "black")

ggplot(plot_data, aes(x = Previous_Treatment, y = Estimate, color = Current_Treatment)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), 
                position = position_dodge(width = 0.5), width = 0.2) +
  facet_grid(Current_Treatment ~ Previous_Treatment, scales = "free_y", space = "free_y") +
  scale_color_manual(values = group_colors, name = "Current Treatment") +
  labs(x = "Previous Treatment", 
       y = "Change in growth rate") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  theme_bw() +
  coord_flip() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    strip.background.x = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.background.y = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.spacing = unit(1, "lines")
  )

# Predicted volumes between Untreated and other Cohorts ####

library(dplyr)
library(tidyr)
library(lmerTest)

model_list <- list(
  Untreated = lmer_Untreated,
  CT = lmer_CT,
  Olaparib = lmer_Olaparib
)

combined_data_list <- list()

for (prev_treatment in names(model_list)) {
  m_t <- model_list[[prev_treatment]]
  
  # Select the appropriate dataset
  if (prev_treatment == "Untreated") {
    tumrs <- Untreated_PreviousTreatment
  } else if (prev_treatment == "CT") {
    tumrs <- CT_Previous_Treatment
  } else if (prev_treatment == "Olaparib") {
    tumrs <- Olaparib_Previous_Treatment
  }
  
  # Create a prediction frame
  dc <- distinct(dplyr:::select(tumrs, Mouse, Cohort))
  pframe <- expand_grid(Mouse = unique(tumrs$Mouse), Time = 84) %>%
    full_join(dc, by = "Mouse") %>%
    mutate(lSize = predict(m_t, newdata = .))
  
  # Combine actual data with predicted data
  comb <- bind_rows(list(data = tumrs, model = pframe), .id = "Type")
  comb$ModelSize <- exp(comb$lSize)
  comb$Previous_Treatment <- prev_treatment
  
  combined_data_list[[prev_treatment]] <- comb
}

final_combined_data <- bind_rows(combined_data_list)

options(scipen = 999)

predicted_final_bxp <- final_combined_data %>%
  filter(Time == 77 & Cohort != "Untreated") %>%
  dplyr:::select(Previous_Treatment, Mouse, Cohort, lSize, ModelSize) %>%
  distinct(Previous_Treatment, Mouse, .keep_all = TRUE)

predicted_final_bxp <- predicted_final_bxp %>%
  rename(Drugs = Cohort, Size = ModelSize)

write.csv(predicted_final_bxp, file = "P1040_revisions_M4.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)

group_colors <- c("CT" = "#C43826", "Olaparib" = "#858FB0")

ggplot(predicted_final_bxp, aes(x = Previous_Treatment, y = lSize)) +
  geom_boxplot(aes(fill = Drugs), position = position_dodge(0.8), width = 0.7, show.legend = TRUE) +
  geom_jitter(aes(color = Drugs), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), alpha = 0.5) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = "Previous Treatment", 
       y = "Predicted log(volume) at the end of the treatment", 
       fill = "Current Treatment",
       color = "Current Treatment",
       title = "Log Tumour Volume Predictions by Previous and Current Treatments") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(min(predicted_final_bxp$lSize) * 0.9, max(predicted_final_bxp$lSize) * 1.1)) +
  scale_y_continuous(breaks = seq(floor(min(predicted_final_bxp$lSize)), ceiling(max(predicted_final_bxp$lSize)), by = 1))


####Adding the Untreated mice predictions ####

model_list <- list(
  Untreated = lmer_Untreated,
  CT = lmer_CT,
  Olaparib = lmer_Olaparib
)

combined_data_list <- list()

for (prev_treatment in names(model_list)) {
  m_t <- model_list[[prev_treatment]]
  
  # Select the appropriate dataset
  if (prev_treatment == "Untreated") {
    tumrs <- Untreated_PreviousTreatment
  } else if (prev_treatment == "CT") {
    tumrs <- CT_Previous_Treatment
  } else if (prev_treatment == "Olaparib") {
    tumrs <- Olaparib_Previous_Treatment
  }
  
  # Create a prediction frame for all Cohorts, including Untreated
  dc <- distinct(dplyr:::select(tumrs, Mouse))
  pframe <- expand_grid(Mouse = unique(tumrs$Mouse), Time = 84, Cohort = c("Untreated", "CT", "Olaparib")) %>%
    full_join(dc, by = "Mouse") %>%
    mutate(lSize = predict(m_t, newdata = .))
  
  # Combine actual data with predicted data
  comb <- bind_rows(list(data = tumrs, model = pframe), .id = "Type")
  comb$ModelSize <- exp(comb$lSize)
  comb$Previous_Treatment <- prev_treatment
  
  combined_data_list[[prev_treatment]] <- comb
}

final_combined_data <- bind_rows(combined_data_list)

options(scipen = 999)

predicted_final_bxp <- final_combined_data %>%
  filter(Time == 84) %>%
  dplyr:::select(Previous_Treatment, Mouse, Cohort, lSize, ModelSize) %>%
  distinct(Previous_Treatment, Mouse, Cohort, .keep_all = TRUE)

predicted_final_bxp <- predicted_final_bxp %>%
  rename(Current_Treatment = Cohort, Size = ModelSize)

write.csv(predicted_final_bxp, file = "P1040_revisions_M4.csv", row.names = FALSE)

predicted_final_bxp <- predicted_final_bxp %>%
  mutate(Current_Treatment = factor(Current_Treatment, levels = c("Untreated", "CT", "Olaparib")),
         Previous_Treatment = factor(Previous_Treatment, levels = c("Untreated", "CT", "Olaparib")))

group_colors <- c("Untreated" = "black", "CT" = "#C43826", "Olaparib" = "#858FB0")

ggplot(predicted_final_bxp, aes(x = Current_Treatment, y = lSize)) +
  geom_boxplot(aes(fill = Current_Treatment), width = 0.7, show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(aes(color = Current_Treatment), width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_grid(. ~ Previous_Treatment, scales = "free_x", space = "free_x") +
  labs(x = "Current Treatment", 
       y = "Predicted log(volume) at day 84", 
       title = "Log Tumour Volume Predictions by Previous and Current Treatments") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#FFE4CB", colour = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  coord_cartesian(ylim = c(min(predicted_final_bxp$lSize) * 0.9, max(predicted_final_bxp$lSize) * 1.1)) +
  scale_y_continuous(breaks = seq(floor(min(predicted_final_bxp$lSize)), ceiling(max(predicted_final_bxp$lSize)), by = 1))
t
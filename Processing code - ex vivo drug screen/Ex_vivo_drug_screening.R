# Script for processing ex vivo high-throughput drug screening data
# Created for publication alongside a manuscript submitted to Nature

# Loading required libraries
library(ggplot2)
library(flux)
library(dplyr)

# Function to perform bootstrapping and generate fitted dose-response curves
getBootSample <- function(data, B=100) {
  # Initialize variables
  unique_doses <- length(unique(data$Concentration))
  bootstrap_results <- matrix(NA, unique_doses, B)
  num_repeats <- nrow(data) / unique_doses
  bootstrap_indices <- matrix(sample(size=num_repeats * unique_doses * B, x=1:num_repeats, replace=TRUE), nrow=B)
  
  # Perform bootstrapping
  for (b in 1:B) {
    bootstrap_data <- NULL
    for (k in 1:num_repeats) {
      bootstrap_data <- c(bootstrap_data, bootstrap_indices[b, ((k-1) * unique_doses + 1):((k-1) * unique_doses + unique_doses)] + num_repeats * (0:(unique_doses-1)))
    }
    bootstrap_data <- data[bootstrap_data, ]
    
    # Fit isotonic regression
    iso_fit <- isoreg(x=bootstrap_data$Concentration, y=bootstrap_data$Response)
    
    # Store the fitted values for each dose
    bootstrap_results[,b] <- iso_fit$yf[seq(1, length.out=unique_doses, by=nrow(bootstrap_data)/unique_doses)]
  }
  
  return(bootstrap_results)
}

# Load input data
input_data <- read.csv("/path/to/input data.csv", stringsAsFactors=FALSE)

# Data preprocessing
# Ensure that concentration is numeric and the data is ordered
input_data <- input_data %>%
  mutate(Concentration = as.numeric(Concentration)) %>%
  arrange(Model, Drug, Previous_treatment, Concentration)

# Prepare an empty dataframe to store results
results <- input_data %>%
  select(Model, Previous_treatment, Drug) %>%
  distinct() %>%
  mutate(AUC = NA, AUC_LI = NA, AUC_UI = NA, AUC_var = NA, IC50 = NA, IC50_LI = NA, IC50_UI = NA)

# List to store fitted models and bootstrapped confidence intervals
model_fits <- list()
lower_intervals <- list()
upper_intervals <- list()

# Predefined x-values for concentrations
x_values <- c(0.01, 0.03, 0.10, 0.30, 1.00, 3.00, 10.00)

# Loop through each unique combination of Model, Previous_treatment, and Drug
for (i in 1:nrow(results)) {
  # Subset data for the current combination
  subset_data <- subset(input_data, 
                        Model == results$Model[i] & 
                          Previous_treatment == results$Previous_treatment[i] & 
                          Drug == results$Drug[i])
  
  # Fit isotonic regression model
  subset_data <- subset_data[order(subset_data$Concentration), ]
  model_fits[[i]] <- isoreg(subset_data$Concentration, subset_data$Response)
  
  # Perform bootstrapping to get variability estimates
  bootstrap_samples <- getBootSample(subset_data)
  lower_intervals[[i]] <- apply(bootstrap_samples, 1, quantile, probs=0.025)
  upper_intervals[[i]] <- apply(bootstrap_samples, 1, quantile, probs=0.975)
  
  # Calculate IC50 using linear interpolation
  approx_fit <- approx(log10(model_fits[[i]]$x), model_fits[[i]]$yf, n=100)
  IC50_value <- 10^(approx_fit$x[which.min(abs(approx_fit$y - 50))])
  
  # Calculate AUC using the flux package
  AUC_value <- auc(log10(model_fits[[i]]$x), model_fits[[i]]$yf, thresh=NULL, dens=100, sort.x=TRUE)
  scaled_AUC <- AUC_value / diff(range(log10(model_fits[[i]]$x)))
  
  # Store AUC and IC50 values in the results dataframe
  results$AUC[i] <- scaled_AUC
  results$IC50[i] <- IC50_value
  
  # Calculate confidence intervals for AUC
  bootstrap_AUC <- apply(bootstrap_samples, 2, function(x) {
    auc(log10(x_values), x, thresh=NULL, dens=100, sort.x=TRUE) / diff(range(log10(model_fits[[i]]$x)))
  })
  
  results$AUC_LI[i] <- quantile(bootstrap_AUC, probs=0.025)
  results$AUC_UI[i] <- quantile(bootstrap_AUC, probs=0.975)
  results$AUC_var[i] <- var(bootstrap_AUC)
  
  # Calculate confidence intervals for IC50
  bootstrap_IC50 <- apply(bootstrap_samples, 2, function(x) {
    tmp <- approx(log10(x_values), x, n=100)
    10^(tmp$x[which.min(abs(tmp$y - 50))])
  })
  
  results$IC50_LI[i] <- quantile(bootstrap_IC50, probs=0.025)
  results$IC50_UI[i] <- quantile(bootstrap_IC50, probs=0.975)
}

# Save the results to a CSV file
write.csv(results, "output_results.csv", row.names=FALSE)

# Plotting AUC values (Figure 6B)
ggplot(results, aes(x=Drug, y=AUC, col=Previous_treatment, ymin=AUC_LI, ymax=AUC_UI)) +
  geom_pointrange(position=position_dodge(width = 0.8), size=0.2) +
  ggtitle("AUC") +
  coord_flip() +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  theme_bw(base_size=14) +
  xlab("Drugs") +
  labs(colour="Cohort") +
  scale_color_manual(values=c("#E59F86", "#C83827","#4A9B87","#50BBD2","#8490B2","#000000")) +
  theme(plot.title=element_text(size=10), axis.title=element_text(size=10),
        axis.text=element_text(size=10), legend.text=element_text(size=10),
        legend.title=element_text(size=10)) +
  facet_wrap(~ Model)

# Save the plot
ggsave("AUC_plot.png")

# Plot dose-response curves (Figure 6C)
par(mfrow=c(2, 2))

# Example: Plotting for specific models/drugs - customize based on actual data
for (i in c(252, 236, 160, 144)) {
  plot(log10(model_fits[[i]]$x), model_fits[[i]]$yf, type="l", col="black", ylim=c(0, 100), 
       main=paste("Model", results$Model[i], ":", results$Drug[i]), 
       xlab="Concentration (uM)", ylab="Response (%)", xaxt="n")
  
  # Add confidence intervals as shaded polygons
  polygon(c(log10(x_values), rev(log10(x_values))), 
          c(upper_intervals[[i]], rev(lower_intervals[[i]])), 
          border=NA, col=rgb(0, 0, 0, 0.2))
  
  # Add the fitted line
  lines(log10(model_fits[[i]]$x), model_fits[[i]]$yf, col="black", lwd=2.5)
  
  # Define the x-axis with specific concentration labels
  axis(1, at=log10(x_values), labels=x_values)
  
  # Add a legend to differentiate between cohorts
  legend("topleft", cex=1, title="Cohort", legend=c("Untreated", "CT"), pch=20, col=c("black", "red"))
}

# Save the plot
dev.copy(png, "Dose_Response_Curves.png")
dev.off()

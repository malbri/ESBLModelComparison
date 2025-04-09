# load optimized model
library(deSolve)
library(dplyr)
library(ggplot2)

rm(list = ls())
#load(C:/Users/SIRModel")/Deterministic.RData")

results_df <- data.frame()

beta_f = Broilers.food + Cattle.food + Pig.food + Broilers.imp + Cattle.imp + Pig.imp + Turkey.imp
beta_p = Cat + Dog + Horse
beta_l = Broilers.livestock + Cattle.livestock + Pig.livestock

parms.list <- c(beta_h=beta_h, beta_f=beta_f, beta_p=beta_p, beta_l=beta_l, lambda=lambda, alpha,
                delta, rho, sigma, gamma, gamma_A, nu_P, mu_P,
                nu_G, mu_G, eta_f, eta_l, eta_p)



Output <- ode(y = init.list, times = time.points, func = ESBL_model, parms = parms.list, method="lsoda", maxsteps=100000)
Output_df <- as.data.frame(Output)


### Population sizes in the last year of simulation
# 1. Filter the data for the last year
result_last_day <- Output_df[which.max(Output_df$time),]

# 3. Calculate the total number of infected in I and IA over the last year (summing over all time points in the last year)
total_I_G <- result_last_day$I_G *365 *gamma
total_IA_G <- result_last_day$IA_G *365 *gamma_A

total_I_F <- result_last_day$I_F *365 *gamma
total_IA_F <- result_last_day$IA_F *365 *gamma_A

total_I_P <- result_last_day$I_P *365 *gamma
total_IA_P <- result_last_day$IA_P *365 *gamma_A

total_Infected<- total_I_G + total_IA_G + total_I_F + total_IA_F + total_I_P + total_IA_P

result <- data.frame(total_I_G, total_IA_G, total_I_F, total_IA_F, total_I_P, total_IA_P, total_Infected)

results_df <- rbind(results_df, result)
print(results_df)


names <- c("Baseline", "Broilers food", "Cattle food", "Pig food", 
           "Broilers import", "Cattle import", "Pig import", "Turkey import",
           "Cat", "Dog", "Horse",
           "Broilers livestock", "Cattle livestock", "Pig livestock",
           "Human","Travel",
           "All food", "All pets","All livestock")

results_df$Removed_variable <- names
results_new <- results_df[c(ncol(results_df), 1:(ncol(results_df) - 1))]

results_new$Category <- factor(c("Baseline", "Food", "Food", "Food", "Import", 
                          "Import", "Import", "Import", "Pets", "Pets", 
                          "Pets", "Livestock", "Livestock", "Livestock", "Human", "Human",
                          "Food", "Pets", "Livestock"), 
                        levels = c("Baseline", "Food", "Import", "Pets", "Livestock", "Human"))

#Calculate the difference from the baseline for each scenario
baseline_value <- results_new$total_Infected[1]
results_new <- results_new %>%
  mutate(
    difference = total_Infected - baseline_value,  # Set difference as NA for baseline
    pct_difference = (difference / baseline_value) * 100  # Set percentage as NA for baseline
  )

#load("C:/Users/SIRModel/Scenarios.RData")


#Plot 
library(scales) 
Source_removal <- ggplot(results_new, aes(x = total_Infected, y = reorder(Removed_variable, -total_Infected), fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Food" = "#f38dd4", "Import" = "#c1aaf6", "Pets" = "#8dd4f3", "Livestock" = "#f3ac8d", "Human" = "#a1f38d")) +
  theme_minimal() +
  labs(x = "Total Infected", y = "Source Removed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  
  scale_x_continuous(breaks = seq(0, max(results_new$total_Infected), by = 10000)) +
  
  #Add total infected labels
  geom_text(aes(label = ifelse(pct_difference == 0,
                               comma(ceiling(total_Infected)),  # Apply comma formatting
                               paste0(comma(ceiling(total_Infected)), " (", gsub("-", "−", round(pct_difference, 0)), "%)" ))), 
            position = position_stack(vjust = 0.5), 
            color = "black")
Source_removal 
#ggsave("C:/Users/SIRModel/Source_removal.png", plot = Source_removal)


library(writexl)
write_xlsx(results_new, "Scenarios.xlsx")

save.image("Scenarios.RData")

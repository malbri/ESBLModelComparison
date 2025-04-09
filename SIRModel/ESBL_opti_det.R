#restarting and setting working directory
rm(list = ls())
par(mfrow=c(1,1))
par(mar=c(5, 4, 4, 8), xpd=TRUE)
set.seed(123)
#setwd("C:/Users/SIRModel")

#loading libraries there are used
library(deSolve)
library(tidyverse)
library(ggplot2)
library(epiR)
library(cowplot)
library(dfoptim)
library(dplyr)  # Ensure dplyr is loaded for data manipulation

ESBL_model <- function(time, state_values, parameters) {
  
  # State_values
  # General
  S_G <- state_values[1]
  C_G <- state_values[2]
  I_G <- state_values[3]
  IA_G <- state_values[4]
  
  # Farmers
  S_F <- state_values[5]
  C_F <- state_values[6]
  I_F <- state_values[7]
  IA_F <- state_values[8]
  
  # Pet Owners
  S_P <- state_values[9]
  C_P <- state_values[10]
  I_P <- state_values[11]
  IA_P <- state_values[12]
  
  with(as.list(c(state_values, parameters)), {
    # All individuals
    N_G <- S_G + C_G + I_G + IA_G
    N_F <- S_F + C_F + I_F + IA_F
    N_P <- S_P + C_P + I_P + IA_P
    
    # Transmission between humans
    omega_G <- beta_h * ((C_G + I_G + IA_G) / N_G)
    omega_F <- beta_h * ((C_F + I_F + IA_F) / N_F)
    omega_P <- beta_h * ((C_P + I_P + IA_P) / N_P)
    
    # Common terms for transmission
    transmission_G <- omega_P + omega_F + omega_G + beta_f * eta_f
    transmission_F <- omega_P + omega_G + omega_F + beta_f * eta_f + beta_l * eta_l
    transmission_P <- omega_G + omega_F + omega_P + beta_f * eta_f + beta_p * eta_p
    
    # Common recovery terms
    recovery_S_G <- sigma * C_G + gamma * I_G + gamma_A * IA_G 
    recovery_S_F <- sigma * C_F + gamma * I_F + gamma_A * IA_F
    recovery_S_P <- sigma * C_P + gamma * I_P + gamma_A * IA_P 
    
    # Equations for General Population
    dS_G <- nu_G - S_G * transmission_G - lambda * S_G + recovery_S_G - mu_G * S_G
    dC_G <- S_G * transmission_G + lambda * S_G - (sigma + (1 - rho) * alpha + rho * theta * alpha + mu_G) * C_G
    dI_G <- alpha * C_G - (gamma + delta + mu_G) * I_G
    dIA_G <- delta * I_G + rho * theta * alpha * C_G - (gamma_A + mu_G) * IA_G
    
    # Equations for Farmers
    dS_F <- - S_F * transmission_F + recovery_S_F
    dC_F <- S_F * transmission_F - (sigma + (1 - rho) * alpha + rho * theta * alpha + mu_G) * C_F
    dI_F <- (1 - rho) * alpha * C_F - (gamma + delta) * I_F
    dIA_F <- delta * I_F + rho * theta * alpha * C_F - gamma_A * IA_F
    
    # Equations for Pet Owners
    dS_P <- nu_P - S_P * transmission_P + recovery_S_P - mu_P * S_P
    dC_P <- S_P * transmission_P - (sigma + (1 - rho) * alpha + rho * theta * alpha + mu_G) * C_P
    dI_P <- (1 - rho) * alpha * C_P - (gamma + delta + mu_P) * I_P
    dIA_P <- delta * I_P + rho * theta * alpha * C_P - (gamma_A + mu_P) * IA_P
    
    # Return a list containing derivatives
    list(c(dS_G, dC_G, dI_G, dIA_G, dS_F, dC_F, dI_F, dIA_F, dS_P, dC_P, dI_P, dIA_P))
    
  })
}

# Populations
N_F <- 24851
N_P <- 2024764 
N_G <- 5969774 - N_F - N_P  
N <- N_G + N_F + N_P

Carriers_G <- 0.037 * N_G
Carriers_F <- 0.11 * N_F
Carriers_P <- 0.042 * N_P

# Parameters
# Source transmissions
Broilers.food <-  1.5*112200/N/365    
Cattle.food <-  5.5*46447/N/365    
Pig.food <-  11.2*96024/N/365    
Broilers.imp <-  13.3*45492/N/365 
Cattle.imp <-  6.5*93687/N/365 
Pig.imp <-  30*195366/N/365 
Turkey.imp <-  52.2*6141/N/365 
Cat <- 4.3*730432/767866/365      
Dog <- 8.5*808519/1112283/365
Horse <- 10.5*174000/144616/365
Broilers.livestock <- 1.5*1462/1715/365
Cattle.livestock <- 5.5*8972/13077/365
Pig.livestock <- 11.2*2131/10059/365

#birth and death rates
nu_P = 157.45*0.339
mu_P = 3.36E-5*0.339
nu_G = 157.45-nu_P
mu_G = 3.36E-5-mu_P

optimPrevalence <- function(x, model, target) {
  # Update the parameter values with the ones being optimized
  parsValue <- list(
    beta_h = x[1],
    beta_f = Broilers.food + Cattle.food + Pig.food + Broilers.imp + Cattle.imp + Pig.imp + Turkey.imp,
    beta_p = Cat + Dog + Horse,
    beta_l = Broilers.livestock + Cattle.livestock + Pig.livestock,
    lambda = 0.0000264,
    alpha = 0.000690506,
    delta = 0,
    rho = (8.86 /1000),
    sigma = x[2],
    gamma = 0.0026,
    gamma_A = 0.0476,
    nu_P = nu_P,
    mu_P = mu_P,
    nu_G = nu_G,
    mu_G = mu_G,
    eta_f = x[3],
    eta_l = x[4],
    eta_p = x[5],
    theta = 3.7
  )
  
  init <- c(
    S_G = N_G * x[6],
    C_G = N_G * (1-x[6]),
    I_G = 0,
    IA_G = 0,
    S_F = N_F - Carriers_F,
    C_F = Carriers_F,
    I_F = 0,
    IA_F = 0,
    S_P = N_P * x[7],
    C_P = N_P * (1-x[7]),
    I_P = 0,
    IA_P = 0
  )
  # Solve the ODE model
  sol <- ode(y = init, times = seq(1, 5*365), func = model, parms = parsValue)
  
  # Extract relevant results and calculate the error
  Carriers_number <- as.data.frame(sol) %>% filter(time > 2*365) %>%
    mutate(month = month(as.Date(time, origin = "2023-01-01")),
           year = year(as.Date(time, origin = "2023-01-01")),
           week = ceiling(time / 7)) %>%
    group_by(week) %>%
    summarise(C_G = mean(C_G), C_F = mean(C_F), C_P = mean(C_P), I_T = mean(I_G)+mean(IA_G)+mean(I_P)+mean(IA_P)+mean(I_F)+mean(IA_F))
  
  
  err <- sum((target[1]*N_G - Carriers_number$C_G)^2 / N_G +
               (target[2]*N_F - Carriers_number$C_F)^2 / N_F +
               (target[3]*N_P - Carriers_number$C_P)^2 / N_P)#+
              #(target[3]-Carriers_number$I_T)^2/N)
  
  print(err)
  
  return(err)
}

# Initial parameter guess
initPar <- c(0.0053, 0.0002, 0.0079, 0.03, 0.03, 0.03, 0.03)
control_list <- list(maxfeval = 10000, tol = 1e-6, trace = 1)

# Optimization using nmkb
cOptim <- nmkb(par = initPar,
               lower = c(0.001, 0.0001, 0.001, 0.001, 0.001, 0.001, 0.001),
               upper = c(0.01,    0.9,     0.2,   0.2,      0.2,     1,     1),
               fn = optimPrevalence,
               model = ESBL_model,
               target = c(0.037,0.11, 0.042),
               control = control_list) 
print(cOptim)
save.image("Optim.RData")

############### Optimized model #######################
#load("C:/Users/SIRModel/Optim.RData")

beta_h = cOptim$par[1]
beta_f = Broilers.food + Cattle.food + Pig.food + Broilers.imp + Cattle.imp + Pig.imp + Turkey.imp
beta_p = Cat + Dog + Horse
beta_l = Broilers.livestock + Cattle.livestock + Pig.livestock
lambda = 0.0000264
alpha = 0.000690506
delta = 0
rho = (8.86 /1000)
sigma = cOptim$par[2]
gamma = 0.0026
gamma_A = 0.0471
eta_f = cOptim$par[3]
eta_l = cOptim$par[4]
eta_p = cOptim$par[5]
theta = 3.7

parms.list <- c(beta_h, beta_f, beta_p, beta_l, lambda, alpha,
           delta, rho, sigma, gamma, gamma_A, nu_P, mu_P,
           nu_G, mu_G, eta_f, eta_l, eta_p, theta)

# Initial conditions
init.list <- c(
  S_G = N_G * cOptim$par[6],
  C_G = N_G * (1-cOptim$par[6]),
  I_G = 0,
  IA_G = 0,
  S_F = N_F - Carriers_F,
  C_F = Carriers_F,
  I_F = 0,
  IA_F = 0,
  S_P = N_P * cOptim$par[7],
  C_P = N_P * (1-cOptim$par[7]),
  I_P = 0,
  IA_P = 0
)

# Simulation  
k <- 5                                       # Years to simulate
time <- 365*k + 1 
# Time points
time.points <- seq(0, time, by=1)

Output <- ode(y = init.list, times = time.points, func = ESBL_model, parms = parms.list, method="lsoda", maxsteps=100000)
Output_df <- as.data.frame(Output)


# Empty plot to plot simulations in
plot(NULL, xlim = c(0, time), ylim = c(0, 2200), ylab = "Number of individuals", xlab = "Days", main = "Transmission of ESBL E. coli in Denmark")

legend("topright", inset = c(-0.3, 0), 1:10, legend = c("General", "Farmers", "Pet Owners", "General", "Farmers", "Pet Owners", "Infected", "Infected_antibiotics"),
       col = c("lightblue", "lightgreen", "lightcoral", "lightblue", "lightgreen", "lightcoral", "grey", "grey"), 
       lty = c(1, 1, 1, 2, 2, 2, 1, 2), cex = 0.6, lwd = 2)

# Plotting the simulation
lines(time.points, Output_df$I_G, col = "lightblue", lwd = 2)
lines(time.points, Output_df$IA_G, col = "lightblue", lty = 2, lwd = 2)

lines(time.points, Output_df$I_F, col = "lightgreen", lwd = 2)
lines(time.points, Output_df$IA_F, col = "lightgreen", lty = 2, lwd = 2)

lines(time.points, Output_df$I_P, col = "lightcoral", lwd = 2)
lines(time.points, Output_df$IA_P, col = "lightcoral", lty = 2, lwd = 2)

## OBS husk at gemme ##


##################### Check Prevalence Carriage ###################
# 1. Calculate the total population for each population (G, F, P) at each time step
Output_df$total_G <- Output_df$S_G + Output_df$C_G + Output_df$I_G + Output_df$IA_G
Output_df$total_F <- Output_df$S_F + Output_df$C_F + Output_df$I_F + Output_df$IA_F
Output_df$total_P <- Output_df$S_P + Output_df$C_P + Output_df$I_P + Output_df$IA_P

# 2. Calculate the percentage of C for each population at each time step
Output_df$percent_C_G <- (Output_df$C_G / Output_df$total_G) * 100
Output_df$percent_C_F <- (Output_df$C_F / Output_df$total_F) * 100
Output_df$percent_C_P <- (Output_df$C_P / Output_df$total_P) * 100

result_prevalence <- Output_df[, c("time", "percent_C_G", "percent_C_F", "percent_C_P")]
tail(result_prevalence)

result_last_day <- Output_df[which.max(Output_df$time),]
# 3. Calculate the total number of infected in I and IA over the last year (summing over all time points in the last year)
total_I_G <- result_last_day$I_G *365 *gamma
total_IA_G <- result_last_day$IA_G *365 *gamma_A

total_I_F <- result_last_day$I_F *365 *gamma
total_IA_F <- result_last_day$IA_F *365 *gamma_A

total_I_P <- result_last_day$I_P *365 *gamma
total_IA_P <- result_last_day$IA_P *365 *gamma_A

total_Infected <- total_I_G + total_IA_G + total_I_F + total_IA_F + total_I_P + total_IA_P

result_df <- data.frame(total_I_G, total_IA_G, total_I_F, total_IA_F, total_I_P, total_IA_P, total_Infected)
print(result_df)

total_C_G <- result_last_day$C_G *365 *sigma
total_C_F <- result_last_day$C_F *365 *sigma
total_C_P <- result_last_day$C_P *365 *sigma

total_C_all <- total_C_G + total_C_F + total_C_P

results_colonization <- data.frame(total_C_G,total_C_F,total_C_P, total_C_all)
print(results_colonization)

baseline_result <- data.frame(total_C_G,total_C_F,total_C_P, total_C_all,total_I_G, total_IA_G, total_I_F, total_IA_F, total_I_P, total_IA_P, total_Infected)
library(writexl)
write_xlsx(baseline_result, "Baseline_final.xlsx")

# Reshape the data for plotting
result_long <- data.frame(
  population = rep(c("General", "Farmers", "Pet Owners"), each = 2),
  infection_type = rep(c("Infected", "Infected with antibiotic usage"), times = 3),
  total = c(result_df$total_I_G, result_df$total_IA_G,
            result_df$total_I_F, result_df$total_IA_F,
            result_df$total_I_P, result_df$total_IA_P)
)

# Plotting
Infected_by_population <- ggplot(result_long, aes(x = population, y = total, fill = infection_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Infected" = "skyblue", 
                               "Infected with antibiotic usage" = "salmon")) +
  labs(title = "Infection by Population",
       x = "Population",
       y = "Total Count",
       fill = "Infection Type") +
  #geom_text(aes(label = round(total)), 
            #position = position_stack(vjust = 0.5), 
            #color = "black", size = 3) +
  theme_minimal()

#ggsave("C:/Users/SIRModel/Infected_by_population.png", plot = Infected_by_population)

save.image("Deterministic.RData")


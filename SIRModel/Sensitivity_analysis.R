rm(list = ls())
#load("C:/Users/SIRModel/Deterministic.RData")

######################################### Sensitivity Analysis ######################################### 
# Finding new intervals for the ODE model parameters

min_beta_h <- beta_h - beta_h * 0.25  
max_beta_h <- beta_h + beta_h * 0.25 

min_beta_f <- beta_f - beta_f * 0.25  
max_beta_f <- beta_f + beta_f * 0.25 

min_beta_p <- beta_p - beta_p * 0.25  
max_beta_p <- beta_p + beta_p * 0.25 

min_beta_l <- beta_l - beta_l * 0.25  
max_beta_l <- beta_l + beta_l * 0.25 

min_lambda <- lambda - lambda * 0.25  
max_lambda <- lambda + lambda * 0.25 

min_alpha <- alpha - alpha * 0.25  
max_alpha <- alpha + alpha * 0.25 

min_delta <- delta - delta * 0.25  
max_delta <- delta + delta * 0.25 

min_rho <- rho - rho * 0.25  
max_rho <- rho + rho * 0.25 

min_sigma <- sigma - sigma * 0.25  
max_sigma <- sigma + sigma * 0.25 

min_gamma <- gamma - gamma * 0.25  
max_gamma <- gamma + gamma * 0.25 

min_gamma_A <- gamma_A - gamma_A * 0.25  
max_gamma_A <- gamma_A + gamma_A * 0.25 

min_eta_f <- eta_f - eta_f * 0.25
max_eta_f <- eta_f + eta_f * 0.25

min_eta_l <- eta_l - eta_l * 0.25
max_eta_l <- eta_l + eta_l * 0.25

min_eta_p <- eta_p - eta_p * 0.25
max_eta_p <- eta_p + eta_p * 0.25

min_theta <- theta - theta * 0.25
max_theta <- theta + theta * 0.25

# Scatterplot
time_steps <- 365 * 10
N.iter <- 1000
x <- 365 * 10

# Generate random parameter values based on the sensitivity analysis
lambda.list <- runif(N.iter, min=min_lambda, max=max_lambda)
beta_f.list <- runif(N.iter, min=min_beta_f, max=max_beta_f)
beta_h.list <- runif(N.iter, min=min_beta_h, max=max_beta_h)
beta_l.list <- runif(N.iter, min=min_beta_l, max=max_beta_l)
beta_p.list <- runif(N.iter, min=min_beta_p, max=max_beta_p)
sigma.list <- runif(N.iter, min=min_sigma, max=max_sigma)
alpha.list <- runif(N.iter, min=min_alpha, max=max_alpha)
rho.list <- runif(N.iter, min=min_rho, max=max_rho)
gamma.list <- runif(N.iter, min=min_gamma, max=max_gamma)
gamma_A.list <- runif(N.iter, min=min_gamma_A, max=max_gamma_A)
delta.list <- runif(N.iter, min=min_delta, max=max_delta)
eta_f.list <- runif(N.iter, min=min_eta_f, max=max_eta_f)
eta_l.list <- runif(N.iter, min=min_eta_l, max=max_eta_l)
eta_p.list <- runif(N.iter, min=min_eta_p, max=max_eta_p)
theta.list <- runif(N.iter, min=min_theta, max=max_theta)

# Initialize storage vectors for results
S_G.pro <- numeric(N.iter)
C_G.pro <- numeric(N.iter)
I_G.pro <- numeric(N.iter)
IA_G.pro <- numeric(N.iter)

S_F.pro <- numeric(N.iter)
C_F.pro <- numeric(N.iter)
I_F.pro <- numeric(N.iter)
IA_F.pro <- numeric(N.iter)

S_P.pro <- numeric(N.iter)
C_P.pro <- numeric(N.iter)
I_P.pro <- numeric(N.iter)
IA_P.pro <- numeric(N.iter)


for ( i in 1:N.iter){
  print(i)
  # simulate the epidemic
  parameter.list <- c(lambda = lambda.list[i], beta_f = beta_f.list[i], beta_h = beta_h.list[i], 
                      beta_l = beta_l.list[i], beta_p = beta_p.list[i], sigma = sigma.list[i], 
                      alpha = alpha.list[i], rho = rho.list[i], gamma = gamma.list[i], 
                      gamma_A = gamma_A.list[i], delta = delta.list[i], eta_f = eta_f.list[i], 
                      eta_l = eta_l.list[i], eta_p = eta_p.list[i], theta = theta.list[i])
  
  # Start values
  S_G_start <- N_G * cOptim$par[6]
  C_G_start <- N_G * (1-cOptim$par[6])
  I_G_start <- 0
  IA_G_start <- 0
  
  S_F_start <- N_F - Carriers_F
  C_F_start <- Carriers_F
  I_F_start <- 0
  IA_F_start <- 0
  
  S_P_start <- N_P * cOptim$par[7]
  C_P_start <- N_P * (1-cOptim$par[7])
  I_P_start <- 0
  IA_P_start <- 0
  
  initial.values <- c(S_G = S_G_start, C_G = C_G_start, I_G = I_G_start, IA_G = IA_G_start,
                      S_F = S_F_start, C_F = C_F_start, I_F = I_F_start, IA_F = IA_F_start,
                      S_P = S_P_start, C_P = C_P_start, I_P = I_P_start, IA_P = IA_P_start)
  
  #setting time points
  time.points <- seq(0,time_steps,by=1)
  
  output <- ode(y = initial.values, times = time.points, func = ESBL_model, parms = parameter.list)
  
  # Extract and store results (assuming SCIS_model returns correct outputs)
  S_G.pro[i] <- output[time_steps, "S_G"]
  C_G.pro[i] <- output[time_steps, "C_G"]
  I_G.pro[i] <- output[time_steps, "I_G"]
  IA_G.pro[i] <- output[time_steps, "IA_G"]
  
  S_F.pro[i] <- output[time_steps, "S_F"]
  C_F.pro[i] <- output[time_steps, "C_F"]
  I_F.pro[i] <- output[time_steps, "I_F"]
  IA_F.pro[i] <- output[time_steps, "IA_F"]
  
  S_P.pro[i] <- output[time_steps, "S_P"]
  C_P.pro[i] <- output[time_steps, "C_P"]
  I_P.pro[i] <- output[time_steps, "I_P"]
  IA_P.pro[i] <- output[time_steps, "IA_P"]
  
}


# Data frame with result from the simulation. Further divided into to sub groups
ESBL.res <- data.frame(lambda = lambda.list, beta_f = beta_f.list, beta_h = beta_h.list,
                       beta_l = beta_l.list, beta_p = beta_p.list, sigma = sigma.list, 
                       alpha = alpha.list, rho = rho.list, gamma = gamma.list, 
                       gamma_A = gamma_A.list, eta_f = eta_f.list, 
                       eta_l = eta_l.list, eta_p = eta_p.list, theta = theta.list,
                       S_G.pro, C_G.pro, I_G.pro, IA_G.pro, 
                       S_F.pro, C_F.pro, I_F.pro, IA_F.pro, 
                       S_P.pro, C_P.pro, I_P.pro, IA_P.pro)

# Subset for General Population
ESBL.general <- ESBL.res[c("lambda", "beta_f", "beta_h", "beta_l", "beta_p", "sigma", "alpha", "rho",
                           "gamma", "gamma_A", "eta_f", "eta_l", "eta_p", "theta",
                           "S_G.pro", "C_G.pro", "I_G.pro", "IA_G.pro")]

# Subset for Farmers
ESBL.farmers <- ESBL.res[c("lambda", "beta_f", "beta_h", "beta_l", "beta_p", "sigma", "alpha", "rho",
                           "gamma", "gamma_A", "eta_f", "eta_l", "eta_p","theta",
                           "S_F.pro", "C_F.pro", "I_F.pro", "IA_F.pro")]

# Subset for Pet Owners
ESBL.pets <- ESBL.res[c("lambda", "beta_f", "beta_h", "beta_l", "beta_p", "sigma", "alpha", "rho",
                        "gamma", "gamma_A", "eta_f", "eta_l", "eta_p","theta",
                        "S_P.pro", "C_P.pro", "I_P.pro", "IA_P.pro")]

# Scatter plots - looking at different settings
pairs(ESBL.general, upper.panel = NULL, pch = 19, cex = 0.1, col = "lightblue")
pairs(ESBL.farmers,upper.panel = NULL, pch = 19, cex = 0.1, col = "lightgreen")
pairs(ESBL.pets, upper.panel = NULL, pch = 19, cex = 0.1, col = "lightpink")

# PRCC Calculation - in relation to C?
prcc.general <- epi.prcc(ESBL.general[,c(1:14,16)], sided.test = 2, conf.level = .95)
# Sort the data by actual PRCC values (not the absolute ones)
prcc.general$abs_value <- abs(prcc.general[,2])
prcc.general <- prcc.general[order(prcc.general$abs_value),]
prcc.general$var <- factor(prcc.general$var, levels = prcc.general$var)

# Create a named vector that maps variables to their Greek symbols
variable_to_greek <- c(
  "lambda" = expression(lambda),
  "beta_f" = expression(beta[f]),
  "beta_h" = expression(beta[h]),
  "beta_l" = expression(beta[l]),
  "beta_p" = expression(beta[p]),
  "sigma" = expression(sigma),
  "alpha" = expression(alpha),
  "rho" = expression(rho),
  "gamma" = expression(gamma),
  "gamma_A" = expression(gamma[A]),
  "eta_f" = expression(eta[f]),
  "eta_l" = expression(eta[l]),
  "eta_p" = expression(eta[p]),
  "theta" = expression(theta)
)

# Ensure the labels follow the reordering
gen_parameter.names <- sapply(prcc.general$var, function(x) variable_to_greek[as.character(x)])

# Create the tornado plot with Greek symbols and increased text size
prcc_general_plot <- ggplot(data = prcc.general, aes(x = var, y = prcc.general[,2])) +
  geom_bar(stat = "identity", fill = "lightblue", width = 0.7) +
  geom_text(aes(label = round(prcc.general[,2], 2)), 
            hjust = ifelse(prcc.general[,2] > 0, -0.2, 1.2), size = 3.5) +  # Adjust text size here
  theme_classic() +
  labs(y = "PRCC - General Population", x = "") +
  ylim(-1, 1) +
  ggtitle(" General Population") +
  coord_flip() +  # Flip coordinates to make it horizontal
  scale_x_discrete(labels = gen_parameter.names) +  # Use reordered Greek symbols for labels
  theme(
    axis.text = element_text(size = 12),  # Adjust axis text size
    plot.title = element_text(size = 18, face = "bold"),  # Adjust title size
    axis.title.y = element_text(size = 14)  # Adjust y-axis label size
  )

# Print the plot
print(prcc_general_plot)

#ggsave("C:/Users/SIRModel/prcc_general_plot.png", plot = prcc_general_plot)

# PRCC for Farmers
prcc.farmers <- epi.prcc(ESBL.farmers[,c(1:14,16)], sided.test = 2, conf.level = .95)

# Sorting the data by the absolute value of PRCC for farmers
prcc.farmers$abs_value <- abs(prcc.farmers[,2])
prcc.farmers <- prcc.farmers[order(prcc.farmers$abs_value),]
prcc.farmers$var <- factor(prcc.farmers$var, levels = prcc.farmers$var)
far_parameter.names <- sapply(prcc.farmers$var, function(x) variable_to_greek[as.character(x)])

# Create the tornado plot for Farmers with Greek symbols and increased text size
prcc_farmers_plot <- ggplot(data = prcc.farmers, aes(x = var, y = prcc.farmers[,2])) +
  geom_bar(stat = "identity", fill = "lightgreen", width = 0.7) +
  geom_text(aes(label = round(prcc.farmers[,2], 2)), 
            hjust = ifelse(prcc.farmers[,2] > 0, -0.2, 1.2), size = 3.5) +  # Adjust text size here
  theme_classic() +
  labs(y = "PRCC - Farmers", x = "") +
  ylim(-1, 1) +
  ggtitle("Farmers") +
  coord_flip() +  # Flip coordinates to make it horizontal
  scale_x_discrete(labels = far_parameter.names) +  # Use reordered Greek symbols for labels
  theme(
    axis.text = element_text(size = 12), 
    plot.title = element_text(size = 18, face = "bold"), #
    axis.title.y = element_text(size = 14)  # Adjust y-axis label size
  )
plot(prcc_farmers_plot)
# Save the tornado plot for farmers
#ggsave("C:/Users/SIRModel/prcc_farmers_plot.png", plot = prcc_farmers_plot)

# Pets #
prcc.pets <- epi.prcc(ESBL.pets[,c(1:14,16)], sided.test = 2, conf.level = .95)

# Sorting the data by the absolute value of PRCC for pet owners
prcc.pets$abs_value <- abs(prcc.pets[,2])
prcc.pets <- prcc.pets[order(prcc.pets$abs_value),]
prcc.pets$var <- factor(prcc.pets$var, levels = prcc.pets$var)
pet_parameter.names <- sapply(prcc.pets$var, function(x) variable_to_greek[as.character(x)])

# Create the tornado plot for pet owners
prcc_pets_plot <- ggplot(data = prcc.pets, aes(x = var, y = prcc.pets[,2])) +
  geom_bar(stat = "identity", fill = "lightpink", width = 0.7) +
  geom_text(aes(label = round(prcc.pets[,2], 2)), 
            hjust = ifelse(prcc.pets[,2] > 0, -0.2, 1.2), size = 3.5) +  # Adjust text size here
  theme_classic() +
  labs(y = "PRCC - Pet Owners", x = "") +
  ylim(-1, 1) +
  ggtitle("Pet Owners") +
  coord_flip() +  # Flip coordinates to make it horizontal
  scale_x_discrete(labels = pet_parameter.names) +  # Use reordered Greek symbols for labels
  theme(
    axis.text = element_text(size = 14),  # Adjust axis text size
    plot.title = element_text(size = 18, face = "bold"),  # Adjust title size
    axis.title.y = element_text(size = 14)  # Adjust y-axis label size
  )
plot(prcc_pets_plot)

# Save the tornado plot for pet owners
#ggsave("C:/Users/SIRModel/prcc_pets_plot.png", plot = prcc_pets_plot)

# Combine the plots into a grid
plot_grid(prcc_general_plot, prcc_farmers_plot, prcc_pets_plot, labels = "AUTO")

save.image(file="Sensitivity_analysis")


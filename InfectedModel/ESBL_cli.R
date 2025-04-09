rm(list = ls())
#set working directory
#setwd("C:/Users/InfectedModel")

# Load packages
library(reshape2)
library(tidyverse)
library(writexl)
library(dplyr)
library(readxl)
library(coda)

set.seed(123)


#make Genenames
All_data <- read_excel("All_data.xlsx")
overall_prev <- read_excel("overall_prev.xlsx")
Genes_list <- read_excel("gene_list.xlsx")$Genes

Genes_data <- All_data #%>% 
  #filter (Genes %in% Genes_list) 

#split data to human and source data
human_data <- Genes_data[Genes_data$Source == "Human_CL", ]
Updated_data <- Genes_data[Genes_data$Source != "Human_CL", ]

# Count the number of each Gene within a source
Updated_data <- Updated_data %>% group_by(Source, Genes) %>% summarise(Count = n())

# Find the amount of source j with resistance genes
bd.qtts <- Updated_data %>% group_by(Source) %>% summarise(Count = sum(Count))

# Find prevalence of Gene i in source j (r)
total_counts <- Updated_data %>% group_by(Genes, Source) %>% summarise(Total_Count = sum(Count)) 
source_totals <- total_counts %>% group_by(Source) %>% summarise(Total_Count = sum(Total_Count))
prevalence_data <- total_counts %>% left_join(source_totals, by = "Source") %>%
  mutate(Prevalence = Total_Count.x / Total_Count.y) %>% select(Genes, Source, Prevalence)

# Create a grid of all possible Gene-source combinations
all_Genes <- unique(Updated_data$Genes)
all_sources <- unique(Updated_data$Source)

bd.prev <- data.frame(Genes = rep(all_Genes, each = length(all_sources)),
                      Source = rep(all_sources, times = length(all_Genes))) %>%
  left_join(prevalence_data, by = c("Genes", "Source")) %>%
  mutate(Prevalence = ifelse(is.na(Prevalence), 0, Prevalence)) %>%
  arrange(Source, Genes)
bd.prev <- bd.prev %>% filter (Genes %in% Genes_list) 
#write_xlsx(bd.prev, "prev.xlsx")

# Find the number of cases/year of gene i
human_data <- human_data %>% group_by(Source, Genes) %>% summarise(Count = n())
bd.obs <- aggregate(Count ~ Genes, data = human_data, sum) %>% filter(Genes %in% unique(bd.prev$Genes))

## model work ##
library(rjags)

writeLines(
  "model{
  for(i in 1:ngene){
    for(j in 1:nsource){
  
      prevs[i, j] <- r[i,j]*pi[j]
  
      lambdaij[i,j] <- prevs[i, j]*a[j]*q[i] 
    }  
    cases[i] ~ dpois(lambdai[i])
    lambdai[i] <- sum(lambdaij[i, 1:nsource])
    }
  
  for (j in 1:nsource) { 
    r[1:ngene, j] ~ ddirch(chi[,j])
  }
  
  for (j in 1:nsource) {
    pi[j] ~ dbeta(alpha[j], beta[j] - alpha[j])
  }
  
  a[1]~dexp(0.002)
  a[2]=a[1]
  a[3]~dexp(0.002)
  a[4]~dexp(0.002)
  a[5]=a[4]
  a[6]~dexp(0.002)
  a[7]~dexp(0.002)
  a[8]~dexp(0.002)
  a[9]~dexp(0.002)
  a[10]=a[9]
  a[11]~dexp(0.002)


  for(j in 1:nsource){
    lambdaj[j] <- sum(lambdaij[1:ngene, j])
  }
  
  for(i in 1:ngene){
    q[i]~dlnorm(0, tau)
  }
  tau ~dgamma(0.01,0.01)
  }
  ", "ESBL_cli.txt")

# Defining the data as a list, like BUGS
ngene = bd.prev %>%summarize(n_distinct(Genes)) #n Genes (?ndice i) 
ngene = ngene$`n_distinct(Genes)`
nsource = bd.prev  %>%summarize(n_distinct(Source)) #n animal source (?ndice j) 
nsource =nsource$`n_distinct(Source)`

chi <- matrix(rep(0.5, ngene * nsource), ncol = nsource)

jagsdata = list(ngene = bd.obs %>% summarize(n_distinct(Genes)),
                nsource = bd.prev  %>%summarize(n_distinct(Source)),
                cases = array(bd.obs$Count),
                r = array(bd.prev$Prevalence, c(ngene, nsource)),
                chi = chi,
                alpha = overall_prev$Positive,
                beta = overall_prev$Tested)



# Defining the initials as a list
inits = list(a=rep(0.002,nsource),tau=0.0001)
inits=list(inits,inits,inits)

# Defining burn-in, number of simulations, thin and parameters for which you want results
burn = 30000
nsim = 100000
nthin = 10
parms = c("a","q", "lambdaij", "lambdai", "lambdaj")


# Starting the model with a chain and 1000 iterations for adaptation
m <- jags.model("ESBL_cli.txt", jagsdata, inits, n.chains=3,n.adapt=1000)

# Updating burn-in
update(m,burn)

# Getting final sample
mcmc <- coda.samples(m, parms, n.iter=nsim,thin=nthin)

# Results - summary and graphs
summary = summary(mcmc) 
par(mar = c(1, 1, 1, 1))
pdf(file= "plot_CL.pdf" )
plot(mcmc)
#gelman.plot(mcmc)

dev.off()

library(coda)
library(MCMCvis)
lambdaij <- MCMCchains(mcmc, params = 'lambdaij')
lambdaij1<-apply(lambdaij, 2, median)
lambdaij2=matrix(lambdaij1, ncol=11, nrow=22,
                 dimnames = list(bd.obs$Genes,
                                 c("Broilers", "Broilers_IMP", "Cat", "Cattle", "Cattle_IMP", "Dog","Horse", "Human_OC", "Pig", "Pig_IMP", "Turkey_IMP")))
summary(lambdaij2)

#observed
pred.obs <- cbind(bd.obs,lambdaij2)
pred.obs <- as.data.frame(pred.obs)
pred.obs <- pred.obs %>%
  mutate(Not_predicted = Count - (Broilers + Broilers_IMP + Cat + Cattle + Cattle_IMP + Dog + Horse + Human_OC + Pig + Pig_IMP + Turkey_IMP)) %>%
  mutate(Not_predicted = pmax(Not_predicted, 0)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

pred.obs <- pred.obs[,-1:-2]
pred.obs$Gene <- rownames(pred.obs)
results <- pred.obs %>%
  pivot_longer(cols = c(-Gene), names_to = "Source", values_to = "Value")

# calculate percentage
results_percentage <- results %>%
  group_by(Gene) %>%
  mutate(Total = sum(Value)) %>%  # Calculate the total value per Gene
  ungroup() %>%
  mutate(Percentage = (Value / Total) * 100) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1))) %>%
  mutate(Group = case_when(
    Source %in% c("Broilers", "Cattle", "Pig") ~ "Food",
    Source %in% c("Cat", "Dog", "Horse") ~ "Pets",
    Source == "Human_OC" ~ "Human",
    Source %in% c("Broilers_IMP", "Cattle_IMP", "Turkey_IMP", "Pig_IMP") ~ "Import",
    Source == "Not_predicted" ~ "Not_predicted",
    TRUE ~ "Other"  
  ))

# Modify factor levels without expression
results_percentage$Source <- factor(results_percentage$Source,
                                    levels = c("Broilers", "Cattle", "Pig",  # Food group
                                               "Cat", "Dog", "Horse",        # Pets group
                                               "Human_OC",                   
                                               "Broilers_IMP", "Cattle_IMP", "Turkey_IMP", "Pig_IMP", # Import group without expression
                                               "Not_predicted")  
)

# Define distinct colors for each Source (same as before)
color_palette <- c(
  "Broilers" = "#f9c6ea",    
  "Cattle" = "#f6aadf",     
  "Pig" = "#f38dd4",       
  "Cat" = "#c6eaf9",         
  "Dog" = "#aadff6",          
  "Horse" = "#8dd4f3",        
  "Human_OC" = "#adf49b",     
  "Broilers_IMP" = "#f0eafd", 
  "Cattle_IMP" = "#e0d4fb",   
  "Pig_IMP" = "#d1bff8",   
  "Turkey_IMP" = "#c1aaf6",
  "Not_predicted" = "#afb1af"
)

custom_y_labels <- c("0%", "20%", "40%", "60%", "80%", "100%")

#Plot with superscript labels for imports
SA_CL_gene_specific <- ggplot(results_percentage, aes(x = Gene, y = Percentage, fill = Source)) +
  geom_bar(stat = "identity", ) +
  scale_y_continuous(
    breaks = c(0, 20, 40, 60, 80, 100),  # Define breaks for the y-axis
    labels = custom_y_labels
  ) +
  scale_fill_manual(values = color_palette, 
                    labels = c("Broilers", "Cattle", "Pig",  # Food group
                               "Cat", "Dog", "Horse",        # Pets group
                               "Human Colonized",                   
                               "Broilers Import", "Cattle Import", "Turkey Import", "Pig Import",  #
                               "Not predicted")) +            # Import group labels with superscript
  labs(y = "Relative probability", x = "ESBL and pAmpC genes") +  # Label the axes
  theme_minimal() +  # Use a minimal theme for the plot
  theme(axis.text.x = element_text(hjust = 1),
        legend.position = "bottom") +
  #guides(fill = guide_legend(nrow = 2)) +
  coord_flip()
# Print plot
SA_CL_gene_specific

# Define the sources for sorting
import_sources <- c("Broilers_IMP", "Cattle_IMP", "Turkey_IMP", "Pig_IMP")
human_source <- "Human_OC"
pet_sources <- c("Cat", "Dog", "Horse")  # Combine Cat, Dog, and Horse as Pets

# Calculate the sum of relevant percentages for sorting
gene_sorting <- results_percentage %>%
  group_by(Gene) %>%
  summarize(
    Turkey_pct = sum(Percentage[Source == "Turkey_IMP"]),
    Import_pct = sum(Percentage[Source %in% import_sources]),
    Human_pct = sum(Percentage[Source == human_source]),
    Pets_pct = sum(Percentage[Source %in% pet_sources])  # Combine Cat, Dog, Horse as Pets
  ) %>%
  # Create a sorting column based on the thresholds
  mutate(Sort_Order = case_when(
    Turkey_pct > 40 ~ 1,   # Turkey > 40%
    Import_pct > 45 ~ 2,   # Other Import > 40%
    Human_pct > 30 ~ 3,    # Human > 40%
    Pets_pct > 40 ~ 4,     # Pets (Cat, Dog, Horse combined) > 40%
    TRUE ~ 5               # All others
  )) %>%
  # Internal sorting within each group based on the relevant percentage
  arrange(Sort_Order, desc(Turkey_pct), desc(Import_pct), desc(Human_pct), desc(Pets_pct))

# Reorder the genes in the results_percentage based on this new sorting
results_percentage$Gene <- factor(results_percentage$Gene, levels = gene_sorting$Gene)

# Plot again with the sorted genes
SA_CL_gene_specific_sorted <- ggplot(results_percentage, aes(x = Gene, y = Percentage, fill = Source)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    breaks = c(0, 20, 40, 60, 80, 100),  # Define breaks for the y-axis
    labels = custom_y_labels
  ) +
  scale_fill_manual(values = color_palette, 
                    labels = c("Broilers", "Cattle", "Pig",  # Food group
                               "Cat", "Dog", "Horse",        # Pets group
                               "Human Colonized",                   
                               "Broilers Import", "Cattle Import", "Turkey Import", "Pig Import",  #
                               "Not predicted")) +            # Import group labels with superscript
  labs(y = "Relative probability", x = "ESBL and pAmpC genes") +  # Label the axes
  theme_minimal() +  # Use a minimal theme for the plot
  theme(axis.text.x = element_text(hjust = 1),
        legend.position = "bottom") +
  coord_flip()

# Print the sorted plot
print(SA_CL_gene_specific_sorted)

#ggsave("C:/Users/InfectedModel/SA_CL_gene_specific.png", plot = SA_CL_gene_specific_sorted)

write_xlsx(pred.obs, "obs_pred_counts_CL.xlsx")


#model output
parameter_names <- rownames(summary$statistics)
means <- round(summary$statistics[, "Mean"], 1)
medians <- round(summary$quantiles[, "50%"], 1)
lower_quantiles <- round(summary$quantiles[, "2.5%"], 1)
upper_quantiles <- round(summary$quantiles[, "97.5%"], 1)
SD <- round(summary$statistics[,"SD"], 1)

summary_df <- data.frame(Parameter = parameter_names,
                         Mean = means,
                         Median = medians,
                         Lower_Quantile = lower_quantiles,
                         Upper_Quantile = upper_quantiles,
                         SD = SD)
write_xlsx(summary_df, "model_output.xlsx")

output <- summary_df[grep("^lambdaj\\[",summary_df$Parameter),]
output$Source <- all_sources

output$Lower_Error <- output$Mean - output$Lower_Quantile
output$Upper_Error <- output$Upper_Quantile - output$Mean
#output$Source <- factor(output$Source, levels = output$Source[order(-output$Mean)])

custom_labels <- c("Broiler", "Broiler Import", "Cat", "Cattle", "Cattle Import", "Dog", "Horse", "Human colonized", "Pig", "Pig Import", "Turkey Import")

SA_CL_output <- ggplot(output, aes(x = Source, y = Mean)) +
  geom_bar(stat = "identity",  fill = "lightblue", color = "black") +
  geom_errorbar(aes(ymin = Mean - Lower_Error, ymax = Mean + Upper_Error), 
                width = 0.2) +
  labs(title = "Mean ESC-EC cases attributed to different sources",
       x = "Source of contamination",
       y = "Mean") +
  scale_x_discrete(labels = custom_labels) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
SA_CL_output

#ggsave("C:/Users/InfectedModel/SA_CL_output.png", plot = SA_CL_output)

#medians = sumli_Gene$quantiles[ ,'50%']
#plot(bd.obs$Count[1:16], medians)
#lines(0:4000,0:4000)

#print(bd.prev)
#remove mcmc to save and open workspace again
rm(mcmc)

save.image(file="DK_CL")



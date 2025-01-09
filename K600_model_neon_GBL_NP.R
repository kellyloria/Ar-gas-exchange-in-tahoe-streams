## Merges clean MIMS data with injection meta data
## Re-worked from Hall and Madinger: Use of argon to measure gas exchange in turbulent mountain streams

## CODE for not pooled sites. 
## So model the sites one at a time:

library(rstan)
library(shinystan)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# daily file info for saving:
set.seed(2021)
rundate <- format(Sys.Date(), "%y%m%d")
file_name <- "ArN_NEON_dist_Kt_GBL"

setwd("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/")
source("/Users/kellyloria/Documents/UNR/Reaeration/AR_code_repo/mims_gas_functions_wHeKr.R")
###

## 1. Read in processed MIMS data for Ar 
files <- list.files(paste("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/processed_dat/Merged_processed_dat",sep=""), full.names = T)
rawdat <-  do.call(rbind, lapply
                   (files, read.csv, as.is=T, header = T))
str(rawdat)
## 2. meta data for each trial location, time,  station distance etc.
meta_dat <- read.csv("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/processed_dat/MIMS_SampleLog_24.csv") %>%
  mutate(SampleID=as.character(sampleID))
str(meta_dat)

## Fix stream velocity to be in m^2 per day
meta_dat$v <- c(meta_dat$v * 86400)

## 3. Merge by SampleID 
ar_data <- rawdat%>%
  full_join(meta_dat, by = c("SampleID")) 

## 4. Select GLB
ar_data_GBL <- ar_data %>%
  filter(site=="GBL")
str(ar_data_GBL)

ar_data_GBL %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = X40.Conc, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

ar_data_GBL$arncalc <- c(ar_data_GBL$X40.Conc/ar_data_GBL$X28.Conc)
# compare with 1 and 2 pt standard curves

ar_data_GBL %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

## 4. Calculate theoretical Ratio of Ar:N2 
ar_data_GBL$arnsat <- arsat(ar_data_GBL$temp_C,ar_data_GBL$pressure_Hg) / nsat(ar_data_GBL$temp_C,ar_data_GBL$pressure_Hg)

ar_data_GBL %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = arnsat, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

## check out theoretical Ar
ar_data_GBL$arsat <- arsat(ar_data_GBL$temp_C,ar_data_GBL$pressure_Hg) 

ar_data_GBL %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

## check out theoretical N2
ar_data_GBL$nsat <- nsat(ar_data_GBL$temp_C, ar_data_GBL$pressure_Hg) 

ar_data_GBL %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = nsat, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)


GLB_data_post<- ar_data_GBL%>%filter(sample_type=='POST')

## Corrections for background Ar:N2
GLB_data_post$arn_corr <- GLB_data_post$arncalc - GLB_data_post$arnsat
GLB_data_post$ar_corr <- GLB_data_post$X40.Conc - GLB_data_post$arsat
GLB_data_post$n_corr <- GLB_data_post$X28.Conc - GLB_data_post$nsat

GLB_data_post %>%
  ggplot(aes(x = station, y = arn_corr, shape=sample_rep, color=as.factor(site))) +
  geom_line()+ geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

GLB_data_postq <- GLB_data_post %>%
  filter(!(trial == 2 & station_no == 1 & sample_rep == "C"))

GLB_data_postq1 <- GLB_data_postq %>%
  filter(!(trial == 2 & station_no == 0 & sample_rep == "C"))

GLB_data_postq2 <- GLB_data_postq1 %>%
  filter(!(trial == 6))

GLB_data_postq3 <- GLB_data_postq2 %>%
  filter(!(trial == 1 & station_no == 0 & sample_rep == "C"))

GLB_data_postq4 <- GLB_data_postq3 %>%
  filter(!(trial == 1 & station_no == 0 & sample_rep == "B"))


GLB_data_postq4 %>%
  ggplot(aes(x = station, y = arn_corr, shape=sample_rep, color=as.factor(site))) +
  geom_line()+ geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

### 5. Normalize for first well mixed station 
GLB_data_post_proc <- GLB_data_postq4 %>%
  filter(station_no>0) %>%
  group_by(trial, site, date) %>%
  mutate(
    norm_arncalc = c(arn_corr / max(arn_corr[station_no == 1])),
    norm_ar_calc = c(ar_corr / max(ar_corr[station_no == 1])), 
    norm_n_calc = c(n_corr / max(n_corr[station_no == 1]))
  ) %>%
  ungroup()

GLB_data_post_proc %>%
  ggplot(aes(x = station, y = norm_arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) 

## might not have worked station 6 looks bad 


###
###

### 5. Normalize for first well mixed station 
GLB_data_post_proc_2 <- GLB_data_postq4 %>%
  filter(trial == 2) %>%
  filter(station_no > 1) %>%
  group_by(trial, site, date) %>%
  mutate(
    norm_arncalc = c(arn_corr / max(arn_corr[station_no == 2])),
    norm_ar_calc = c(ar_corr / max(ar_corr[station_no == 2])), 
    norm_n_calc = c(n_corr / max(n_corr[station_no == 2]))
  ) %>%
  ungroup()

GLB_data_post_proc_2 %>%
  ggplot(aes(x = station, y = norm_arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) 



### 5. Normalize for first well mixed station 

GLB_data_t6 <- GLB_data_post %>%
  filter((trial == 6 & station_no > 2))

GLB_data_t6 %>%
  ggplot(aes(x = station, y = arn_corr, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) 

GLB_data_post_proc_6 <- GLB_data_t6 %>%
  group_by(trial, site, date) %>%
  mutate(
    norm_arncalc = c(arn_corr / max(arn_corr[station_no == 3])),
    norm_ar_calc = c(ar_corr / max(ar_corr[station_no == 3])), 
    norm_n_calc = c(n_corr / max(n_corr[station_no == 3]))
  ) %>%
  ungroup()

GLB_data_post_proc_6 %>%
  ggplot(aes(x = station, y = norm_arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) 

################
## Stan model
sink("NEON_NP_Ar.stan")

cat("data {
  int<lower = 1> N; // Number of observations
  vector[N] dist;   // Distances for the single trial
  vector[N] Ar;     // Argon measurements (log-transformed if needed)
  real Q;           // Discharge for the trial
  real V;           // Velocity for the trial
  real temp;        // Temperature for the trial (used for converting to K600)
}

parameters {
  real logK600;           // log-transformed K600 parameter
  real intercept;         // Intercept parameter, close to 0
  real<lower = 0> sigma;  // Standard deviation of the residuals
  real a;                 // Coefficient for the relationship with log(Q)
  real b;                 // Slope for the relationship with log(Q)
}

transformed parameters {
  real<lower=0> Kd;       // Decay rate per distance
  real<lower=0> KAr;      // Argon exchange coefficient

  // Calculate KAr based on temperature and logK600
  KAr = exp(logK600) / ((600 / (1759.7 - 117.37 * temp + 3.6959 * temp^2 - 0.046527 * temp^3))^-0.5);

  // Calculate Kd as KAr divided by velocity
  Kd = KAr / V;  // KAr (1/day) / V (m/day) = Kd (1/m)
}

model {
  // Likelihood for the Argon measurements
  for (i in 1:N) {
    Ar[i] ~ normal(intercept + -Kd * dist[i], sigma);
  }

  // Prior for the logK600 parameter
  logK600 ~ normal(a + b * log(Q), 1); // Use a prior that reflects the relationship with log(Q)

  // Priors for the other parameters
  a ~ normal(0, 10);
  b ~ normal(0, 1);
  sigma ~ normal(0, 0.2);     // Prior on sigma
  intercept ~ normal(0, 0.1); // Prior on intercept
}
"
,fill=TRUE)
sink()



#@@@@@###

### 6. Predict Ar based on exponential decay for downstream stations 
###################################
## Stan model where experiments are no longer pooled :
sink("NEON_Ar_no_pool.stan")

cat("
data {
  int<lower = 1> N;               // Number of observations for the current experiment
  int<lower = 1> nexpt;           // Current experiment index (should always be 1 in this setup)
  vector[N] dist;                 // Distance observations
  vector[N] Ar;                   // Observed Ar (e.g., log-transformed)
  real Q;                         // Discharge for the current experiment
  real V;                         // Velocity for the current experiment
  real temp;                      // Stream temperature (for K600 calculation)
  real w;                         // Other data related to the experiment
}

parameters {
  real logK600;                   // Log-transformed K600 for the current experiment
  real intercept;                 // Intercept for the Ar ~ distance relationship
  real<lower = 0> sigma;          // Standard deviation for likelihood
  real a;                         // Regression parameter for log(Q)
  real b;                         // Regression slope for log(Q)
  real<lower = 0> sigma_expt;     // Standard deviation for experiment-specific logK600
}

transformed parameters {
  real<lower = 0> Kd;             // Kd (1/m) for the current experiment
  real<lower = 0> KAr;            // KAr (1/day) for the current experiment

  KAr = exp(logK600) / ((600 / (1759.7 - (117.37 * temp) + (3.6959 * temp^2) - (0.046527 * temp^3)))^-0.5);
  Kd = KAr / V;                   // KAr (1/day) divided by V (m/day)
}

model {
  // Likelihood for Ar observations
  for (i in 1:N) {
    Ar[i] ~ normal(intercept + -Kd * dist[i], sigma);
  }

  // Priors
  logK600 ~ normal(a + b * log(Q), sigma_expt);
  a ~ normal(0, 10);
  b ~ normal(0, 1);
  sigma_expt ~ normal(0, 2);
  sigma ~ normal(0, 0.2);
  intercept ~ normal(0, 0.1);
}
    "
,fill=TRUE)
sink()



#################

## Format data for model:
GLB_data_post_norm <- GLB_data_post_proc %>%
  group_by(trial) %>%  
  mutate(dist = station - station[station_no == 1])

unique(GLB_data_post_norm$trial)

GW_data_post_1 <- GLB_data_post_norm %>%
  group_by(trial) %>%  
  filter(trial==1)  %>%
  mutate(trial = case_when(
    trial == 1 ~ 1,  
    TRUE ~ NA_real_ ))

GW_data_post_2 <- GLB_data_post_norm %>%
  group_by(trial) %>%  
  filter(trial==2)  %>%
  mutate(trial = case_when(
    trial == 2 ~ 1,  
    TRUE ~ NA_real_ ))


GW_data_post_4 <- GLB_data_post_norm %>%
  group_by(trial) %>%  
  filter(trial==4)  %>%
  mutate(trial = case_when(
    trial == 4 ~ 1,  
    TRUE ~ NA_real_ ))

GW_data_post_8 <- GLB_data_post_norm %>%
  group_by(trial) %>%  
  filter(trial==8)  %>%
  mutate(trial = case_when(
    trial == 8 ~ 1,  
    TRUE ~ NA_real_ ))


GLB_data_norm_2 <- GLB_data_post_proc_2 %>%
  group_by(trial) %>%
  mutate(dist = station - station[station_no == 2])

GB_data_post_2 <- GLB_data_norm_2 %>%
  group_by(trial) %>%
  filter(trial==2)  %>%
  mutate(trial = case_when(
    trial == 2 ~ 1,
    TRUE ~ NA_real_ ))

# GLB_data_post_proc_6


GLB_data_norm_6 <- GLB_data_post_proc_6 %>%
  group_by(trial) %>%
  mutate(dist = station - station[station_no == 3])

GB_data_post_6 <- GLB_data_norm_6 %>%
  group_by(trial) %>%
  filter(trial==6)  %>%
  mutate(trial = case_when(
    trial == 6 ~ 1,
    TRUE ~ NA_real_ ))


stan_data3 <- list(
  N = nrow(GB_data_post_6),
  nexpt = 1,
  dist = GB_data_post_6$dist,
  Ar = GB_data_post_6$norm_arncalc,
  Q = first(GB_data_post_6$Q_cms),
  V = first(GB_data_post_6$v),
  temp = first(GB_data_post_6$temp_C),
  w = first(GB_data_post_6$w)
)



# Run the model:
arfit <- stan(file = "NEON_Ar_no_pool.stan", data = stan_data3,
              iter = 5000, chains = 3,
              warmup = 2500, thin = 1,
              control = list(adapt_delta = 0.95))
##########
#@@@@@###

#################

fit_summary <- summary(arfit, probs=c(0.025,0.5,0.975))$summary %>% 
  {as_tibble(.) %>%
      mutate(var = rownames(summary(arfit)$summary))}

plot(arfit)




# Extract model parameters from output
stan_samples <- rstan::extract(arfit)
logK600_est <- apply(stan_samples$logK600, 1, mean)
Kd_est <- apply(stan_samples$Kd, 1, mean)
KAr_est <- apply(stan_samples$KAr, 1, mean)


stan_summary <- as.data.frame(summary(arfit)$summary)

# Create a dataframe for the Stan results
stan_results <- data.frame(
  trial = 1:length(logK600_est),  # assuming trial is indexed from 1 to nexpt
  logK600 = logK600_est,
  K600 = exp(logK600_est),
  Kd = Kd_est,
  KAr = KAr_est
)

# Summarize metadata by trial
trial_metadata <- GB_data_post_6 %>%
  group_by(trial) %>%
  summarize(
    date = first(date),
    site = first(site),
    Q_cms = first(Q_cms),
    temp_C = first(temp_C),
    w = first(w),
    v = first(v)
  )

# Merge the Stan results with metadata
merged_results <- trial_metadata %>%
  left_join(stan_results, by = "trial")

# save fits?
## Export model fit
##==================================
# export path
# getwd()
output_path <- paste0("")
output_path_sum <- paste0("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/model_data/K_Model_sum/")
output_path_fit <- paste0("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/model_data/K_Model_fits")

# extra label for pooling 
trial_name <- c("GBL_Trial_6")
# save model full output
# saveRDS(arfit, paste0(output_path_fit,"/",file_name,trial_name,"_K_estimates_",rundate,"_.rds"))
# write_csv(fit_summary, paste0(output_path_sum,"/",file_name,trial_name,"_K_estimates_",rundate,".csv"))
# write_csv(merged_results, paste0(output_path_sum,"/",file_name,trial_name,"_K_metadata_",rundate,".csv"))




##################
# Step 1: Extract model parameters from arfit
stan_samples <- rstan::extract(arfit)
logK600_est <- apply(stan_samples$logK600, 1, mean)  # Mean logK600 per trial
Kd_est <- apply(stan_samples$Kd, 1, mean)            # Mean Kd per trial
KAr_est <- apply(stan_samples$KAr, 1, mean)          # Mean KAr per trial
intercept_est <- mean(stan_samples$intercept)        # Mean intercept
sigma_est <- mean(stan_samples$sigma)                # Mean sigma

# Step 2: Organize the data by site and trial
unique_trials <- unique(GB_data_post_6[, c("site", "trial")])

# Function to predict Ar values based on the distance and other parameters
predict_Ar <- function(initial_Ar, Kd, dist_diff) {
  # Predict the downstream Ar concentration based on the model equation
  Ar_predicted <- initial_Ar * exp(-Kd * dist_diff)
  return(Ar_predicted)
}

# List to store results
results <- list()

# Loop through unique trials to process one site and trial at a time
for (i in 1:nrow(unique_trials)) {
  current_site <- unique_trials$site[i]
  current_trial <- unique_trials$trial[i]
  
  # Filter the data for the current site and trial
  trial_data <- GB_data_post_6 %>%
    filter(site == current_site, trial == current_trial) %>%
    arrange(.data$dist) # Explicitly reference column `dist`
  
  # Initialize predicted Ar values
  trial_data$Ar_predicted <- NA
  initial_Ar <- trial_data$norm_arncalc[1]  # Assuming the first station provides the initial Ar
  
  # Get estimated Kd for the current trial
  Kd <- Kd_est[current_trial]
  
  # Predict Ar values for each station sequentially
  trial_data$Ar_predicted[1] <- initial_Ar
  for (j in 2:nrow(trial_data)) {
    dist_diff <- trial_data$dist[j] - trial_data$dist[j - 1]  # Distance between stations
    # Predict Ar based on the previous station's Ar value
    trial_data$Ar_predicted[j] <- predict_Ar(trial_data$Ar_predicted[j - 1], Kd, dist_diff)
  }
  
  # Extract posterior predictive samples to calculate the 95% credible interval
  Ar_samples <- matrix(NA, nrow = nrow(stan_samples$Kd), ncol = nrow(trial_data))
  Ar_samples[, 1] <- initial_Ar
  for (j in 2:nrow(trial_data)) {
    dist_diff <- trial_data$dist[j] - trial_data$dist[j - 1]
    Ar_samples[, j] <- Ar_samples[, j - 1] * exp(-stan_samples$Kd * dist_diff)
  }
  
  # Calculate the 95% credible interval
  trial_data$Ar_lower <- apply(Ar_samples, 2, quantile, probs = 0.025) # Note column-based quantiles
  trial_data$Ar_upper <- apply(Ar_samples, 2, quantile, probs = 0.975)
  
  # Store the results for this trial
  results[[paste(current_site, current_trial, sep = "_")]] <- trial_data
}

######

# Step 4: Combine the results into a single dataframe
all_results <- do.call(rbind, results)

## exponential decay for observed plots 
all_results1 <- all_results %>%
  group_by(trial) %>%
  mutate(
    # Fit the exponential decay model: norm_arncalc = a * exp(-b * dist)
    exp_model = list(nls(norm_arncalc ~ a * exp(-b * dist), 
                         start = list(a = max(norm_arncalc), b = 0.01))),
    # Get fitted values from the model
    norm_arncalc_fitted = predict(exp_model[[1]])
  ) %>%
  ungroup()

# Calculate RMSE for each trial
rmse_results <- all_results1 %>%
  group_by(trial) %>%
  summarize(
    RMSE = round(sqrt(mean((norm_arncalc - Ar_predicted)^2, na.rm = TRUE)),3)
  )

all_results2 <- all_results1 %>%
  left_join(rmse_results, by = "trial")

# Step 5: Plot observed vs. predicted Argon concentrations with 95% credible interval
GB_plot <- ggplot(all_results2, aes(x = dist, y = norm_arncalc)) +
  geom_point(aes(color = "Observed")) +
  geom_line(aes(y = norm_arncalc_fitted, color = "Observed")) +
  geom_point(aes(y = Ar_predicted, color = "Predicted")) +
  geom_line(aes(y = Ar_predicted, color = "Predicted")) +
  geom_ribbon(aes(ymin = Ar_lower, ymax = Ar_upper), alpha = 0.3, fill = "#b81212") + 
  xlim(0,110) +
  facet_wrap(~ trial, scales = "free_x") +
  labs(title = "Glenbrook creek exp 6",
       x = "Distance (m)",
       y = "Normalized argon concentration") +
  theme_bw() +
  scale_color_manual(values = c("Observed" = "#1a4791", "Predicted" = "#b81212"))+
  # Add RMSE annotations
  geom_text(
    data = rmse_results,
    aes(x = Inf, y = Inf, label = paste("RMSE:", RMSE)),
    hjust = 1.1, vjust = 1.1,
    size = 3,
    color = "black",
    inherit.aes = FALSE
  )

GB_plot

## ggsave("/Users/kellyloria/Documents/UNR/Reaeration/GBL_neon_fit_T06.png", plot = GB_plot, width = 5, height = 3.75, units = "in")

summary(all_results$Ar_lower)
summary(all_results$Ar_upper)
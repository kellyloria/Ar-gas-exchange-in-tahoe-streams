## Merges clean MIMS data with injection meta data
## Re-worked from Hall and Madinger: Use of argon to measure gas exchange in turbulent mountain streams

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
  filter(!(trial == 6))

GLB_data_t6 <- GLB_data_postq %>%
  filter((trial == 6 & station_no > 3))

GLB_data_t61 <- GLB_data_t6 %>%
  mutate(station_no = case_when(
    station_no == 4 ~ 1,
    station_no == 5 ~ 2,
    station_no == 6 ~ 3,
    station_no == 7 ~ 4,
    TRUE ~ NA_real_  # To handle any unexpected values
  ))

ar_data_post_q1 <- rbind(GLB_data_postq1, GLB_data_t61)

ar_data_post_q1 %>%
  ggplot(aes(x = station, y = arn_corr, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) 

### 5. Normalize for first well mixed station 
GLB_data_post_q2 <- ar_data_post_q1 %>%
  filter(station_no>0) %>%
  group_by(trial, site, date) %>%
  mutate(
    norm_arncalc = c(arn_corr / mean(arn_corr[station_no == 1])),
    norm_ar_calc = c(ar_corr / mean(ar_corr[station_no == 1])), 
    norm_n_calc = c(n_corr / mean(n_corr[station_no == 1]))
  ) %>%
  ungroup()

GLB_data_post_q2 %>%
  ggplot(aes(x = station, y = norm_arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) 

## might not have worked station 6 looks bad 


################
## Stan model
sink("NEON_multi_Ar.stan")

cat("
data {
  int<lower = 1> N;
  int<lower = 1> nexpt;
  int<lower = 1> exptID[N];
  vector[N] dist;
  vector[N] Ar; // need logAr
  vector[nexpt] Q;
  vector[nexpt] V;
  vector [nexpt] temp; // one value for each stream temp to convert to K600
  vector[nexpt] w; 
  
}



parameters {
  vector [nexpt] logK600;
  real intercept; //very close to 0
  real<lower = 0> sigma; //  standard deviation
  real a; // real <lower = 0> a;
  real b;
  real <lower = 0> sigma_expt;
  
  
}


transformed parameters{
  vector <lower=0> [nexpt] Kd;
  vector <lower=0> [nexpt] KAr;
  
  for (j in 1: nexpt){  // make a loop here // very these smitts 
    KAr[j] = exp(logK600[j]) /  ((600/(1759.7-(117.37*temp[j])+(3.6959*temp[j]^2)-(0.046527*temp[j]^3)))^-0.5); 
  }
  
  Kd= KAr ./ V;  //KAr (1/day) / V (m/day) = Kd (1/m)
  
  
}

model {
  for (i in 1:N){
    Ar[i] ~ normal(intercept + -Kd[exptID[i]]*dist[i], sigma); // likelihood 
    // if you know travel time you solve for travel time here instead of distance 
  }
  
  for (j in 1:nexpt){
    logK600[j]~normal( a + b*log(Q[j]) , sigma_expt);
  }
  
  a ~ normal(0,10);
  b ~ normal(0,1);
  sigma_expt ~ normal (0,2);
  sigma ~ normal (0,0.2);  // added prior on sigma
  intercept ~ normal (0, 0.1); // working with proportions so strong prior on intercept
}
    "
,fill=TRUE)
sink()


#################

GB_data_post_norm <- GLB_data_post_q2 %>%
  group_by(trial) %>%  
  mutate(dist = station - station[station_no == 1])


unique(GB_data_post_norm$trial)

# re-name the trial 1-5 
GB_data_post_norm <- GB_data_post_norm %>%
  mutate(trial = case_when(
    trial == 1 ~ 1,
    trial == 2 ~ 2,
    trial == 4 ~ 3,
    trial == 6 ~ 4,
    trial == 8 ~ 5,
    TRUE ~ NA_real_  # To handle any unexpected values
  ))
ar_data_post_q <- GB_data_post_norm

arstandata <- list(
  N = nrow(ar_data_post_q),  # total number of observations
  nexpt = length(unique(ar_data_post_q$trial)),  # number of experiments (trials)
  exptID = ar_data_post_q$trial,  # experiment IDs
  dist = ar_data_post_q$dist,  # distance of each station
  Ar = ar_data_post_q$norm_arncalc,  # normalized argon proportion
  Q = ar_data_post_q %>%
    group_by(trial) %>%
    summarize(Q = first(Q_cms)) %>%
    pull(Q),  # discharge per experiment
  V = ar_data_post_q %>%
    group_by(trial) %>%
    summarize(V = first(v)) %>%
    pull(V),  # water velocity per experiment
  temp = ar_data_post_q %>%
    group_by(trial) %>%
    summarize(temp = first(temp_C)) %>%
    pull(temp),  # temperature per experiment
  w = ar_data_post_q %>%
    group_by(trial) %>%
    summarize(w = first(w)) %>%
    pull(w) 
)
# Check the list structure
str(arstandata)

# Run the model:
arfit <- stan(file = "NEON_multi_Ar.stan", data = arstandata,
              iter = 5000, chains = 3,
              warmup = 2500, thin = 1,
              control = list(adapt_delta = 0.95))

fit_summary <- summary(arfit, probs=c(0.025,0.5,0.975))$summary %>% 
  {as_tibble(.) %>%
      mutate(var = rownames(summary(arfit)$summary))}

plot(arfit)
plot(arfit, pars = c("Kd[1]", "Kd[2]", "Kd[3]", "Kd[4]", "Kd[5]"))

plot(arfit, pars = c("KAr[1]", "KAr[2]", "KAr[3]", "KAr[4]", "KAr[5]"))

# Extract model parameters from output
stan_samples <- rstan::extract(arfit)
logK600_est <- apply(stan_samples$logK600, 2, mean)
Kd_est <- apply(stan_samples$Kd, 2, mean)
KAr_est <- apply(stan_samples$KAr, 2, mean)



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
trial_metadata <- GB_data_post_norm %>%
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
# trial 6 is now trial 4 

# save fits?
## Export model fit
##==================================
# export path
# getwd()
output_path <- paste0("")
output_path_sum <- paste0("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/model_data/K_Model_sum/")
output_path_fit <- paste0("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/model_data/K_Model_fits")

# save model full output
# saveRDS(arfit, paste0(output_path_fit,"/",file_name,"_K_estimates_",rundate,"_.rds"))
# write_csv(fit_summary, paste0(output_path_sum,"/",file_name,"_K_estimates_",rundate,".csv"))
# write_csv(merged_results, paste0(output_path_sum,"/",file_name,"_K_metadata_",rundate,".csv"))


####
# end of temp edit 
####

##################
# Step 1: Extract model parameters from arfit
stan_samples <- rstan::extract(arfit)
logK600_est <- apply(stan_samples$logK600, 2, mean)  # Mean logK600 per trial
Kd_est <- apply(stan_samples$Kd, 2, mean)            # Mean Kd per trial
KAr_est <- apply(stan_samples$KAr, 2, mean)          # Mean KAr per trial
intercept_est <- mean(stan_samples$intercept)        # Mean intercept
sigma_est <- mean(stan_samples$sigma)                # Mean sigma

# Step 2: Organize the data by site and trial
unique_trials <- unique(ar_data_post_q[, c("site", "trial")])

# Function to predict Ar values based on the distance and other parameters
predict_Ar <- function(initial_Ar, Kd, dist_diff) {
  # Predict the downstream Ar concentration based on the model equation
  Ar_predicted <- initial_Ar * exp(-Kd * dist_diff)
  return(Ar_predicted)
}

# Step 3: Generate predictive Ar values iteratively for each station in each trial
results <- list()
for (i in 1:nrow(unique_trials)) {
  current_site <- unique_trials$site[i]
  current_trial <- unique_trials$trial[i]
  
  # Filter the data for the current site and trial
  trial_data <- ar_data_post_q %>%
    filter(site == current_site, trial == current_trial) %>%
    arrange(dist)
  
  # Initialize predicted Ar values
  trial_data$Ar_predicted <- NA
  initial_Ar <- trial_data$norm_arncalc[1]  # Assuming the first station provides the initial Ar
  
  # Get estimated Kd for the current trial
  Kd <- Kd_est[current_trial]
  
  # Predict Ar values for each station sequentially
  trial_data$Ar_predicted[1] <- initial_Ar
  for (j in 2:nrow(trial_data)) {
    dist_diff <- trial_data$dist[j] - trial_data$dist[j-1]  # Distance between stations
    # Predict Ar based on the previous station's Ar value
    trial_data$Ar_predicted[j] <- predict_Ar(trial_data$Ar_predicted[j-1], Kd, dist_diff)
  }
  
  # Extract posterior predictive samples to calculate the 95% credible interval
  Ar_samples <- matrix(NA, nrow = length(stan_samples$Kd[, current_trial]), ncol = nrow(trial_data))
  Ar_samples[, 1] <- initial_Ar
  for (j in 2:nrow(trial_data)) {
    dist_diff <- trial_data$dist[j] - trial_data$dist[j-1]
    Ar_samples[, j] <- Ar_samples[, j-1] * exp(-stan_samples$Kd[, current_trial] * dist_diff)
  }
  
  # Calculate the 95% credible interval
  trial_data$Ar_lower <- apply(Ar_samples, 2, quantile, probs = 0.025)
  trial_data$Ar_upper <- apply(Ar_samples, 2, quantile, probs = 0.975)
  
  # Store the results for this trial
  results[[paste(current_site, current_trial, sep = "_")]] <- trial_data
}

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
  xlim(0,100) +
  facet_wrap(~ trial, scales = "free_x") +
  labs(title = "Glenbrook creek",
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

## ggsave("/Users/kellyloria/Documents/UNR/Reaeration/GBL_neon_fit.png", plot = GB_plot, width = 8, height = 5.8, units = "in")


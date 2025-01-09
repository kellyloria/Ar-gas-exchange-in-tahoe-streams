##==================================
## Models for the 5 trials at BWL:
##==================================
## Merges clean MIMS data with injection meta data
##  Re-worked from Hall and Madinger: Use of argon to measure gas exchange in turbulent mountain streams
##  "NP" refers to not pooled where each trial at the Reach BWL was modeled independently

packages <- c("rstan", "shinystan", "ggplot2", "readr", "tidyr")
lapply(packages, library, character.only = TRUE)

###
## Daily file info for saving and model comparison:
set.seed(2021)
rundate <- format(Sys.Date(), "%y%m%d")
file_name <- "ArN_NEON_dist_Kt_BWL_NP"

source("./AR_code_repo/mims_gas_functions_wHeKr.R")
###

##==================================
## 1. Read in processed MIMS data for Ar 
files <- list.files(paste("./MIMS_dat/processed_dat/Merged_processed_dat",sep=""), full.names = T)
rawdat <-  do.call(rbind, lapply
                   (files, read.csv, as.is=T, header = T))
str(rawdat)

## 2. Read in meta data for each trial location, time,  station distance etc.
meta_dat <- read.csv("./MIMS_dat/processed_dat/MIMS_SampleLog_24.csv") %>%
  mutate(SampleID=as.character(sampleID))
str(meta_dat)

## Fix stream velocity to be in m^2 per day
meta_dat$v <- c(meta_dat$v * 86400)

## 3. Merge by SampleID 
ar_data <- rawdat%>%
  full_join(meta_dat, by = c("SampleID")) 

##==================================
## 4. Select BWL
ar_data_BWL <- ar_data %>%
  filter(site=="BWL" & X40.Conc<2) # filter some outliers for this reach. 
str(ar_data_BWL)

## Visualize Ar at different trials
ar_data_BWL %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = X40.Conc, shape=sample_rep, color=as.factor(site))) +
  geom_line() + geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) # station 3 might be the better well mixed zone

## 4. Calculate observed Ratio of Ar:N2 
ar_data_BWL$arncalc <- c(ar_data_BWL$X40.Conc/ar_data_BWL$X28.Conc)

## 5. Calculate theoretical Ratio of Ar:N2 
ar_data_BWL$arnsat <- arsat(ar_data_BWL$temp_C,ar_data_BWL$pressure_Hg) / nsat(ar_data_BWL$temp_C,ar_data_BWL$pressure_Hg)

ar_data_BWL %>%
  filter(sample_type=="POST" & station_no >2)%>%
  ggplot(aes(x = station, y = arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) +   geom_line()+ theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

ar_data_BWL %>%
  filter(sample_type=="POST" & station_no >2)%>%
  ggplot(aes(x = station, y = arnsat, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) +  geom_line()+ theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

## Calculate theoretical Ar 
ar_data_BWL$arsat <- arsat(ar_data_BWL$temp_C,ar_data_BWL$pressure_Hg) 
##  Calculate theoretical N2
ar_data_BWL$nsat <- nsat(ar_data_BWL$temp_C, ar_data_BWL$pressure_Hg) 

ar_data_BWL %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = nsat, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + geom_line()+ theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

## get some summary to get at % changes in Ar
ar_data_BWL_sum <- ar_data_BWL%>%
  filter(sample_type=="POST")%>%
  group_by(trial, site, date, sample_type) %>%
  summarise(
    arnsat_m = mean(arnsat, na.rm=T),
    nsat_m = mean(nsat, na.rm=T),
    arsat_m = mean(arsat, na.rm=T)
  )

BW_data_post<- ar_data_BWL%>%filter(sample_type=='POST')

## Corrections for background Ar:N2
BW_data_post$arn_corr <- BW_data_post$arncalc - BW_data_post$arnsat
BW_data_post$ar_corr <- BW_data_post$X40.Conc - BW_data_post$arsat
BW_data_post$n_corr <- BW_data_post$X28.Conc - BW_data_post$nsat

BW_data_post %>%
  filter(sample_type=="POST" & station_no>2)%>%
  ggplot(aes(x = station, y = arn_corr, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + geom_line()+ theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

##==================================
## Adjust trials for replicate noise and distance of first "well-mixed" station
# Possibly not enough AR at the first two Blackwood trials 
BW_data_post1<- BW_data_post %>% filter(station_no>2) 

BW_data_postq <- BW_data_post1 %>%
  filter(!(trial == 9 & station_no == 7 & sample_rep == "A"))

BW_data_postq1 <- BW_data_postq %>%
  filter(!(trial == 9 & station_no == 3 & sample_rep == "C"))

BW_data_postq2 <- BW_data_postq1 %>%
  filter(!(trial == 10 & station_no == 3 & sample_rep == "B"))

BW_data_postq3 <- BW_data_postq2 %>%
  filter(!(trial == 10 & station_no == 3 & sample_rep == "C"))

BW_data_postq4 <- BW_data_postq3 %>%
  filter(!(trial == 5))

BW_data_postq5 <- BW_data_postq4 %>%
  filter(!(trial == 3))

##==================================
### 6. Normalize for first well mixed station 
BW_data_post_norm1 <- BW_data_postq5 %>% ## NOT trusting group_by...
  group_by(trial, site, date) %>%
  mutate(
    norm_arncalc = c(arn_corr / mean(arn_corr[station_no == 3], na.rm = TRUE)),
    norm_ar_calc = c(ar_corr / mean(ar_corr[station_no == 3],na.rm = TRUE)), 
    norm_n_calc = c(n_corr / mean(n_corr[station_no == 3],na.rm = TRUE))
  ) %>%
  ungroup()

BW_data_post_norm1 %>%
  ggplot(aes(x = station, y = norm_arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_line() +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)


BW_data_trial3 <- BW_data_post %>%
  filter((trial == 3))

BW_data_trial3a <- BW_data_trial3 %>%
  filter(!(station_no == 0 & sample_rep == "A"))

BW_data_trial3b <- BW_data_trial3a %>%
  filter(!(station_no == 1))
BW_data_trial3c <- BW_data_trial3b %>%
  filter(!(station_no == 2))

### Quick check for 3 
BW_data_post_normtemp <- BW_data_trial3c %>% ## NOT trusting group_by...
  group_by(trial, site, date) %>%
  mutate(
    norm_arncalc = c(arn_corr / mean(arn_corr[station_no == 0], na.rm = TRUE)),
    norm_ar_calc = c(ar_corr / mean(ar_corr[station_no == 0],na.rm = TRUE)), 
    norm_n_calc = c(n_corr / mean(n_corr[station_no == 0],na.rm = TRUE))
  ) %>%
  ungroup()

BW_data_post_normtemp %>%
  ggplot(aes(x = station, y = norm_arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_line() +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

##==================================
### 7. Predict Ar based on exponential decay for downstream stations 
## model "NEON_Ar_no_pool.stan" developed in Aho et al. 2024: Gas exchange velocities (k600), gas 
##    exchange rates (K600), and hydraulic geometries for streams and rivers 
##    derived from the NEON Reaeration field and lab collection data product 
##==================================
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

##==================================
## Format data for model:
BW_data_post_norm1 <- BW_data_post_norm1 %>%
  group_by(trial) %>%  
  mutate(dist = station - station[station_no == 3])

unique(BW_data_post_norm1$trial)

## trial==7
BW_data_post_7 <- BW_data_post_norm1 %>%
  group_by(trial) %>%  
  filter(trial==7)  %>%
  mutate(trial = case_when(
    trial == 7 ~ 1,  
    TRUE ~ NA_real_ ))

## trial==9
BW_data_post_9 <- BW_data_post_norm1 %>%
  group_by(trial) %>%  
  filter(trial==9)  %>%
  mutate(trial = case_when(
    trial == 9 ~ 1,  
    TRUE ~ NA_real_ ))
## trial==10
BW_data_post_10 <- BW_data_post_norm1 %>%
  group_by(trial) %>%  
  filter(trial==10)  %>%
  mutate(trial = case_when(
    trial == 10 ~ 1,  
    TRUE ~ NA_real_ ))

# see if we can slavage trial==3
BW_data_post_normtemp <- BW_data_post_normtemp %>%
  group_by(trial) %>%  
  mutate(dist = station - station[station_no == 0])

BW_data_post_3 <- BW_data_post_normtemp %>%
  group_by(trial) %>%  
  filter(trial==3)  %>%
  mutate(trial = case_when(
    trial == 3 ~ 1,  
    TRUE ~ NA_real_ ))

##==================================
## get data in list form for stan:

stan_data3 <- list(
  N = nrow(BW_data_post_3), # replace with trial no. 
  nexpt = 1,
  dist = BW_data_post_3$dist,
  Ar = BW_data_post_3$norm_arncalc,
  Q = first(BW_data_post_3$Q_cms),
  V = first(BW_data_post_3$v),
  temp = first(BW_data_post_3$temp_C),
  w = first(BW_data_post_3$w)
)

##==================================

## 8. Run the model:
# arfit <- stan(file = "NEON_Ar_no_pool.stan", data = stan_data3,
#               iter = 5000, chains = 3,
#               warmup = 2500, thin = 1,
#               control = list(adapt_delta = 0.95))

##==================================

## Get model summary
fit_summary <- summary(arfit, probs=c(0.025,0.5,0.975))$summary %>% 
  {as_tibble(.) %>%
      mutate(var = rownames(summary(arfit)$summary))}

## Rough plot of fit:
plot(arfit)

## Extract model parameters from output
stan_samples <- rstan::extract(arfit)
logK600_est <- apply(stan_samples$logK600, 1, mean)
Kd_est <- apply(stan_samples$Kd, 1, mean)
KAr_est <- apply(stan_samples$KAr, 1, mean)

##  Create a dataframe for the Stan results
stan_results <- data.frame(
  trial = 1:length(logK600_est),  # trial is indexed from 1 to nexpt
  logK600 = logK600_est,
  K600 = exp(logK600_est),
  Kd = Kd_est,
  KAr = KAr_est
)

## Summarize metadata by trial
trial_metadata <- BW_data_post_3 %>%
  group_by(trial) %>%
  summarize(
    date = first(date),
    site = first(site),
    Q_cms = first(Q_cms),
    temp_C = first(temp_C),
    w = first(w),
    v = first(v)
  )

## Merge the Stan results with metadata
merged_results <- trial_metadata %>%
  left_join(stan_results, by = "trial")
##==================================
## save fits?
## Export model fit
##==================================
# export path
# getwd()
output_path <- paste0("")
output_path_sum <- paste0("./MIMS_dat/model_data/K_Model_sum/")
output_path_fit <- paste0("./MIMS_dat/model_data/K_Model_fits")

## extra label for pooling 
trial_name <- c("Trial3")

# save model full output
# saveRDS(arfit, paste0(output_path_fit,"/",file_name,trial_name,"_K_estimates_",rundate,"_.rds"))
# write_csv(fit_summary, paste0(output_path_sum,"/",file_name,trial_name,"_K_estimates_",rundate,".csv"))
# write_csv(merged_results, paste0(output_path_sum,"/",file_name,trial_name,"_K_metadata_",rundate,".csv"))


#################################################################
## Plot predicted verse observed Ar:N2 at stations 
##  to validate model fit
#################################################################
###
## Extract model parameters from arfit
stan_samples <- rstan::extract(arfit)
logK600_est <- apply(stan_samples$logK600, 1, mean)  # Mean logK600 per trial
Kd_est <- apply(stan_samples$Kd, 1, mean)            # Mean Kd per trial
KAr_est <- apply(stan_samples$KAr, 1, mean)          # Mean KAr per trial
intercept_est <- mean(stan_samples$intercept)        # Mean intercept
sigma_est <- mean(stan_samples$sigma)                # Mean sigma

## Organize the data by site and trial
unique_trials <- unique(BW_data_post_3[, c("site", "trial")])

## Function to predict Ar values based on the distance and other parameters
predict_Ar <- function(initial_Ar, Kd, dist_diff) {
  # Predict the downstream Ar concentration based on the model equation
  Ar_predicted <- initial_Ar * exp(-Kd * dist_diff)
  return(Ar_predicted)
}

## List to store results
results <- list()

## Loop through unique trials to process one site and trial at a time
for (i in 1:nrow(unique_trials)) {
  current_site <- unique_trials$site[i]
  current_trial <- unique_trials$trial[i]
  
  ## Filter the data for the current site and trial
  trial_data <- BW_data_post_3 %>%
    filter(site == current_site, trial == current_trial) %>%
    arrange(.data$dist) # Explicitly reference column `dist`
  
  ## Initialize predicted Ar values
  trial_data$Ar_predicted <- NA
  initial_Ar <- trial_data$norm_arncalc[1]  # the first station provides the initial Ar
  
  ## Get estimated Kd for the current trial
  Kd <- Kd_est[current_trial]
  
  ## Predict Ar values for each station sequentially
  trial_data$Ar_predicted[1] <- initial_Ar
  for (j in 2:nrow(trial_data)) {
    dist_diff <- trial_data$dist[j] - trial_data$dist[j - 1]  # Distance between stations
    # Predict Ar based on the previous station's Ar value
    trial_data$Ar_predicted[j] <- predict_Ar(trial_data$Ar_predicted[j - 1], Kd, dist_diff)
  }
  
  ## Extract posterior predictive samples to calculate the 95% credible interval
  Ar_samples <- matrix(NA, nrow = nrow(stan_samples$Kd), ncol = nrow(trial_data))
  Ar_samples[, 1] <- initial_Ar
  for (j in 2:nrow(trial_data)) {
    dist_diff <- trial_data$dist[j] - trial_data$dist[j - 1]
    Ar_samples[, j] <- Ar_samples[, j - 1] * exp(-stan_samples$Kd * dist_diff)
  }
  
  ## Calculate the 95% credible interval
  trial_data$Ar_lower <- apply(Ar_samples, 2, quantile, probs = 0.025) 
  trial_data$Ar_upper <- apply(Ar_samples, 2, quantile, probs = 0.975)
  
  ## Store the results for this trial
  results[[paste(current_site, current_trial, sep = "_")]] <- trial_data
}

######
## Combine the results into a single dataframe
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

# Plot observed vs. predicted Argon concentrations with 95% credible interval
BW_plot <- ggplot(all_results2, aes(x = dist, y = norm_arncalc)) +
  geom_point(aes(color = "Observed")) +
  geom_line(aes(y = norm_arncalc_fitted, color = "Observed")) +
  geom_point(aes(y = Ar_predicted, color = "Predicted")) +
  geom_line(aes(y = Ar_predicted, color = "Predicted")) +
  geom_ribbon(aes(ymin = Ar_lower, ymax = Ar_upper), alpha = 0.3, fill = "#b81212") + 
  xlim(0,110) +
  facet_wrap(~ trial, scales = "free_x") +
  labs(title = "Blackwood creek exp 3",
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

## ggsave("./BWL_neon_fit_T03_bad.png", plot = BW_plot, width = 5, height = 3.75, units = "in")
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
file_name <- "ArN_W1_dist_Kt_BWL"

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

## 4. Select BWL
ar_data_BWL <- ar_data %>%
  filter(site=="BWL" & X40.Conc<2)
str(ar_data_BWL)

ar_data_BWL$arncalc <- c(ar_data_BWL$X40.Conc/ar_data_BWL$X28.Conc)

## 4. Calculate theoretical Ratio of Ar:N2 
ar_data_BWL$arnsat <- arsat(ar_data_BWL$temp_C,ar_data_BWL$pressure_Hg) / nsat(ar_data_BWL$temp_C,ar_data_BWL$pressure_Hg)

ar_data_BWL %>%
  filter(sample_type=="POST" & station_no >2)%>%
  ggplot(aes(x = station, y = X40.Conc, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

ar_data_BWL %>%
  filter(sample_type=="POST" & station_no >2)%>%
  ggplot(aes(x = station, y = arnsat, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

## check out theoretical Ar
ar_data_BWL$arsat <- arsat(ar_data_BWL$temp_C,ar_data_BWL$pressure_Hg) 

ar_data_BWL %>%
  filter(sample_type=="POST" & station_no >2)%>%
  ggplot(aes(x = station, y = arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

## check out theoretical N2
ar_data_BWL$nsat <- nsat(ar_data_BWL$temp_C, ar_data_BWL$pressure_Hg) 

ar_data_BWL %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = nsat, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
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


BW_data_post1<- BW_data_post %>% filter(station_no>2) 

BW_data_postq <- BW_data_post1 %>%
  filter(!(trial == 9 & station_no == 7 & sample_rep == "A"))

BW_data_postq1 <- BW_data_postq %>%
  filter(!(trial == 9 & station_no == 3 & sample_rep == "C"))

BW_data_postq2 <- BW_data_postq1 %>%
  filter(!(trial == 10 & station_no == 3 & sample_rep == "B"))

BW_data_postq3 <- BW_data_postq2 %>%
  filter(!(trial == 10 & station_no == 3 & sample_rep == "C"))

### 5. Normalize for first well mixed station 
BW_data_post_norm1 <- BW_data_postq3 %>% ## NOT trusting group_by...
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

### 6. Predict Ar based on exponential decay for downstream stations 
###################################
## Stan model
sink("W1_dist_Kt.stan")

cat("
      
data {
  int<lower = 1> N;
  int<lower = 1> nexpt;
  int<lower = 1> exptID[N]; // indicator variable for expt
  vector[N] dist;
  vector[N] Ar; //proportion relative to upstreammost sampling station
  vector[nexpt] Q;
  vector[nexpt] V;
  vector [nexpt] temp;/// one value for each stream temp to convert to K600
  //vector[nexpt] w; 
  
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
  
  for (j in 1: nexpt){  // make a loop here
  KAr[j] = exp(logK600[j]) /  ((600/(1759.7-(117.37*temp[j])+(3.6959*temp[j]^2)-(0.046527*temp[j]^3)))^-0.5); 
  }
  
  Kd= KAr ./ V;  //KAr (1/day) / V (m/day) = Kd (1/m)
  
  
}

model {
  for (i in 1:N){
    Ar[i] ~ normal(intercept * exp(-Kd[exptID[i]]*dist[i]), sigma); // normal likelihood represnets MIMS error well
  }
  
  for (j in 1:nexpt){
    logK600[j]~normal( a + b*log(Q[j]) , sigma_expt); //next level model. Power model.
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



## Format data for model:
BW_data_post_norm1 <- BW_data_post_norm1 %>%
  group_by(trial) %>%  
  mutate(dist = station - station[station_no == 3])

unique(BW_data_post_norm1$trial)

BW_data_post_norm1 <- BW_data_post_norm1 %>%
  mutate(trial = case_when(
    trial == 3 ~ 1,
    trial == 5 ~ 2,
    trial == 7 ~ 3,
    trial == 9 ~ 4,
    trial == 10 ~ 5,
    TRUE ~ NA_real_  # To handle any unexpected values
  ))

ar_data_post_q <- BW_data_post_norm1 
  #rbind(GB_data_post_norm, BW_data_post_norm1)

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
    pull(temp)  # temperature per experiment
)
# Check the list structure
str(arstandata)

# Run the model:
arfit <- stan(file = "W1_dist_Kt.stan", data = arstandata,
              iter = 5000, chains = 3,
              warmup = 2500, thin = 1,
              control = list(adapt_delta = 0.95))

fit_summary <- summary(arfit, probs=c(0.025,0.5,0.975))$summary %>% 
  {as_tibble(.) %>%
      mutate(var = rownames(summary(arfit)$summary))}


# 10 trials so Kd
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
trial_metadata <- BW_data_post_norm1 %>%
  group_by(trial) %>%
  summarize(
    date = first(date),
    site = first(site),
    Q_cms = first(Q_cms),
    temp_C = first(temp_C),
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

# save model full output
# saveRDS(arfit, paste0(output_path_fit,"/",file_name,"_K_estimates_",rundate,"_.rds"))
# write_csv(fit_summary, paste0(output_path_sum,"/",file_name,"_K_estimates_",rundate,".csv"))
# write_csv(merged_results, paste0(output_path_sum,"/",file_name,"_K_metadata_",rundate,".csv"))


####
# end of temp edit 
####

## Extra figures
hist(ar_data_post$norm_arncalc)


plot_n2AR<- ar_data_post %>%
  ggplot(aes(x = dist, y = log(norm_arncalc+1), shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) +  
  labs(x = "station distance (m)",
       y = 'Ar:N2 concnetrations normalized to station 1') +
  scale_color_manual(values = c("#0b2549", "#daa520")) 

# ggsave("/Users/kellyloria/Documents/UNR/Reaeration/ARN2_trials.png", plot = plot_n2AR, width = 10, height = 6, units = "in")


plot_AR<- ar_data_post %>%
  ggplot(aes(x = dist, y = norm_arncalc, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial) +  
  labs(x = "station distance (m)",
       y = 'Ar concnetrations normalized to station 1') +
  scale_color_manual(values = c("#0b2549", "#daa520")) 

# ggsave("/Users/kellyloria/Documents/UNR/Reaeration/ARN2_trials.png", plot = plot_AR, width = 10, height = 6, units = "in")



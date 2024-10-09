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
file_name <- "Ar_prelim_fit"

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
            
## 3. Merge by SampleID 
ar_data <- rawdat%>%
  full_join(meta_dat, by = c("SampleID"))

## Look at some values with distance from injection:
ar_data %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = X40.Conc, shape=sample_rep, color=as.factor(trial))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ site)


ar_data %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = X40.Conc, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

ar_data %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = N2.Ar.Conc, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

# Still trying to work out what "arncalc" should be for these trials. 

ar_data$arncalc <- c(ar_data$X40.Conc)
# compare with 1 and 2 pt standard curves

## 4. Calculate theoretical Ratio of Ar:N2 
ar_data$arnsat <- arsat(ar_data$Temp,ar_data$Pressure) / nsat(ar_data$Temp,ar_data$Pressure)

ar_data %>%
  filter(sample_type=="POST")%>%
  ggplot(aes(x = station, y = arnsat, shape=sample_rep, color=as.factor(site))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ trial)

ar_data_pre<- ar_data%>%filter(sample_type=='PRE ')
ar_data_post<- ar_data%>%filter(sample_type=='POST')

## some stuff around SPC corrections, possibly for GW?
# ## %. Now for calculations where we subtract background and correct for conductivity
# ar_data_post$arcorr <- ar_data_post$arncalc - ar_data_post$arnsat
# 
# ar_data_post %>%
#   filter(sample_type=="POST")%>%
#   ggplot(aes(x = station, y = arcorr, shape=sample_rep, color=as.factor(site))) +
#   geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
#   facet_wrap(~ trial)

# ## Calc mean of all pre cond by station and by site
# ardataprecond <- ar_data_pre %>%  
#   group_by(trial, Station) %>%  # not really sure why station cor instead of station ?
#   summarise(precond=mean(SPC..uScm., na.rm=T)) # cond = SPC
# 
# ardatapost<-merge(ar_data_post,ardataprecond)
# 
# ardatapost$condcor<- ardatapost$SPC..uScm. - ardatapost$precond
# hist(ardatapost$condcor)
# 
# ardatapost %>% 
#   ggplot(aes(x = Station, y = condcor , color = sample_type, shape=sample_rep)) +
#   geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
#   facet_wrap(~ site)
# ar_data_post$arcond<- ar_data_post$arcorr / ar_data_post$condcor
#hist(ar_data_post$arcond)

### 5. Normalize for first well mixed station 
ar_data_post <- ar_data_post %>%
  group_by(trial, site, date) %>%
  mutate(norm_arncalc = arncalc / arncalc[station_no == 1]) %>%
  ungroup()

hist(ar_data_post$norm_arncalc)

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
##   filter data to stations greater than 1 and correct distance 
ar_data_post <- ar_data_post %>%
  arrange(trial, station) %>%
  filter(station > 0) %>%
  group_by(trial) %>%  
  mutate(dist = station - station[station_no == 1]) %>% 
  ungroup() 

arstandata <- list(
  N = nrow(ar_data_post),  # total number of observations
  nexpt = length(unique(ar_data_post$trial)),  # number of experiments (trials)
  exptID = ar_data_post$trial,  # experiment IDs
  dist = ar_data_post$dist,  # distance of each station
  Ar = ar_data_post$norm_arncalc,  # normalized argon proportion
  Q = ar_data_post %>%
    group_by(trial) %>%
    summarize(Q = first(Q_cms)) %>%
    pull(Q),  # discharge per experiment
  V = ar_data_post %>%
    group_by(trial) %>%
    summarize(V = first(v)) %>%
    pull(V),  # water velocity per experiment
  temp = ar_data_post %>%
    group_by(trial) %>%
    summarize(temp = first(temp_C)) %>%
    pull(temp)  # temperature per experiment
)
# Check the list structure
str(arstandata)

# ## Run the model:
arfit <- stan(file = "W1_dist_Kt.stan", data = arstandata,
              iter = 3000, chains = 3,
              warmup = 1500, thin = 1)

fit_summary <- summary(arfit, probs=c(0.025,0.5,0.975))$summary %>% 
  {as_tibble(.) %>%
      mutate(var = rownames(summary(arfit)$summary))}

plot(arfit)

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
  Kd = Kd_est,
  KAr = KAr_est
)

# Summarize metadata by trial
trial_metadata <- ar_data_post %>%
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


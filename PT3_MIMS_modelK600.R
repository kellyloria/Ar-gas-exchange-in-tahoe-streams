# ROUGH DRAFT ** script to try to run gas exchange model:
#' @description stan model W1_dist_Kt.stan
#' @param will be listed in ...
#' 
#' @return Returns .r file for running metabolism model 
#' @export 

##==============================================================================
## Notes: 
## Created  09/03/2024 by KAL
#===============================================================================

rm(list=ls())
getwd()

set.seed(2021)
# stan settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rundate <- format(Sys.Date(), "%y%m%d") # usefull fo iterating model runs

## KAL's temporary path reminders: 
## setwd("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat")

# load packages
library(tidyverse)
library(rstan)
library(loo)
library(patchwork)
library(lubridate)
## Example from lake metabolism 
# source("/Users/kellyloria/Desktop/Castle2018/LittoralAlder3M/stan_utility.R")
# data <- read_rdump(paste("./ModelInputs/F/",lake,"_",min(year),"_",max(year),"_5ms_fourthlake_Offset_sonde_list.R",sep=""))


# read data
# data needs to be in list that runs:
  # data {
  #   int<lower = 1> N;
  #   int<lower = 1> nexpt; # ? I think this the no. / count of experiments in the data array 
  #   int<lower = 1> exptID[N]; # experiments id as a number 
  #   vector[N] dist; # ? hmmm  for each Ar sample? 
  #   vector[N] Ar; # proportion relative to upstreammost sampling station
  #   vector[nexpt] Q; # ? probably mean for the entire reach
  #   vector[nexpt] V; # ? probably mean for the entire reach
  #   vector[nexpt] temp; # mean for entire reach 
  #   vector[nexpt] w; # mean for entire reach 
  # }

model_data <- readRDS("./model_data/meta_datq.rds") %>%
  mutate(date = as.Date((date), format ="%m/%d/%y"))

# Filter data for one experiment:
filtered_data <- model_data %>%
  filter(sample_type == "POST", 
         site == "GBL" , 
         date == as.Date("2024-07-10")) 

## need to fix the meta data workflow for Q, V, and W so just plugging them in here 
experiment_summary <- model_data %>%
  group_by(sample_type,site, date) %>%
  summarize(
    Q = c(0.016), #cms
    V = c(0.126), # 
    temp = c(13.1), #mean(temp..C.),
    w = c(1.77)
  )

# something is broken here with reading the the Q,V,temp,w in the correct array to match the unique "nexpt" length
# Prepare the list structure for Stan
data_list <- list(
  N = nrow(model_data),
  nexpt = n_distinct(model_data$site), # need to fix for and create label for exp in .rds 
  exptID = as.numeric(as.factor(model_data$site)), # Convert to numeric IDs
  dist = filtered_data$distance, # Distance vector
  Ar = log(filtered_data$X40/filtered_data$X28), # Proportion relative to upstream-most station
  Q = unique(experiment_summary$Q), # Mean flow by experiment
  V = unique(experiment_summary$V), # Mean velocity by experiment
  temp = unique(experiment_summary$temp) # Mean temp by experiment
  #w = unique(experiment_summary$w) # Mean width by experiment
)

data_list <- na.omit(data_list) 

# call the stan based model 
model <- "W1_dist_Kt.stan" #Steele 2 param inhibition
model_path <- paste0("./stan/",model) # edit for your file path

# complie model:
compiled_model <- stan_model(file = model_path)

# fit stan model to data:
fit <- sampling( 
  object = compiled_model,
  data = data_list,
  iter = 2000,        # Number of iterations
  chains = 3,         # Number of chains
  warmup = 1000,      # Number of warm-up iterations (1/2 iterations)
  thin = 1            # Thinning interval
)

print(fit, pars = c("logK600", "intercept", "sigma", "a", "b", "sigma_expt"))

stan_trace(fit, pars = c("logK600", "a", "b"))
stan_dens(fit, pars = c("logK600", "a", "b"))



##==================================
##==================================
############################
## *¡
## STOP .. 
## this is old code from metabolism stuff that might be nice for data export and visualization once edited
###########################
##==================================
##==================================



##==================================
## Assess model fit
##==================================

fit_summary <- summary(stanfit, probs=c(0.025,0.5,0.975))$summary %>% 
  {as_tibble(.) %>%
      mutate(var = rownames(summary(stanfit)$summary))}

# specific summaries
check_n_eff(stanfit)
check_rhat(stanfit)
check_div(stanfit)
check_treedepth(stanfit,max_treedepth)
check_energy(stanfit)

##==================================
## Export model fit
##==================================
# export path
output_path <- paste0("./ModelOutput/")
output_path_sum <- paste0("./ModelOutput/summary/")
output_path_fit <- paste0("./ModelOutput/fits/")


# save model full output
saveRDS(stanfit, paste0(output_path_fit,"/",lake,"_5ms_fourthlake_Offset_fit_",rundate,"_.rds"))

fit_clean <- fit_summary %>%
  dplyr:: rename(lower = '2.5%', middle = `50%`,upper = '97.5%')  %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         index = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         day = ifelse(name %in% c("GPP","ER","NEP","AIR","Flux","GPP_m2","ER_m2","NEP_m2","b0","r","i0","b"), 
                      index, 
                      data$map_days[index])) %>%
  select(name, index, day, middle,lower,upper)


##==================================
## Rejoin output with meta dat
##==================================
# * CHECK THIS FILE NAME 

sonde_data <- read_csv(paste("./ModelInputMeta/F/","sonde_dat_5ms_fourthlake_Offset_",lake,"_",min(year),"_",max(year),".csv",sep=""))


out <- fit_clean %>%
  filter(name %in% c("GPP","ER","NEP")) %>%
  dplyr::rename(unique_day = day) %>% 
  left_join(sonde_data %>% select(unique_day,yday,year) %>% distinct()) %>% 
  full_join(sonde_data %>% expand(year,yday,name=c("GPP","ER","NEP"))) %>% 
  mutate(middle = ifelse(name=="ER",-middle,middle),
         lower = ifelse(name=="ER",-lower,lower),
         upper = ifelse(name=="ER",-upper,upper),
         name = factor(name, levels=c("GPP","NEP","ER")))

out2 <- fit_clean %>%
  filter(name %in% c("GPP_m2","ER_m2","NEP_m2"))%>%
  dplyr::rename(unique_day = day) %>% 
  left_join(sonde_data %>% select(unique_day,yday,year) %>% distinct()) %>% 
  full_join(sonde_data %>% expand(year,yday,name=c("GPP_m2","ER_m2","NEP_m2"))) %>% 
  mutate(middle = ifelse(name=="ER_m2",-middle,middle),
         lower = ifelse(name=="ER_m2",-lower,lower),
         upper = ifelse(name=="ER_m2",-upper,upper),
         name = factor(name, levels=c("GPP_m2","NEP_m2","ER_m2")))

out3 <- rbind(out,out2)

out4 <- out3 %>%
  filter(year==2023)

##==================================
## Export model results
##==================================


write_csv(out4, paste0(output_path,"/",lake,"_","_daily_5ms_fourthlake_Offset_",rundate,".csv"))
write_csv(fit_clean, paste0(output_path_sum,"/",lake,"_","_summary_5ms_fourthlake_Offset_",rundate,".csv"))


##==================================
## Plot model results
##==================================
#plot primary parameters

figure_path <- paste0("./ModelOutput/figures/")

p1 <- fit_clean %>%  
  filter(name=="b0" | name == "r" | name == "i0") %>% 
  dplyr::rename(unique_day = index) %>% 
  left_join(sonde_data %>% filter(year==2023)%>%select(unique_day,yday,year) %>% distinct()) %>%  drop_na(year)%>%
  ggplot(aes(x=yday,y=middle,color=factor(year))) + 
  geom_point(size=0.5) +
  geom_line(alpha=0.5) +
  facet_wrap(vars(name, year),ncol=3,scales="free_y") +
  theme_bw() +
  labs(y="Mean Estimated Value",color="year",x="Day of Year") +
  coord_cartesian(ylim = c(0, 0.6)) 

p1

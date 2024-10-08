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
proc_dat <- read.csv("2024_KAL_prelimdat_processed_09.csv")%>%
  dplyr::select("Samp", "SampleID", "Pressure", "N2.Ar", 
                "WatDens", "O2Sat", "N2Sat", "ArSat",
                "O2.ArSat", "N2.ArSat", "X40", "X28", "X32")
str(proc_dat)
## 2. meta data for each trial location, time,  station distance etc.
meta_dat <- read.csv("kelly_prelim_sample_metadata.csv") %>%
  mutate(SampleID=as.character(SampleID))
str(meta_dat)
            
## 3. Merge by SampleID 
ar_data <- proc_dat%>%
  full_join(meta_dat, by = c("SampleID"))

## Look at some values with distance from injection:
ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = (ArSat) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 



## Look at some values with distance from injection:
ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = (X40.Conc) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = (N2.Ar.Conc) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 



ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = 1/(N2.Ar) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = 1/log(N2.Ar) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

## ?? HELP ??
# Still trying to work out what "arncalc" should be. 
#  Possibly it's ArSat but it may need to be calculated as the inverse?
#  link: "./SOP_literature/bg-15-3085-2018-supplement/Argon%20supplement/Arcode.html"

ar_data$arncalc <- c(ar_data$N2.Ar.Conc)
# compare with 1 and 2 pt standard curves


## 4. Correct for background saturation conditions and ground water:
## Ratio of Ar:N2 
ar_data$arnsat <- arsat(ar_data$temp..C.,ar_data$Pressure) / nsat(ar_data$temp..C.,ar_data$Pressure)

ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = arnsat , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

ar_data_pre<- ar_data%>%filter(sample_type=='PRE ')
ar_data_post<- ar_data%>%filter(sample_type=='POST')

## Now for calculations where we subtract background and correct for conductivity
ar_data_post$arcorr <- ar_data_post$arncalc - ar_data_post$arnsat

## Calc mean of all pre cond by station and by site
ardataprecond <- ar_data_pre %>%  
  group_by(trial, Station) %>%  # not really sure why station cor instead of station ?
  summarise(precond=mean(SPC..uScm., na.rm=T)) # cond = SPC

ardatapost<-merge(ar_data_post,ardataprecond)

ardatapost$condcor<- ardatapost$SPC..uScm. - ardatapost$precond
hist(ardatapost$condcor)

ardatapost %>% 
  ggplot(aes(x = Station, y = condcor , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ site)


ardatapost$arcond<- ardatapost$arcorr / ardatapost$condcor
hist(ardatapost$arcorr)
hist(ardatapost$arcond)

### 5. Normalize for first well mixed station 
ardatapost_0<-ardatapost[ardatapost$station.no.==1,]
ardata_0sum <- ardatapost_0 %>% 
  group_by(trial) %>% 
  summarise(arcond_0=mean(arcond), arn_enrich=mean(ArSat/arnsat))

#join with ardatapost
ardatapost<-merge(ardatapost,ardata_0sum)
ardatapost$arcondnorm<-ardatapost$arcond/ardatapost$arcond_0

ardatapost %>% 
  ggplot(aes(x = Station, y = arcondnorm , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  facet_wrap(~ site)

## ? HELP ## will I need to adjust distance relationships to 
## subtract first station meter (12) to the other station distances?

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
arstandata=list("Ar"= ardatapost$arcondnorm, "distar"=ardatapost$Station, "N"=length(ardatapost$Station), 
                "sites"=max(ardatapost$trial),  "sitear" = ardatapost$trial)

# ## Run the model:
# arfit <- stan(file = "W1_dist_Kt.stan", data = arstandata,
#               iter = 5000, chains = 3,
#               warmup = 2500, thin = 1)


print(arfit, pars=c("a","mu_a","sigma_ar", "k", "Ar0", "sigma_a" ))

arfitsum<- summary(arfit)$summary

asum<-arfitsum[1:8,c("2.5%", "50%", "97.5%")]

# Extract posterior summary for k (decline rate) and rescale:
ksum<-(arfitsum[12:19,c("2.5%", "50%", "97.5%")])*0.01
# Extract summary for ar_tilde (predicted normalized Ar concentration)
dat_tildesum<- summary(arfit, pars="ar_tilde", probs=0.5)$summary
ar_tildesum<-dat_tildesum[,4]
k_tildesum<-dat_tildesum[,3]

plot(ardatapost$arcondnorm, ar_tildesum, xlab="Measured normalized Ar concentration", ylab="Predicted normalized Ar concentration", pch=16, cex.lab=0.8, cex.axis=0.8)

plot(ardatapost$Station, ar_tildesum, xlab="Station", ylab="Predicted normalized Ar concentration", pch=16, cex.lab=0.8, cex.axis=0.8)

plot(ardatapost$Station, k_tildesum, xlab="Station", ylab="Predicted k", pch=16, cex.lab=0.8, cex.axis=0.8)

# save fits?
## Export model fit
##==================================
# export path
# getwd()
output_path <- paste0("")
output_path_sum <- paste0("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/model_data/K_Model_sum/")
output_path_fit <- paste0("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/model_data/K_Model_fits")

# save model full output
saveRDS(arfit, paste0(output_path_fit,"/",file_name,"_sample_run_JAK_",rundate,"_.rds"))


#####
## upload or change in meta data file for stream geomorphology and hydro
# Calculate derived variables such as K(per time) and k600

streamslope<-c(0.21)
# example dat from GBL 07/10/24
    # Q: 0.0243 cms 0.86 cfs 

# z <-(Q)/(w*v)
# v <- Q/(w*z)

Q <- c(0.0243) # discharge in cms
w <- c(1.7) # width in m 
z <- c(0.43) # depth in m 
v <- c(0.0332) # velocity in m per s

ardataposttemp<-ardatapost %>% 
  group_by(trial) %>% 
  summarise(temp=mean(temp..C., na.rm=T))

## NEED to double check some of this...
k<- ksum*Q*1440/w # 1440 is coming from ...?
# where k is gas exchange per unit time (per day)
k600<-K600fromAr(ardataposttemp, k)
# where k600 is gas exchange scaled at specific temperatures 


z<-(Q)/(w*v)

Kt<-v*ksum[,2]*1440


### save metadata merge file for list format:
# saveRDS(meta_datq, file = "./model_data/meta_datq.rds")


####
# end of temp edit 
####




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
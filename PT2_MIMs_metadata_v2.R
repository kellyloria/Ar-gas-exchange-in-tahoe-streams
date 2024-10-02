## Merges clean MIMS data with injection meta data
## Re-worked from Hall and Madinger: Use of argon to measure gas exchange in turbulent mountain streams

library(rstan)
library(shinystan)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)


setwd("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/")

source("/Users/kellyloria/Documents/UNR/Reaeration/AR_code_repo/mims_gas_functions_wHeKr.R")
###
proc_dat <- read.csv("2024_KAL_prelimdat_processed.csv")%>%
  dplyr::select("Samp", "SampleID", "Pressure", "N2.Ar", 
                "WatDens", "O2Sat", "N2Sat", "ArSat",
                "O2.ArSat", "N2.ArSat", "X40", "X28", "X32")

proc_dat$arncalc <- c(1/proc_dat$N2.Ar)

meta_dat <- read.csv("kelly_prelim_sample_metadata.csv") %>%
  mutate(Samp=as.integer(MIMs_label))
            
## 2. Merge with metadata 
ar_data <- proc_dat%>%
  full_join(meta_dat, by = c("Samp"))

## Look at some values with distance from injection:
ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = (ArSat) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = (N2.ArSat) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 


ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = 1/(N2.Ar) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

ar_data %>%
  filter(trial=="4")%>%
  ggplot(aes(x = Station, y = 1/(X28/X40) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 



# Station 3 D looking weird 
# Station 7 B also not great lets trim those 
### Fix this 

# ar_data <-ar_data %>%
#   filter(ArSat< 0.60)

# Still trying to work out what "arncalc" should be. Possibly it's ArSat but it may need to be calculated as the inverse?
# c(N2.Ar) coming from X28/X40 in other data sheet 
# link: "./SOP_literature/bg-15-3085-2018-supplement/Argon%20supplement/Arcode.html"


## Estimates K600 for KO2 at a given temperature. From Wanninkhof (1992).
K600fromO2<-function (temp, KO2) {
  ((600/(1800.6 - (120.1 * temp) + (3.7818 * temp^2) - (0.047608 * temp^3)))^-0.5) * KO2
}

## Estimates K600 for KAr at a given temperature. From Raymond et al  (2012).
K600fromAr<-function (temp, KAr) {
  ((600/(1799 - (106.96 * temp) + (2.797 * temp^2) - (0.0289 * temp^3)))^-0.5) * KAr
}

##
##
ar_data$arnsat <- arsat(ar_data$temp..C.,ar_data$Pressure) / nsat(ar_data$temp..C.,ar_data$Pressure)
ar_data_pre<- ar_data%>%filter(sample_type=='PRE ')
ar_data_post<- ar_data%>%filter(sample_type=='POST')

# Now for calculations where we subtract background and correct for conductivity
ar_data_post$arcorr<- ar_data_post$arncalc - ar_data_post$arnsat

#2.2.  Calc mean of all pre cond by station and by site
ardataprecond <- ar_data_pre %>% 
  group_by(trial, Station) %>%  # not really sure why station cor instead of station ?
  summarise(precond=mean(SPC..uScm., na.rm=T)) # cond = SPC


ardatapost<-merge(ar_data_post,ardataprecond)
hist(ardatapost$SPC..uScm.)

ardatapost$condcor<- ardatapost$SPC..uScm. - ardatapost$precond
hist(ardatapost$condcor)

ardatapost$arcond<- ardatapost$arcorr / ardatapost$condcor
hist(ardatapost$arcorr)
hist(ardatapost$arcond)


ardatapost_0<-ardatapost[ardatapost$Station==12,]
ardata_0sum <- ardatapost_0 %>% 
  group_by(trial) %>% 
  summarise(arcond_0=mean(arcond), arn_enrich=mean(ArSat/arnsat))

#join with ardatapost
ardatapost<-merge(ardatapost,ardata_0sum)
ardatapost$arcondnorm<-ardatapost$arcond/ardatapost$arcond_0



###################################
#Stan model
#Here is the Stan model that we used for the analysis. Details on priors etc are in the text.
sink("combineargonmodel_KL.stan")

cat("
    
    data {
    
    int <lower=1> N; //number of data points, ar
    int <lower=1> sites; //number of sites
    int <lower=0> sitear[N]; //stream number  for ar indexing vector

    
    vector [N] ar;//conductivity corrected Ar data
    vector [N] distar;//
    
    }
    
    parameters {
    
    vector <lower=0> [sites] a;
    real mu_a;   //take out for hyperprior info  // mean prior
    real<lower=0, upper=0.5> sigma_ar; // error ar

    vector[sites] k; // decline
    real <lower=0, upper=2> Ar0;
    
    //real d;
    //real b;
    real <lower=0> sigma_a; // mean prior
    }
    
    model { 
    
    //priors. 
    k ~ normal(0, 10);
  
    a~normal (mu_a,sigma_a); // mean prior

    Ar0 ~normal(1,0.05);

    mu_a~normal(1.35, 1); // mean prior
     sigma_a ~ normal(0, 2);
    
    //likelihood        
    for (i in 1:N) {
    ar[i] ~ normal( Ar0 * exp(-k[sitear[i]]*distar[i]*0.01) , sigma_ar); 
    }
    
    }

generated quantities {   //These estimate the posterior predicted

    vector [N] ar_tilde;

    for (n in 1:N) ar_tilde[n] = normal_rng( Ar0 * exp(-k[sitear[n]]*distar[n]*0.01) , sigma_ar);

}

    "
,fill=TRUE)
sink()


arstandata=list("ar"= ardatapost$arcondnorm, "distar"=ardatapost$Station, "N"=length(ardatapost$Station), 
                "sites"=max(ardatapost$trial),  "sitear" = ardatapost$trial)

arfit <- stan(file = "combineargonmodel_KL.stan", data = arstandata,
              iter = 5000, chains = 3,
              warmup = 2500, thin = 1)


print(arfit, pars=c("a","mu_a","sigma_ar", "k", "Ar0", "sigma_a" ))

arfitsum<- summary(arfit)$summary

asum<-arfitsum[1:8,c("2.5%", "50%", "97.5%")]

ksum<-(arfitsum[12:19,c("2.5%", "50%", "97.5%")])*0.01# the 0.01 here rescales the k estimate

#. Essentially we are using the model output to peredict a new set of data which we compare with our actual data.

dat_tildesum<- summary(arfit, pars="ar_tilde", probs=0.5)$summary
ar_tildesum<-dat_tildesum[,4]
k_tildesum<-dat_tildesum[,3]

plot(ardatapost$arcondnorm, ar_tildesum, xlab="Measured normalized Ar concentration", ylab="Predicted normalized Ar concentration", pch=16, cex.lab=0.8, cex.axis=0.8)

plot(ardatapost$Station, ar_tildesum, xlab="Station", ylab="Predicted normalized Ar concentration", pch=16, cex.lab=0.8, cex.axis=0.8)

plot(ardatapost$Station, k_tildesum, xlab="Station", ylab="Predicted k", pch=16, cex.lab=0.8, cex.axis=0.8)


#####
## upload or change in meta data file 

# Calculate derived variables such as K(per time) and k600
streamslope<-c(0.21)
Q <- c(0.097)*60
w <- c(0.7)

ardataposttemp<-ardatapost %>% 
  group_by(trial) %>% 
  summarise(temp=mean(temp..C., na.rm=T))

## NEED to double check some of this...
k<- ksum*Q*1440/w # 1440 is coming from ...?
k600<-K600fromAr(ardataposttemp$temp, k)
v<-c(6.3,5.2)

z<-(Q)/(w*v)

Kt<-v*ksum[,2]*1440

### visualize

ardatapost %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = Station, y = N2Sat , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 
  #scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)

ardatapost %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = Station, y = O2Sat , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

ardatapost %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = Station, y = ArSat , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

ardatapost %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = Station, y = N2.ArSat , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 


ardatapost %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = Station, y = arnsat , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 


ardatapost %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = Station, y = arcondnorm , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 


### save metadata merge file for list format:
# saveRDS(meta_datq, file = "./model_data/meta_datq.rds")


####
# end of temp edit 
####

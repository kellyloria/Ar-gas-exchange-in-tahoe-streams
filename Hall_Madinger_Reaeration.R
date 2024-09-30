# Supplement and code to Use of argon to measure gas exchange 
# in turbulent mountain streams
# link: file:///Users/kellyloria/Documents/UNR/Reaeration/SOP_literature/bg-15-3085-2018-supplement/Argon%20supplement/Arcode.html

library(rstan)
library(shinystan)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
ardata <- read_csv("/Users/kellyloria/Documents/UNR/Reaeration/SOP_literature/bg-15-3085-2018-supplement/Argon\ supplement/ardata.csv")
sf6data <- read_csv("/Users/kellyloria/Documents/UNR/Reaeration/SOP_literature/bg-15-3085-2018-supplement/Argon\ supplement/sf6data.csv")

# First need to calculate saturation concentrations. 
# These are based on Hamme and Emerson (2004)
watdens<-function(temp){
  
  t<-temp
  
  A <- 7.0132e-5
  B <- 7.926295e-3 
  C <-  -7.575477e-5 
  D<- 7.314701e-7
  E <-  -3.596363e-9
  to<- 3.9818
  
  dens<- (999.97358- (A*(t-to) + B*(t-to)^2 +C*(t-to)^3 + D*(t-to)^4+E*(t-to)^5) ) -4.873e-3 + 1.708e-4*t - 3.108e-6 * t^2
  dens/1000
}

nsat<- function(temp, bp) {
  
  
  ts<-log((298.15-temp) / (273.15 + temp))
  a0<-6.42931
  a1<-2.92704
  a2<-4.32531
  a3<-4.69149
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  satn<-(exp(a0 + a1*ts + a2*ts^2 + a3*ts^3))*((bp-u)/(760-u))
  watdens(temp)*satn*(28.014/1000)##converts umol/kg to mg/L
}

arsat<- function(temp, bp) {
  
  
  ts<-log((298.15-temp) / (273.15 + temp))
  a0<-2.79150
  a1<-3.17609
  a2<-4.13116
  a3<-4.90379
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  satar<-(exp(a0 + a1*ts + a2*ts^2 + a3*ts^3))*((bp-u)/(760-u))
  watdens(temp)*satar*(39.948/1000)##converts umol/kg to mg/L
}


## Estimates K600 for KO2 at a given temperature. From Wanninkhof (1992).
K600fromO2<-function (temp, KO2) {
  ((600/(1800.6 - (120.1 * temp) + (3.7818 * temp^2) - (0.047608 * temp^3)))^-0.5) * KO2
}


## Estimates K600 for KAr at a given temperature. From Raymond et al  (2012).
K600fromAr<-function (temp, KAr) {
  ((600/(1799 - (106.96 * temp) + (2.797 * temp^2) - (0.0289 * temp^3)))^-0.5) * KAr
}

###
ardata$arnsat <- arsat(ardata$temp,ardata$pressure) / nsat(ardata$temp,ardata$pressure)
ardatapre<- ardata[ardata$type=='pre', ]
ardatapost<- ardata[ardata$type=='post', ]
sf6datapost<- sf6data[sf6data$type=='post', ]

#Plot of raw Ar:N2. 
#Note the high variability in some of the pre Ar data. 
#We used the calculated Ar:N2 since it was in the middle of our samples, but not nearly as variable. 
#The two lines on the Ar sat are the post and pre calculations; they vary because of temperature. 
#The upper line is the warmer, pleateau collection temperature.

plot <- ggplot(data = ardata, aes(stationcorr, arncalc, color = as.factor(type)))+
  geom_point(size= 1, alpha = 0.6)+
  theme_bw() + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = 'right', 
        legend.direction = "vertical") + facet_grid(Trial~.)


# Now for calculations where we subtract background and correct for conductivity

ardatapost$arcorr<- ardatapost$arncalc - ardatapost$arnsat

#2.2.  Calc mean of all pre cond by station and by site
ardataprecond <- ardatapre %>% 
  group_by(Trial, stationcorr) %>%  # not really sure why station cor instead of station ?
  summarise(precond=mean(cond, na.rm=T)) # cond = SPC

#join with the post cond
ardatapost<-merge(ardatapost,ardataprecond)
hist(ardatapost$cond)

ardatapost$condcor<- ardatapost$cond - ardatapost$precond
hist(ardatapost$condcor)

ardatapost$arcond<- ardatapost$arcorr / ardatapost$condcor
hist(ardatapost$arcorr)
hist(ardatapost$arcond)

##calc the mean for station 0
ardatapost_0<-ardatapost[ardatapost$stationcorr==0,]
ardata_0sum <- ardatapost_0 %>% 
  group_by(Trial) %>% 
  summarise(arcond_0=mean(arcond), arn_enrich=mean(arncalc/arnsat))


#join with ardatapost
ardatapost<-merge(ardatapost,ardata_0sum)

ardatapost$arcondnorm<-ardatapost$arcond/ardatapost$arcond_0

############################
# Prep SF6 in same way
#2.2.  Calc mean of all pre cond by station and by site
ardataprecond <- ardatapre %>% 
  group_by(Trial, stationcorr) %>% 
  summarise(precond=mean(cond, na.rm=T))

#join with the post cond
sf6datapost<-merge(sf6datapost,ardataprecond) #yes, using pre from Ar



sf6datapost$condcor<- sf6datapost$cond - sf6datapost$precond

sf6datapost$sf6cond<- sf6datapost$sf6 / sf6datapost$condcor

##calc the mean for station 0
sf6datapost_0<-sf6datapost[sf6datapost$stationcorr==0,]
sf6data_0sum <- sf6datapost_0 %>% group_by(Trial) %>% summarise(sf6cond_0=mean(sf6cond))


#join with sf6datapost
sf6datapost<-merge(sf6datapost,sf6data_0sum)

sf6datapost$sf6condnorm<-sf6datapost$sf6cond/sf6datapost$sf6cond_0

###################################
#Stan model
#Here is the Stan model that we used for the analysis. Details on priors etc are in the text.
sink("combineargonmodelb.stan")

cat("
    
    data {
    
    int <lower=1> N; //number of data points, ar
    int <lower=1> NS; //number of data points, sf6
    int <lower=1> sites; //number of sites
    int <lower=0> sitear[N]; //stream number  for ar indexing vector
    int <lower=0> sitesf6[NS]; //stream number for sf6 indexing vector
    
    
    vector [N] ar;//conductivity corrected Ar data
    vector [N] distar;//
    
    vector [NS] sf6;//conductivity corrected sf6 data
    vector [NS] distsf6;// 
    
    }
    
    parameters {
    
    vector <lower=0> [sites] a;
    real mu_a;   //take out for hyperprior info  // mean prior
    real<lower=0, upper=0.5> sigma_ar; // error ar
    real<lower=0, upper=0.5> sigma_sf6; // error sf6
    
    vector[sites] k; // decline
    real <lower=0, upper=2> Ar0;
    real <lower=0, upper=2> SF60;
    
    //real d;
    //real b;
    real <lower=0> sigma_a; // mean prior
    }
    
    model { 
    
    //priors. 
    k ~ normal(0, 10);
  
    a~normal (mu_a,sigma_a); // mean prior

    Ar0 ~normal(1,0.05);
    SF60~normal(1,0.05);

    mu_a~normal(1.35, 1); // mean prior
     sigma_a ~ normal(0, 2);
    
    //likelihood        
    for (i in 1:N) {
    ar[i] ~ normal( Ar0 * exp(-k[sitear[i]]*distar[i]*0.01) , sigma_ar); 
    }
    for (j in 1:NS) {
    sf6[j] ~ normal( SF60 * exp(-k[sitesf6[j]]*distsf6[j]*0.01/a[sitesf6[j]]) , sigma_sf6); 
    }
    
    }

generated quantities {   //These estimate the posterior predicted

    vector [N] ar_tilde;
    vector [NS] sf6_tilde;
    
    for (n in 1:N) ar_tilde[n] = normal_rng( Ar0 * exp(-k[sitear[n]]*distar[n]*0.01) , sigma_ar);
    for (p in 1:NS) sf6_tilde[p] = normal_rng(SF60 * exp(-k[sitesf6[p]]*distsf6[p]*0.01/a[sitesf6[p]]) , sigma_sf6);
    
}

    "
,fill=TRUE)
sink()


arstandata=list("ar"= ardatapost$arcondnorm, "distar"=ardatapost$stationcorr, "N"=length(ardatapost$stationcorr), 
                "sf6"= sf6datapost$sf6condnorm, "distsf6"=sf6datapost$stationcorr, "NS"=length(sf6datapost$stationcorr), 
                "sites"=max(ardatapost$Trial),  "sitear" = ardatapost$Trial, 
                "sitesf6" = sf6datapost$Trial)

arfit <- stan(file = "combineargonmodelb.stan", data = arstandata,
              iter = 1000, chains = 3)





print(arfit, pars=c("a","mu_a","sigma_ar", "sigma_sf6", "k", "Ar0", "SF60", "sigma_a" ))

arfitsum<- summary(arfit)$summary

asum<-arfitsum[1:8,c("2.5%", "50%", "97.5%")]

ksum<-(arfitsum[12:19,c("2.5%", "50%", "97.5%")])*0.01# the 0.01 here rescales the k estimate

#. Essentailly we are using the model output to peredict a new set of data which we compare with our actual data.

ar_tildesum<- summary(arfit, pars="ar_tilde", probs=0.5)$summary
ar_tildesum<-ar_tildesum[,4]

sf6_tildesum<- summary(arfit, pars="sf6_tilde", probs=0.5)$summary
sf6_tildesum<-sf6_tildesum[,4]


par( mfrow=c(2,1),  mai=c(0.7,0.7,0.2,0.3), mgp=c(2,0.7,0) )
plot(ardatapost$arcondnorm, ar_tildesum, xlab="Measured normalized Ar concentration", ylab="Predicted normalized Ar concentration", pch=16, cex.lab=0.8, cex.axis=0.8)

plot(sf6datapost$sf6condnorm, sf6_tildesum, xlab="Measured normalized SF6 concentration", ylab="Predicted normalized SF6 concentration", pch=16, cex.lab=0.8, cex.axis=0.8)

# Calculate derived variables such as K(per time) and k600
streamslope<-c( 0.00703, 0.00703,0.015,0.015, 0.06, 0.11,0.12,0.12)
Q<- c(0.084,0.07,0.02,0.02,0.021,0.097,0.022,0.022)*60
w<- c(2.3,1.6,0.8,0.8,0.9,3.3,0.7,1.3)

ardataposttemp<-ardatapost %>% 
  group_by(Trial) %>% 
  summarise(temp=mean(temp, na.rm=T))
k<- ksum*Q*1440/w
k600<-K600fromAr(ardataposttemp$temp, k)
v<-c(12,15.4,9.5,9.5,3.1,4,6.3,5.2)

z<-(Q)/(w*v)

Kt<-v*ksum[,2]*1440

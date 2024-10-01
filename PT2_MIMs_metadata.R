## Merges MIMS data with meta data

## Goal create 1 data table with columns for 
## (Site, Date, Trial, Station, Rep, arncalc, type, cond, pressure, temp,	arconc,	nconc)

setwd("/Users/kellyloria/Documents/UNR/Reaeration/MIMS_dat/")
###
proc_dat <- read.csv("2024_KAL_prelimdat_processed.csv")%>%
  dplyr::select("Samp", "SampleID", "Pressure", "N2.Ar", )

meta_dat <- read.csv("kelly_prelim_sample_metadata.csv") %>%
  mutate(Samp=as.integer(MIMs_label))

head(proc_dat)
# 
# ## Possible cleaning steps
# ## 1. Average by SampleID and Samp
# proc_datq <- proc_dat %>%
#   group_by(SampleID, Calibnum) %>%
#   summarise(Pressure_m=mean(Pressure, na.rm=T),
#             Temp_m=mean(Temp, na.rm=T),
#             WatDens_m=mean(WatDens, na.rm=T),
#             O2Sat_m=mean(O2Sat, na.rm=T),
#             N2Sat_m=mean(N2Sat, na.rm=T),
#             ArSat_m=mean(ArSat, na.rm=T),
#             O2.ArSat_m=mean(O2.ArSat, na.rm=T),
#             N2.ArSat_m=mean(N2.ArSat, na.rm=T),
#             HeSat_m=mean(HeSat, na.rm=T),
#             #KrSat_83=mean(KrSat_83, na.rm=T),
#             #KrSat_84=mean(KrSat_84, na.rm=T),
#             He.ArSat_m=mean(He.ArSat, na.rm=T),
#             #He.ArSat=mean(Ar.KrSat_83, na.rm=T),
#             #Ar.KrSat_84_m=mean(Ar.KrSat_84, na.rm=T),
#             X28_m=mean(X28, na.rm=T),
#             X32_m=mean(X32, na.rm=T),
#             X40_m=mean(X40, na.rm=T),
#             X99_m=mean(X99, na.rm=T),
#             X29_m=mean(X29, na.rm=T),
#             X34_m=mean(X34, na.rm=T),
#             X30_m=mean(X30, na.rm=T),
#             N2.Ar_m=mean(N2.Ar, na.rm=T),
#             O2.Ar_m=mean(O2.Ar, na.rm=T),
#             X29.28_m=mean(X29.28, na.rm=T),
#             X34.32_m=mean(X34.32, na.rm=T),
#             X30.28_m=mean(X30.28, na.rm=T))
#             
            
## 2. Merge with metadata 
meta_datq <- proc_dat%>%
  full_join(meta_dat, by = c("Samp"))

## Look at some values with distance from injection:
meta_datq %>%
  filter(site=="GBL" & sample_type=="POST")%>%
  ggplot(aes(x = Station, y = c(X40/X28) , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 



########
# You are here:
######## 

meta_datq$Ar.N2_KL <- (meta_datq$X40/meta_datqX28)

meta_datq %>%
  filter(site=="BWL")%>%
  ggplot(aes(x = distance, y = N2.ArSat_m , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 
  #scale_shape_manual(values=c(19,3)) + facet_grid(order_group~.)

meta_datq %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = distance, y = O2Sat_m , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

meta_datq %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = distance, y = N2Sat_m , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

meta_datq %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = distance, y = ArSat_m , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

meta_datq %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = distance, y = O2.ArSat_m , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 

meta_datq %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = distance, y = X40 , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") +
  geom_smooth(se=F, lty=2)

meta_datq %>%
  filter(site=="GBL")%>%
  ggplot(aes(x = distance, y = X28 , color = sample_type, shape=sample_rep)) +
  geom_point(size=1, alpha=0.75) + theme_bw() + theme(legend.position = "right") 


### save metadata merge file for list format:
# saveRDS(meta_datq, file = "./model_data/meta_datq.rds")

## FORMAT into list for stan :
stan_data_list <- list(
  N = nrow(meta_datq), # total number of observations
  temp = mean(na.omit(meta_datq$temp..C.)), # temperature data
  DO = meta_datq$DO..mgL.,  # dissolved oxygen
  SPC = meta_datq$SPC..uScm., # specific conductance
  pH = meta_datq$pH, # pH values
  # Add other required variables
  # e.g., barometric pressure, station number, etc.
  baro = meta_datq$baro..mmHg.,
  distance = meta_datq$distance
)


# data {
#   int<lower = 1> N;
#   int<lower = 1> nexpt; # no. of experiments 
#   int<lower = 1> exptID[N]; # experiments id 
#   vector[N] dist;
#   vector[N] Ar; # proportion relative to upstreammost sampling station
#   vector[nexpt] Q; 
#   vector[nexpt] V;
#   vector[nexpt] temp; # mean for entire reach 
#   vector[nexpt] w; # mean for entire reach 
# }

# If you need to convert factors or dates
stan_data_list$date <- as.numeric(as.Date(meta_datq$date, format="%m/%d/%y")) # Example conversion

# Make sure other categorical variables are converted to numeric if needed for Stan
stan_data_list$site <- as.numeric(factor(meta_datq$site))



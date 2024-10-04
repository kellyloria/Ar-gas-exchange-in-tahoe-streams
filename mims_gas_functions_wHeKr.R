
############################################################################################################
# Gas functions to estimate saturation concentrations
# Bob Hall  
# Sept 2015

# This code includes calculations for barometric pressure (bp), water density, oxygen saturation, nitrogen saturation, and argon saturation.
# If you have bp in mmHg from a handheld sensor you do not need to calculate the bp again. 
# Execute these functions first before using MIMS_data_function_example.R 
############################################################################################################

###########################################
# Correct BAROMETRIC PRESSURE for altitude
###########################################	

### bpcalc function estimates barometric pressure (BP) given altitude and standard BP
### This is based on the barometric formula
### temp = degC for the water
### alt = elevation in m
### bpst = standard barometric pressure in inches of Hg (usually given by weather websites).  
### Temp is usually relative to a standard, 15 degC, and I have hard coded it at 15.  
### Source: John Colt, Dissolved Gas Concentration in Water: Computation as Functions of Temperature, Salinity and Pressure, 2012

# bpcalc<- function(bpst, alt) {
# 
# bpst*25.4*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+15)))
# 
# }
# 
# 
# bpcalc(bpst=29.92,alt=2400)
# ##or as 
# bpcalc(29.92, 2400)


###########################################
# Calculate WATER DENSITY w/ Temp and BP
###########################################	

### Water density of air saturated water given water temperature in degC.  
### Source: Paterson and Morris 1994, Meterologia
watdens<-function(temp){
	
	t<-temp
	# coeffs to model temperature and water density
A <- 7.0132e-5
B <- 7.926295e-3 
C <-  -7.575477e-5 
D<- 7.314701e-7
E <-  -3.596363e-9
# max density reference temperature
to<- 3.9818

# calculates the density of water (in kg/m³) based on the temperature difference to reference temp 
dens<- (999.97358- (A*(t-to) + B*(t-to)^2 +C*(t-to)^3 + D*(t-to)^4+E*(t-to)^5) ) -4.873e-3 + 1.708e-4*t - 3.108e-6 * t^2
dens/1000
}

##call as
watdens(20)


###########################################
# Calculate O2 SATURATION w/ Temp and BP
###########################################	

### Oxygen saturation. Alternative formulation from Garcia and Gordon (umol/kg), which is converted to mg/L and corrected for water density.  
### This function gives the same values as from Colt and is the one a MIMSer should use.
### u is the vapor pressure of water
### Ending units mg/L
osat1<- function(temp, bp) {
	u<-10^(8.10765-(1750.286/(235+temp)))
	ts<-log((298.15-temp) / (273.15 + temp))
	a0<-5.80871
	a1<-3.20291
	a2<-4.17887
	a3<-5.1006
	a4<- -9.86643e-2
	a5<-3.88767
	
	u<-10^(8.10765-(1750.286/(235+temp)))
	sato<-(exp(a0 + a1*ts + a2*ts^2 + a3*ts^3 + a4*ts^4)+a5*ts^5)*((bp-u)/(760-u))
	watdens(temp)*sato*(31.9988/1000)##converts umol/kg to mg/L
	
	}

##call as
osat1(20, 590.2)


###########################################
# Calculate N2 SATURATION w/ Temp and BP
###########################################	
### Nitrogen saturation. From Hamme and Emerson 2004 Deep Sea Res.
### Ending units mg/L

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

##call as
nsat(20,760)


###########################################
# Calculate Ar SATURATION w/ Temp and BP
###########################################	
### Argon saturation. From Hamme and Emerson 2004 Deep Sea Res.
### Checked against values in Colts table and correct.
### Units mg/L   

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

##call as
arsat(10,760)


###########################################
# Calculate Kr SATURATION w/ Temp and BP
###########################################	
### Krypton saturation. 

### Units mg/L   


krsat_83<- function(temp, bp) {
  
  # Temperature Transformation into scaled logarithmic form (ts)
  T<-(273.15 + temp)
  # Solubility Coefficients for Argon (ax)
  a0<- -122.4694
  a1<- 153.5654
  a2<-  70.1969
  a3<- -8.5224
  
  # Water vapor pressure calculation (u)
  # u reduces the amount of Ar that can be dissolved in water because it takes up some space in the gas phase
  u<-10^(8.10765-(1750.286/(235+temp)))
  # Ar saturation: the concentration of Ar that would be dissolved in water at a given temperature, barometric pressure, and u 
  satkr<-(exp(a0 + a1*(100/T) + a2*log(T/100) + a3*(T/100)  ))*((bp-u)/(760-u))
  watdens(temp)*satkr*(83/0.001)#### converts Ar from mol/kg to mg/L
}


krsat_84<- function(temp, bp) {
  
  
  T<-(273.15 + temp)
  a0<- -122.4694
  a1<- 153.5654
  a2<-  70.1969
  a3<- -8.5224
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  satkr<-(exp(a0 + a1*(100/T) + a2*log(T/100) + a3*(T/100)  ))*((bp-u)/(760-u))
  watdens(temp)*satkr*(84/0.001)##converts mol/kg to mg/L
}

krsat_83(20,760)


Temp<-seq(1:20)

plot(arsat(Temp, 760) / krsat_83(Temp, 760))


###########################################
# Calculate He SATURATION w/ Temp and BP
###########################################	
### Helium saturation.
### Equation from Jenkins et al. 2019
### Units mg/L   

hesat<- function(temp, bp) {
  
  
  
  
  T<- (273.15 + temp)
  a1<- -178.1424
  a2<-217.5991
  a3<-140.7506
  a4<- -23.01954
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  sathe<-(exp(a1 + a2*(100/T) + a3*log(T/100) + a4*(T/100) ))*((bp-u)/(760-u))##Concentration in mol/kg
  watdens(temp)*sathe*(4/0.001)##converts mol/kg to mg/L
}

hesat(20,760)

Temp<-seq(1:20)

plot(hesat(Temp, 760) / arsat(Temp, 760))








###########################################	
### Examples of plotting saturation info. 
### Plot to see how O2/Ar changes with temp. Note that these are mass and not molar.
temp<- seq(from=0, to=25, by=1)
#plot(temp, osat1(temp,760)/arsat(temp,760))

#plot(temp, (nsat(temp,571)/28)/(arsat(temp,571)/40)  )
#plot(temp, (nsat(temp,571))/(arsat(temp,571))  )

#(nsat(20,571)/28.014)/(arsat(20,571)/39.948)


#points(10,38.41, pch=16)
#points(16,38.297, pch=16)





###########################################
# GAS EXCHANGE (K) functions
###########################################	

### The gas exchange functions here are not required for the MIMS_datasheet_function_example.R code, but are useful in many other related projects. 

### Gas exchange is sensitive to temperature. 
### This code converts K600 to a value at a specific temp by Schmidt number scaling. From J?hne et al. (1987).
### K600 to KO2
Kcor<-function (temp,K600) {
	K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
	}
	
### K600 to KN2
KcorN<-function (temp,K600) {
	K600/(600/(1970.7-(temp*131.45)+(4.139*temp^2)-(0.052106*temp^3)))^-0.5
	}

### K600 to KAr
KcorAr<-function (temp,K600) {
	K600/(600/(1759.7-(temp*117.37)+(3.6959*temp^2)-(0.046527*temp^3)))^-0.5
	}


### An example of using the above functions: 
### Fort Lupton nsat, elevation 1496 m
#bpcalc(bpst=29.92,alt=1496)
#636.4546

fltemp<-c(10,11,12,13,14)
#nsat(fltemp,638.4)*1000/28
# [1] 547.5234 535.8124 524.5701 513.7720 503.3951


## Estimates K600 for KO2 at a given temperature. From Wanninkhof (1992).
K600fromO2<-function (temp, KO2) {
  ((600/(1800.6 - (120.1 * temp) + (3.7818 * temp^2) - (0.047608 * temp^3)))^-0.5) * KO2
}

## Estimates K600 for KAr at a given temperature. From Raymond et al  (2012).
K600fromAr<-function (temp, KAr) {
  ((600/(1799 - (106.96 * temp) + (2.797 * temp^2) - (0.0289 * temp^3)))^-0.5) * KAr
}

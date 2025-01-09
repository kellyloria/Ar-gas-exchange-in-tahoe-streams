# Ar-gas-exchange-in-tahoe-streams
Estimating first order decay of Ar:N2 gas injections in streams to get K600 

## PT1 
Processes raw MIMS files and sample data files to calculate saturation and concentration values of N2 and Ar
Workflow developed by Bob Hall, Tyler Tappenbeck and others at FBLS 

## PT2 
Uses an expontential decay relationship to model K600 using Stan model for Ar decline for W1 developed by Bob Hall. 
- Gas exchange model developed in Aho et al. 2024 (link - https://essd.copernicus.org/articles/16/5563/2024/)
- Where the model estimates log(Kt), which is per time rate of exchange. 
- We chose to not pool multiple Ar releases.
    -  And so there are two scripts for modeling the multiple releases at BWL or GBL to adjust for the location of the first best mixed station.
- Ar must be entered as proportion based on mean value from first site. 
- Code will work with Ar or Ar/N2 since it takes a proprtion. 

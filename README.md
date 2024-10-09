# Ar-gas-exchange-in-tahoe-streams
Estimating first order decay of Ar:N2 gas injections in streams to get K600 

## PT1 
Processes raw MIMS files and sample data files to calculate saturation and concentration values of N2 and Ar
Workflow developed by Bob Hall, Tyler Tappenbeck and others at FBLS 

## PT2 
Uses an expontential decay relationship to model K600 using Stan model for Ar decline for W1 developed by Bob Hall. 
- Where the model estimates log(Kt), which is per time rate of exchange, 
- Pools multiple releases based on distance and assumed mean velocity.
- Ar must be entered as proportion based on mean value from first site. 
- Code will work with Ar or Ar/N2 since it takes a proprtion. 

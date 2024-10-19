
data {
  int<lower = 1> N;
  int<lower = 1> nexpt;
  int<lower = 1> exptID[N];
  vector[N] dist;
  vector[N] Ar; // need logAr
  vector[nexpt] Q;
  vector[nexpt] V;
  vector [nexpt] temp; // one value for each stream temp to convert to K600
  vector[nexpt] w; 
  
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
  
  for (j in 1: nexpt){  // make a loop here // very these smitts 
    KAr[j] = exp(logK600[j]) /  ((600/(1759.7-(117.37*temp[j])+(3.6959*temp[j]^2)-(0.046527*temp[j]^3)))^-0.5); 
  }
  
  Kd= KAr ./ V;  //KAr (1/day) / V (m/day) = Kd (1/m)
  
  
}

model {
  for (i in 1:N){
    Ar[i] ~ normal(intercept + -Kd[exptID[i]]*dist[i], sigma); // likelihood 
    // if you know travel time you solve for travel time here instead of distance 
  }
  
  for (j in 1:nexpt){
    logK600[j]~normal( a + b*log(Q[j]) , sigma_expt);
  }
  
  a ~ normal(0,10);
  b ~ normal(0,1);
  sigma_expt ~ normal (0,2);
  sigma ~ normal (0,0.2);  // added prior on sigma
  intercept ~ normal (0, 0.1); // working with proportions so strong prior on intercept
}
    

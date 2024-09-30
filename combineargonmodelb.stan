
    
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

    

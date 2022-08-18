## RW with Intercept for Pois Distributed Variable###
data{
  int <lower=1> T; //Time Index
  int <lower=1> A; // Age Index
  int <lower=1> R; //Region Index
  int <lower=0> y[T*A*R]; // count outcomes
  vector<lower=0>[T*A*R] E;// exposure
  
  // Indicies for Data 
  int TInd[T*A*R];  // Index of Time in Data
  int AInd[T*A*R];  // Age Index in Data
  int ReInd[T*A*R]; // Index of Random Effects
  int RInd[T*A*R]; // Index Region Effects
  
  int <lower=0> TFor; //Number of years to forecast 
  
  // Neighbor for Region Effects
  int<lower=0> N_edges;
  int<lower=1, upper=R> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=R> node2[N_edges];  // and node1[i] < node2[i]
  
  real<lower=0> scaling_factor; // scales the Variances of the spatial effects
  
} transformed data {
  int N = T*A*R; //Total Number of Observations
  vector[N] log_E = log(E); // log of exposure

  int<lower = 1> L; // size of prediction vector
  L=TFor*A*R; //  
}
parameters{
  real Intercept; // intercept
  real drift; // drift Parameter
  
  vector[A] beta_age; // RW Age Coefficients
  
  vector[T] kappa_time; // RW Time Coefficients 

  //Variance Parameters
  real<lower=0> sigma_time;    // Variance of Time
  real<lower=0> sigma_age;     // Variance of Age
  real<lower=0> sigma_eps;     //Variance of epsilon
  
  real<lower=0> sigma_BYM2; //overall Standard deviation in BYM2 Model
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  
  vector[R] v;       // heterogeneous effects (non spatial)
  vector[R] u;         // spatial effects
  
  vector[N] eps;     // Vector of overdispersion effects
  vector[N] eps_new; // new vector of overdisp effects for PPC
} 
transformed parameters {

  vector[N] InterceptVek = rep_vector(Intercept,N);
  
  vector[R] convolved_re; // Total BYM2 Effect
  vector[R] phi_region; // Total regional effect
  
  // variance of each component should be approximately equal to 1
  convolved_re =  sqrt(1 - rho) * v + sqrt(rho / scaling_factor) * u; // see Riebler et. al (2017)
  phi_region = convolved_re*sigma_BYM2; // (non centered Version)

} 
model {
   // Muss ganz am Anfang stehen (Stan only supports variable definitions at top of the block)
   vector[N] mu; // Vector of random effects
  mu[ReInd]= kappa_time[TInd]+beta_age[AInd]+phi_region[RInd];
  
  //Varianzen
  target += student_t_lpdf(sigma_age|5,0,1);
  target += student_t_lpdf(sigma_time|5,0,1);
  target += student_t_lpdf(sigma_eps|5,0,1);
  target += student_t_lpdf(sigma_BYM2|5,0,1);
  
  //Other Effects
  target += normal_lpdf(Intercept|-5,5);
  target += normal_lpdf(drift|0,2);
  target += normal_lpdf(eps|0,1);
  
  target += std_normal_lpdf(eps_new); // 
  
  //Prior for rho
  target += beta_lpdf(rho|0.5,0.5); // prior as proposed by Riebler et. al 
  
  
  // Random Walk with Drift Prior (gives invariant Forecasts, see Huang et. al(2008))
  target += normal_lpdf(kappa_time[1]|drift,sigma_time);
  target += normal_lpdf(kappa_time[2:T]|drift+kappa_time[1:(T- 1)],sigma_time);    // Random walk with drift prior
  
  sum(kappa_time) ~ normal(0,0.001*T); //soft sum to zero constraint  // Random walk with drift prior
  
  //Age Effect
  target += normal_lpdf(beta_age[1]|0,sigma_age); // First Effects seperate prior (see e.g. Fosse(2018))
  target += normal_lpdf(beta_age[2]|0,sigma_age);
  //RW(2) Prior
  target+= normal_lpdf(beta_age[3:A]| 2*beta_age[2:(A-1)]-beta_age[1:(A-2)], sigma_age);
  
  sum(beta_age) ~ normal(0,0.001*A); //soft sum to zero constraint
  
   //Region Effect (ICAR Model)
  target += -0.5 * dot_self(u[node1] - u[node2]);
  
  // soft sum-to-zero constraint on u)
  sum(u) ~ normal(0, 0.001 * R);  // equivalent to mean(u) ~ normal(0,0.001)
  
  target += normal_lpdf(v|0,1);

  target += poisson_log_lpmf(y|log_E + InterceptVek + mu + eps*sigma_eps);

  
} generated quantities { // posterior
  
  vector[N] log_like_y; // create Vector to store log likelihood (for waic (see Vethari, Gabry (2019)))
  
  vector[N] lambdahat = exp(log_E + InterceptVek + //Log plus Intercept
                         kappa_time[TInd]+beta_age[AInd]+phi_region[RInd]+ //mu
                         eps_new*sigma_eps); // Error Term
  
  vector[TFor] kappa_time_pred;
  vector[L] mufor; // predicted death rates
  int pos_f = 1;
  int pos_L1 =1; 

  
   kappa_time_pred[1] = drift+kappa_time[T]+sigma_time * normal_rng(0,1);
  // Check if Number of Years to Forecast is greater than 1
   if(TFor > 1){
     //RW (1) - Drift Cohort Index
    for (t in 2:TFor) kappa_time_pred[t] = drift+kappa_time_pred[t - 1] + sigma_time * normal_rng(0,1);
 }
 
  // Forecast Log Mortality Rates
  for (t in 1:TFor) for (r in 1:R) for (a in 1:A)  {
    mufor[pos_f] =  Intercept+
                    beta_age[a]+kappa_time_pred[t]+phi_region[r]+
                    sigma_eps*normal_rng(0,1);
    pos_f += 1;
  }
   // Log Likelihood
  for (t in 1:T) for (r in 1:R) for (a in 1:A){
    log_like_y[pos_L1] = poisson_log_lpmf(y[pos_L1]| log_E[pos_L1] + Intercept + 
                                                 kappa_time[t]+beta_age[a]+phi_region[r]+
                                                 eps[pos_L1]*sigma_eps);
    pos_L1 += 1;
  }
  
}

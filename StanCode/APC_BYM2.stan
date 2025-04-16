// APC_BYM2 Model
data{
  int <lower=1> T; //Time Index
  int <lower=1> A; // Age Index
  int <lower=1> R; //Region Index
  int <lower=0> y[T*A*R]; // count outcomes
  vector<lower=0>[T*A*R] E;// exposure
  int <lower=1> M; // How much wider is Age interval to Time 
  
  // Indicies for Data 
  int TInd[T*A*R];  // Index of Time in Data
  int AInd[T*A*R];  // Age Index in Data
  int ReInd[T*A*R]; // Index of Random Effects
  int RInd[T*A*R]; // Index Region Effects
  int CInd[T*A*R]; // Index for Cohort Effects
  
  int <lower=0> TFor; //Number of years to forecast 
  
  // Neighbor for Region Effects
  int<lower=0> N_edges;
  int<lower=1, upper=R> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=R> node2[N_edges];  // and node1[i] < node2[i]
  
  real<lower=0> scaling_factor; // scales the Variances of the spatial effects
  
} transformed data {
  int N = T*A*R; //Total Number of Observations
  vector[N] log_E = log(E); // log of exposure
  
  int<lower=1> C;// Cohort Index
  int<lower = 1> L; // size of prediction vector
  
  C = (M* (A-2)) + T; // max Cohort Index
  L=TFor*A*R; // Length of prediction vector
  
}
parameters{
  real Intercept; // intercept
  real drift; // drift Parameter
  
  vector[A] alpha_age; //  Age Coefficients
  
  vector[C] gamma_cohort;// Cohort Coefficients
  
  vector[T] kappa_time; // Time Coefficients

  //Variance Parameters
  real<lower=0> sigma_time;    // Standard deviation of Time
  real<lower=0> sigma_age;     // Standard deviation of Age
  real<lower=0> sigma_eps;     // Standard deviation of epsilon
  real<lower=0> sigma_cohort;  // Standard deviation of Cohort Effect
  
  real<lower=0> sigma_BYM2; //overall Standard deviation in BYM2 Model
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  
  vector[R] v;       // heterogeneous effects (non spatial)
  vector[R] u;         // spatial effects
  
  vector[N] eps;     // Vector of overdispersion effects
  vector[N] eps_new; // new vector of overdisp effects for PPC
} 
transformed parameters {
  vector[N] InterceptVek = rep_vector(Intercept,N); //Vector of Intercept
  
  vector[R] convolved_re; // Total BYM2 Effect
  vector[R] phi_region; // Total regional effect
  
  // variance of each component should be approximately equal to 1
  convolved_re =  sqrt(1 - rho) * v + sqrt(rho / scaling_factor) * u; // see Riebler et. al (2017)
  phi_region = convolved_re*sigma_BYM2; // (non centered Version)
} 
model {
   // Needs to be at the beginning of the section (Stan only supports variable definitions at top of the block)
  vector[N] mu; // Vector of random effects
  mu[ReInd]= kappa_time[TInd]+alpha_age[AInd]+phi_region[RInd]+
             gamma_cohort[CInd];
  
  //Standard deviations (t distributed with 5 degrees of freedom)
  target += student_t_lpdf(sigma_age|5,0,1);
  target += student_t_lpdf(sigma_time|5,0,1);
  target += student_t_lpdf(sigma_eps|5,0,1);
  target += student_t_lpdf(sigma_cohort|5,0,1);
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
  
  sum(kappa_time) ~ normal(0,0.001*T); //soft sum to zero constraint
  
  //Age Effect
  target += normal_lpdf(alpha_age[1]|0,sigma_age); // First Effects seperate prior 
  target += normal_lpdf(alpha_age[2]|0,sigma_age);
  //RW(2) Prior
  target+= normal_lpdf(alpha_age[3:A]| 2*alpha_age[2:(A-1)]-alpha_age[1:(A-2)], sigma_age);
  
  sum(alpha_age) ~ normal(0,0.001*A); //soft sum to zero constraint
  
  // Cohort Effect
  target+=normal_lpdf(gamma_cohort|0,sigma_cohort); //IID normal prior
  
  sum(gamma_cohort) ~ normal(0,0.001*C); // Soft sum to zero constraint
  
  //Region Effect (ICAR Model)
  target += -0.5 * dot_self(u[node1] - u[node2]);
  
  // soft sum-to-zero constraint on u)
  sum(u) ~ normal(0, 0.001 * R);  // equivalent to mean(u) ~ normal(0,0.001)
  
  target += normal_lpdf(v|0,1);

  target += poisson_log_lpmf(y|log_E + InterceptVek + mu + eps*sigma_eps);
  
} generated quantities { 
  
  // posterior predictive checks
  vector[N] lambdahat = exp(log_E + InterceptVek + //Log plus Intercept
                         kappa_time[TInd]+alpha_age[AInd]+phi_region[RInd]+ gamma_cohort[CInd]+ //mu
                         eps_new*sigma_eps); // Error Term (new random draws)
  
  // Insample Rates
  vector[N] MHat = exp(InterceptVek + //Log plus Intercept
                         kappa_time[TInd]+alpha_age[AInd]+phi_region[RInd]+ gamma_cohort[CInd]+ //mu
                         eps*sigma_eps); // InSample Fit of estimated Mortality Rates (for Life Expectancy)
  
  vector[TFor] kappa_time_pred;
  vector[TFor] gamma_cohort_pred;
  vector[C+TFor] gamma_cohort_final;
  
  vector[L] mufor; // predicted death rates
  int pos_f = 1;
  int pos_L1 =1; 
  int <lower=0> K; // Kohort Index for Loop
  
  kappa_time_pred[1] = drift+kappa_time[T]+sigma_time * normal_rng(0,1);
  
  // Check if Number of Years to Forecast is greater than 1
   if(TFor > 1){
     //RW (1) - Drift Cohort Index
    for (t in 2:TFor) kappa_time_pred[t] = drift+kappa_time_pred[t - 1] + sigma_time * normal_rng(0,1);
 }
 // Forcast for Cohort Index
  for(t in 1:TFor){
    gamma_cohort_pred[t] = sigma_cohort*normal_rng(0,1);
  }
  // Append Forecasted Cohort Index on Existing Index
  gamma_cohort_final = append_row(gamma_cohort, gamma_cohort_pred);
  
  // Forecast Log Mortality Rates
  for (t in 1:TFor) for (r in 1:R) for (a in 1:A)  {
    if(a != 1 ){
      K = M*((A-1)-(a-1))+(T+t);  
    } else {
      K = M*((A-1)-(a))+(T+t); 
    }
    mufor[pos_f] =  Intercept+alpha_age[a]+kappa_time_pred[t]+
                    phi_region[r]+
                    gamma_cohort_final[K]+
                    sigma_eps*normal_rng(0,1);
    pos_f += 1;
  }
  
}

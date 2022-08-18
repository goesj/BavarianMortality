## RH Model with SAR###
functions{
  real normal_lagsar_lpdf(vector y, vector mu, real sigma, 
                          real rho, data matrix W, data vector eigenW) { 
    int N = rows(y);
    real inv_sigma2 = inv_square(sigma);
    matrix[N, N] W_tilde = add_diag(-rho * W, 1);
    vector[N] half_pred;
    real log_det;
    half_pred = W_tilde * y - mu;
    log_det = sum(log1m(rho * eigenW));
    return  0.5 * N * log(inv_sigma2) + log_det - 
      0.5 * dot_self(half_pred) * inv_sigma2;
  }
}
data{
  int <lower=1> T; //Time Index
  int <lower=1> A; // Age Index
  int <lower=1> R; //Region Index
  int <lower=0> y[T*A*R]; // count outcomes
  vector<lower=0>[T*A*R] E;// exposure
  int <lower=1> M; // How much wider is Age interval to Time 
  
  // Neighbor for Region Effects
  matrix<lower=0>[R,R] W; // Matrix of Weights
  vector[R] eigenWsar;  // eigenvalues of Msar
  
  int<lower = 1> TFor;                  // number of forecast years
  

  
} transformed data {
  int N = T*A*R; //Total Number of Observations
  
  vector[N] log_E = log(E); // log of exposure
  
  int<lower=1> C;// Maximum Cohort Index Value
  int<lower = 1> Pred; // size of prediction vector
  int<lower=1> k[N]; // Kohort Index Vector
  int pos_k=1; // Helper for Cohort Index
  
  // Calculate Cohort Index for each Observation
  for(t in 1:T) for(r in 1:R) for (a in 1:A){
    if(a==1){
      k[pos_k] = M*((A-1)-a)+t;
    } else { 
      k[pos_k] = M*((A-1)-(a-1))+t;
    }
    pos_k += 1;
  }
  
  C = (M* (A-2)) + T; // max Cohort Index
  Pred=TFor*A*R; //
}

parameters{

  vector[A] alpha_age; // Age Specific Isntercept
  simplex[A] beta_age1; //Interaction Term (sums to 1 and stricly positive)
  simplex[A] beta_age2; // Interaction Term for Cohort Effect
  
  vector[C] var_cohort;// Cohort Index (non Centered Para)
  
  vector[T] kappa_time; //  for corner constraint
  real drift; // drift Parameter
  
  //Variances 
  real<lower=0> sigma_time;    // Variance of Time
  real<lower=0> sigma_eps;     //Variance of epsilon
  real<lower=0> sigma_cohort; //
  real<lower=0> sigma_alpha; //
  
  real<lower=0> sigma_region; //overall Standard deviation in SAR Model
  real<lower=-1,upper=1> Xi; // Parameter in SAR Model (Matrix is row standardized, so that values are between 0 and 1)
  
  
  //Spatial Model (SAR)
  vector[R] phi;         // spatial effects
  
  
  vector[N] eps; // Vector of overdispersion effects
} transformed parameters {
  vector[C] gamma_cohort = var_cohort * sigma_cohort; // Non centered version
}
model {
   // Muss ganz am Anfang stehen (Stan only supports variable definitions at top of the block)
  vector[N] mu; // Vector of random effects
  int pos = 1;
  for(t in 1:T) for(r in 1:R) for (a in 1:A){
    mu[pos]=alpha_age[a]+beta_age1[a]*kappa_time[t]+beta_age2[a]*gamma_cohort[k[pos]]+
    phi[r];
    pos += 1; //pos = pos + 1 
   }
  
  //Varianzen
  target +=student_t_lpdf(sigma_time|5,0,1);
  target +=student_t_lpdf(sigma_eps|5,0,1);
  target +=student_t_lpdf(sigma_region|5,0,1);
  target +=student_t_lpdf(sigma_cohort|5,0,1);
  target +=student_t_lpdf(sigma_alpha|5,0,1);
  
  //Other Effects
  target += normal_lpdf(drift|0,2); // prior Intercept (medium vague)
  target += normal_lpdf(eps|0,1); //Epsilon Effects
  
   //Prior for Xi
  target += normal_lpdf(Xi|0,1);
  
  //TIME Effect
  // Random Walk with Drift Prior (see https://github.com/kabarigou/StanMoMo/blob/50ff0a4e5b2288d308c0b8f57cff5b21e72ddca6/inst/stan/leecarter.stan)
  target += normal_lpdf(kappa_time[1]|drift,sigma_time);
  target += normal_lpdf(kappa_time[2:T]|drift+kappa_time[1:(T- 1)],sigma_time);    // Random walk with drift prior
  
  target += normal_lpdf(sum(kappa_time)|0,0.001*R); // soft sum to zero 
  
  //Age Effect
  target += normal_lpdf(alpha_age[1:2]|0,sigma_alpha);
  target += normal_lpdf(alpha_age[3:A]|2*alpha_age[2:(A-1)]-alpha_age[1:(A-2)],sigma_alpha); //RW2 Prior
  
  target += dirichlet_lpdf(beta_age1|rep_vector(1, A));// Prior on beta_x
  target += dirichlet_lpdf(beta_age2|rep_vector(1, A));// Prior on beta_x
  
  //Cohort Effect (IID)
  target += std_normal_lpdf(var_cohort);
  
  //Region Effect (Besag Model)
  target += normal_lagsar_lpdf(phi | rep_vector(0.0, R), sigma_region, Xi, W, eigenWsar);
  
  // soft sum-to-zero constraint on phi)
  sum(phi) ~ normal(0, 0.001 * R);  // equivalent to mean(phi) ~ normal(0,0.001)
  
  
  target += poisson_log_lpmf(y| log_E + mu + eps*sigma_eps);
  
} generated quantities {
  
  // Quantities for InSample Fit
  vector[N] log_like_y; 
  vector[N] lambdahat; // Vector of random effects
  int pos = 1;
  
  // Quantites for Forecast
  vector[TFor] kappa_time_pred;
  vector[TFor] gamma_cohort_pred;
  vector[C+TFor] gamma_cohort_final;
  vector[Pred] mufor; // predicted death rates
  int pos_f = 1;
  int <lower=1> kFor; // Cohort Index for Forecast
  
  // get mu for insample fit
  for(t in 1:T) for(r in 1:R) for (a in 1:A){
    lambdahat[pos]=exp(alpha_age[a]+beta_age1[a]*kappa_time[t]+
                        beta_age2[a]*gamma_cohort[k[pos]]+
                        phi[r]+
                        log_E[pos]+normal_rng(0,1)*sigma_eps);
                        
    log_like_y[pos] = poisson_lpmf(y[pos]|exp(alpha_age[a]+beta_age1[a]*kappa_time[t]+
                                              beta_age2[a]*gamma_cohort[k[pos]]+
                                              phi[r]+eps[pos]*sigma_eps+
                                              log_E[pos]));
    pos += 1; //pos = pos + 1 
   }
  
  kappa_time_pred[1] = drift+kappa_time[T]+sigma_time * normal_rng(0,1);
  //gamma_cohort_pred[1] = 2*gamma_cohort[C-1]-gamma_cohort[C-2]+sigma_cohort*normal_rng(0,1);
  
  // Check if Forecast Period is greater than 1
  if(TFor > 1){
  //   // RW(2) FC Cohort Index
  //   gamma_cohort_pred[2] = 2*gamma_cohort_pred[1]-gamma_cohort[C-1]+sigma_cohort*normal_rng(0,1);
  //   for(t in 3:TFor){
  //     gamma_cohort_pred[t]= 2*gamma_cohort_pred[t-1]-gamma_cohort_pred[t-2]+sigma_cohort*normal_rng(0,1);
  //   }
  //   
  //   //RW (1) - Drift Cohort Index
   for (t in 2:TFor) kappa_time_pred[t] = drift+kappa_time_pred[t - 1] + sigma_time * normal_rng(0,1);
   }
  for(t in 1:TFor){
    gamma_cohort_pred[t] = sigma_cohort*normal_rng(0,1);
  }
  // Append Forecasted Cohort Index on Existing Index
  gamma_cohort_final = append_row(gamma_cohort, gamma_cohort_pred);
  
    // See that it has the same structure as loop in model block (else predictions are not comparable)
  for (t in 1:TFor) for (r in 1:R) for (a in 1:A)  {
     if(a==1){
      kFor = M*((A-1)-a)+(T+t);
    } else { 
      kFor = M*((A-1)-(a-1))+(T+t);
    }
    mufor[pos_f] = alpha_age[a] + 
                   beta_age1[a] * kappa_time_pred[t]+
                   beta_age2[a]*gamma_cohort_final[kFor]+
                   phi[r]+
                   sigma_eps*normal_rng(0,1);
    pos_f += 1;
  }

}

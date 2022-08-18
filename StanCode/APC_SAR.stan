// APC Model as described in Riebler,Held(2017) & Smith, Wakefield (2016) & Kuang et. al (2008)
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
  
  // Indicies for Data 
  int TInd[T*A*R];  // Index of Time in Data
  int AInd[T*A*R];  // Age Index in Data
  int ReInd[T*A*R]; // Index of Random Effects
  int RInd[T*A*R]; // Index Region Effects
  int CInd[T*A*R]; // Index for Cohort Effects
  
  int <lower=0> TFor; //Number of years to forecast 
  
  // Neighbor for Region Effects
  matrix<lower=0>[R,R] W; // Matrix of Weights
  vector[R] eigenWsar;  // eigenvalues of Msar
  
} transformed data {
  int N = T*A*R; //Total Number of Observations
  vector[N] log_E = log(E); // log of exposure

  
  int<lower=1> C;// Cohort Index
  int<lower = 1> L; // size of prediction vector
  
  C = (M* (A-2)) + T; // max Cohort Index
  L=TFor*A*R; //
  

}
parameters{
  real Intercept; // intercept
  real drift; // drift Parameter
  
  vector[A] beta_age; // RW Age Coefficients
  
  vector[C] var_cohort;// Cohort Index
  
  vector[T] kappa_time; // RW Time Coefficients (only T-1 since Corner constraint is applied)

  //Variance Parameters
  real<lower=0> sigma_time;    // Variance of Time
  real<lower=0> sigma_age;     // Variance of Age
  real<lower=0> sigma_eps;     //Variance of epsilon
  real<lower=0> sigma_cohort; //
  
  real<lower=0> sigma_region; //overall Standard deviation in SAR Model
  real<lower=-1,upper=1> xi; // Parameter in SAR Model (Matrix is row standardized, so that values are between -1 and 1)
  
  vector[R] phi_region;         // spatial effects
  
  vector[N] eps; // Vector of overdispersion effects
  vector[N] eps_new; // new vector of overdisp effects for PPC
} 
transformed parameters {

  vector[N] InterceptVek = rep_vector(Intercept,N);

  vector[C] gamma_cohort = var_cohort*sigma_cohort; // Total effect gamma cohort
} 
model {
   // Muss ganz am Anfang stehen (Stan only supports variable definitions at top of the block)
   vector[N] mu; // Vector of random effects
  mu[ReInd]= kappa_time[TInd]+beta_age[AInd]+phi_region[RInd]+
             gamma_cohort[CInd];
  
  //Varianzen
  target += student_t_lpdf(sigma_age|5,0,1);
  target += student_t_lpdf(sigma_time|5,0,1);
  target += student_t_lpdf(sigma_eps|5,0,1);
  target += student_t_lpdf(sigma_cohort|5,0,1);
  target += student_t_lpdf(sigma_region|5,0,1);
  
  //Other Effects
  target += normal_lpdf(Intercept|-5,5);
  target += normal_lpdf(drift|0,2);
  target += normal_lpdf(eps|0,1);
  target += std_normal_lpdf(eps_new);
  
  //Prior for rho
  target += normal_lpdf(xi|0,1); //trunc norm prior
  
  
  // Random Walk with Drift Prior (gives invariant Forecasts, see Huang et. al(2008))
  target += normal_lpdf(kappa_time[1]|drift,sigma_time);
  target += normal_lpdf(kappa_time[2:T]|drift+kappa_time[1:(T- 1)],sigma_time);    // Random walk with drift prior
  
  sum(kappa_time) ~ normal(0,0.001*T); //soft sum to zero constraint
  
  //Age Effect
  target += normal_lpdf(beta_age[1]|0,sigma_age); // First Effects seperate prior (see e.g. Fosse(2018))
  target += normal_lpdf(beta_age[2]|0,sigma_age);
  //RW(2) Prior
  target+= normal_lpdf(beta_age[3:A]| 2*beta_age[2:(A-1)]-beta_age[1:(A-2)], sigma_age);

  //Cohort Effect (IID)
  target+=std_normal_lpdf(var_cohort); //IID standard normal prior (for non centered)
  
  sum(beta_age) ~ normal(0,0.001*A); //soft sum to zero constraint
  
  
  //Region Effect (SAR Model)
  target += normal_lagsar_lpdf(phi_region | rep_vector(0.0, R), sigma_region, xi, W, eigenWsar);
  
  // soft sum-to-zero constraint on phi_region) to make mean identifiable
  sum(phi_region) ~ normal(0, 0.001 * R);  // equivalent to mean(phi_region) ~ normal(0,0.001)
  
  target += poisson_log_lpmf(y|log_E + InterceptVek + mu + eps*sigma_eps);
  
} generated quantities { // posterior
  
  vector[N] lambdahat  = exp(log_E + InterceptVek + //Log plus Intercept
                         kappa_time[TInd]+beta_age[AInd]+phi_region[RInd]+ gamma_cohort[CInd]+ //mu
                         eps_new*sigma_eps); // Error Term
                         
  vector [N] log_like_y; //= poisson_lpmf(y|lambda); // create Vector to store log likelihood (for waic (see Vethari, Gabry (2019)))
  
  vector[TFor] kappa_time_pred;
  vector[TFor] gamma_cohort_pred;
  vector[C+TFor] gamma_cohort_final;
  
  vector[L] mufor; // predicted death rates
  int pos_f = 1;
  int pos_L1 =1; 
  int <lower=0> K; // Kohort Index for Loop
  
  
  kappa_time_pred[1] = drift+kappa_time[T]+sigma_time * normal_rng(0,1);
  
  // FC of New Cohort Index
  //gamma_cohort_pred[1] = 2*gamma_cohort[C-1]-gamma_cohort[C-2]+sigma_cohort*normal_rng(0,1);
  
  // Check if Number of Years to Forecast is greater than 1
   if(TFor > 1){
  //   // RW(2) FC Cohort Index
  //   gamma_cohort_pred[2] = 2*gamma_cohort_pred[1]-gamma_cohort[C-1]+sigma_cohort*normal_rng(0,1);
  //   for(k in 3:TFor){
  //     gamma_cohort_pred[k]= 2*gamma_cohort_pred[k-1]-gamma_cohort_pred[k-2]+sigma_cohort*normal_rng(0,1);
  //   }
  //   //RW (1) - Drift Cohort Index
    for (t in 2:TFor) kappa_time_pred[t] = drift+kappa_time_pred[t - 1] + sigma_time * normal_rng(0,1);
 }
  for(t in 1:TFor){
    gamma_cohort_pred[t] = sigma_cohort*normal_rng(0,1);
  }
  // Append Forecasted Cohort Index on Existing Index
  gamma_cohort_final = append_row(gamma_cohort, gamma_cohort_pred);
  
  // Forecast Log Mortality Rates (attention may be slow due to if statement in loop)
  for (t in 1:TFor) for (r in 1:R) for (a in 1:A)  {
    if(a != 1 ){
      K = M*((A-1)-(a-1))+(T+t);  
    } else {
      K = M*((A-1)-(a))+(T+t); 
    }
    mufor[pos_f] =  Intercept+beta_age[a]+kappa_time_pred[t]+
                    phi_region[r]+
                    gamma_cohort_final[K]+
                    sigma_eps*normal_rng(0,1);
    pos_f += 1;
  }
   // Log Likelihood
  for (t in 1:T) for (r in 1:R) for (a in 1:A){
    log_like_y[pos_L1] = poisson_log_lpmf(y[pos_L1]| log_E[pos_L1] + Intercept +
                                                 kappa_time[t]+beta_age[a]+phi_region[r]+
                                                 gamma_cohort[CInd[pos_L1]]+
                                                 eps[pos_L1]*sigma_eps);
    pos_L1 += 1;
  }

}


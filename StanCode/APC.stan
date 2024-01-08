// APC Model
data{
  int <lower=1> T; //Time Index
  int <lower=1> A; // Age Index
  int <lower=1> R; // Region Index (not used for Estimation only indexing!!!)
  int <lower=0> y[T*A*R]; // count outcomes
  vector<lower=0>[T*A*R] E;// exposure
  int <lower=1> M; // How much wider is Age interval to Time 
  
  // Indicies for Data 
  int TInd[T*A*R];  // Index of Time in Data
  int AInd[T*A*R];  // Age Index in Data
  int ReInd[T*A*R]; // Index of Random Effects
  int CInd[T*A*R]; // Index for Cohort Effects
  
  int <lower=0> TFor; //Number of years to forecast 
  
} transformed data {
  int N = T*A*R; //Total Number of Observations
  vector[N] log_E = log(E); // log of exposure
  
  int<lower=1> C;// Cohort Index
  int<lower = 1> L; // size of prediction vector
  
  C = (M* (A-2)) + T; // max Cohort Index
  L=TFor*A*R; // Length of Forecasting Vector
  

}
parameters{
  real Intercept; // intercept
  real drift; // drift Parameter
  
  vector[A] alpha_age; // RW Age Coefficients
  
  vector[C] gamma_cohort;// Cohort Index
  
  vector[T] kappa_time; // RW Time Coefficients (only T-1 since Corner constraint is applied)

  //Variance Parameters
  real<lower=0> sigma_time;    // Variance of Time
  real<lower=0> sigma_age;     // Variance of Age
  real<lower=0> sigma_eps;     //Variance of epsilon
  real<lower=0> sigma_cohort; //Variance of Cohort Index
  
  vector[N] eps;     // Vector of overdispersion effects
  vector[N] eps_new; // new vector of overdisp effects for PPC
} 
transformed parameters {
  vector[N] InterceptVek = rep_vector(Intercept,N); // Vector of mu's
} 
model {
  vector[N] mu; // Vector of random effects
  mu[ReInd] = kappa_time[TInd]+alpha_age[AInd]+
              gamma_cohort[CInd];
  
  //Varianzen
  target += student_t_lpdf(sigma_age|5,0,1);
  target += student_t_lpdf(sigma_time|5,0,1);
  target += student_t_lpdf(sigma_eps|5,0,1);
  target += student_t_lpdf(sigma_cohort|5,0,1);

  
  //Other Effects
  target += normal_lpdf(Intercept|-5,5);
  target += normal_lpdf(drift|0,2);
  target += normal_lpdf(eps|0,1);
  
  target += std_normal_lpdf(eps_new); //
  
  // Random Walk with Drift Prior 
  target += normal_lpdf(kappa_time[1]|drift,sigma_time);
  target += normal_lpdf(kappa_time[2:T]|drift+kappa_time[1:(T- 1)],sigma_time);  // Random walk with drift prior
  
  sum(kappa_time)~normal(0,0.001*T);
  
  //Age Effect
  target += normal_lpdf(alpha_age[1:2]|0,sigma_age); // First Effects seperate prior 
  
  //RW(2) Prior
  target+= normal_lpdf(alpha_age[3:A]| 2*alpha_age[2:(A-1)]-alpha_age[1:(A-2)], sigma_age);
  sum(beta_age) ~ normal(0,0.001*A); //soft sum to zero constraint
  
  //Cohort Effect (IID Prior)
  target+=normal_lpdf(gamma_cohort|0,sigma_cohort); //IID normal prior
  
  sum(gamma_cohort)~ normal(0,0.001*C); // Soft sum to zero constraint
  
  
  target += poisson_log_lpmf(y|log_E + InterceptVek + mu + eps*sigma_eps);
  
} generated quantities { // posterior
  vector[N] lambdahat = exp(log_E + InterceptVek + //Log plus Intercept
                         kappa_time[TInd]+alpha_age[AInd]+ gamma_cohort[CInd]+ //mu
                         eps_new*sigma_eps); // Error Term (new random draws)
  
  
  // Insample Rates
  vector[N] MHat = exp(InterceptVek + //Log plus Intercept
                      kappa_time[TInd]+alpha_age[AInd]+ gamma_cohort[CInd]+ //mu
                       eps*sigma_eps); 
                       
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
    for (t in 2:TFor) kappa_time_pred[t] = drift+kappa_time_pred[t - 1] + sigma_time * normal_rng(0,1);
 }
  for(t in 1:TFor){
    gamma_cohort_pred[t] = sigma_cohort*normal_rng(0,1);
  }
 // Append Forecasted Cohort Index on Existing Index
  gamma_cohort_final = append_row(gamma_cohort, gamma_cohort_pred);
  
  // Forecast Log Mortality Rates
  for (t in 1:TFor) for (r in 1:R) for (a in 1:A)  {
    if(a != 1 ){
      K = M*((A-1)-(a-1))+(T+t);  //Cohort Index
    } else {
      K = M*((A-1)-(a))+(T+t); //Cohort Index 
    }
    mufor[pos_f] =  Intercept+kappa_time_pred[t]+alpha_age[a]+
                    gamma_cohort_final[K]+
                    sigma_eps*normal_rng(0,1);
    pos_f += 1;
  }
  
}

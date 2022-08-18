//RH_BYM2 Model
data{
  int <lower=1> T; //Time Index
  int <lower=1> A; // Age Index
  int <lower=1> R; //Region Index
  int <lower=0> y[T*A*R]; // count outcomes
  vector<lower=0>[T*A*R] E;// exposure
  int <lower=1> M; // How much wider is Age interval to Time 
  
  int<lower=0> N_edges;
  int<lower=1, upper=R> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=R> node2[N_edges];  // and node1[i] < node2[i]
  
  int<lower = 1> TFor;                  // number of forecast years
  
  real<lower=0> scaling_factor; // scales the Variances of the spatial effects
  
} transformed data {
  int N = T*A*R; //Total Number of Observations
  
  vector[N] log_E = log(E); // log of exposure
  
  int<lower=1> C;// Maximum Cohort Index Value
  int<lower = 1> Pred; // size of prediction vector
  int<lower=1> k[N]; // Kohort Index Vector
  int pos_k=1; // Helper for Cohort Index
  // 
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
  // 
  // 
  vector[C] gamma_cohort;// Cohort Index
  
  vector[T] kappa_time; //  for corner constraint
  real drift; // drift Parameter
  
  //Variances 
  real<lower=0> sigma_time;    // Variance of Time
  real<lower=0> sigma_eps;     //Variance of epsilon
  real<lower=0> sigma_BYM2;    //overall Standard deviation in BYM2 Model
  real<lower=0> sigma_cohort; //
  real<lower=0> sigma_alpha; //
  
  
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  
  
  //Spatial Model (BYM2)
  vector[R] v;       // heterogeneous effects (non spatial)
  vector[R] u;         // spatial effects
  
  
  vector[N] eps; // Vector of overdispersion effects
} 
transformed parameters {
  vector[R] convolved_re; // BYM2 Region Effect
  vector[R] phi_region; 
  convolved_re =  sqrt(1 - rho) * v + sqrt(rho / scaling_factor) * u;
  phi_region = convolved_re*sigma_BYM2;
} 
model {
  vector[N] mu; // Vector of random effects
  int pos = 1;
  for(t in 1:T) for(r in 1:R) for (a in 1:A){
    mu[pos]=alpha_age[a]+
    beta_age1[a]*kappa_time[t]+beta_age2[a]*gamma_cohort[k[pos]]+
    phi_region[r];
    pos += 1; //pos = pos + 1 
   }
  
  //Varianzen
  target +=student_t_lpdf(sigma_time|5,0,1);
  target +=student_t_lpdf(sigma_eps|5,0,1);
  target +=student_t_lpdf(sigma_BYM2|5,0,1);
  target +=student_t_lpdf(sigma_cohort|5,0,1);
  target +=student_t_lpdf(sigma_alpha|5,0,1);
  
  //Other Effects
  target += normal_lpdf(drift|0,2); // prior Intercept (medium vague)
  target += normal_lpdf(eps|0,1); //Epsilon Effects
  
   //Prior for rho
  target += beta_lpdf(rho|0.5,0.5);
  
  //TIME Effect
  target += normal_lpdf(kappa_time[1]|drift,sigma_time);
  target += normal_lpdf(kappa_time[2:T]|drift+kappa_time[1:(T- 1)],sigma_time);    // Random walk with drift prior
  
  sum(kappa_time) ~ normal(0,0.001*T) ; // soft sum to zero constraint on time 
  
  //Age Effect
  //target += normal_lpdf(alpha_age|0,sigma_alpha); // Prior on alpha_x
  target += normal_lpdf(alpha_age[1:2]|0,sigma_alpha);
  target += normal_lpdf(alpha_age[3:A]|2*alpha_age[2:(A-1)]-alpha_age[1:(A-2)],sigma_alpha); //RW2 Prior
  // 
  target += dirichlet_lpdf(beta_age1|rep_vector(1, A));// Prior on beta_x
  target += dirichlet_lpdf(beta_age2|rep_vector(1, A));// Prior on beta_x
  
  // Cohort Effect
  target += normal_lpdf(gamma_cohort|0,sigma_cohort);
  
  sum(gamma_cohort)~ normal(0,0.001*C); // Soft sum to zero constraint
  
  
  //Region Effect (Besag Model)
  target += -0.5 * dot_self(u[node1] - u[node2]);
  
  // soft sum-to-zero constraint on u)
  sum(u) ~ normal(0, 0.001 * R);  // equivalent to mean(u) ~ normal(0,0.001)
  
  target += normal_lpdf(v|0,1);
  
  target += poisson_log_lpmf(y| log_E + mu + eps*sigma_eps);
  
} generated quantities {
  
  // Quantities for InSample Fit
  vector[N] lambdahat;
  vector[N] MHat; // Vector of random effects
  int pos = 1;
  
  // Quantites for Forecast
  vector[TFor] kappa_time_pred;
  vector[TFor] gamma_cohort_pred;
  //Append Forecasted Cohort Index on Existing Index
  vector[C+TFor] gamma_cohort_final = append_row(gamma_cohort, gamma_cohort_pred);
  
  vector[Pred] mufor; // predicted death rates
  int pos_f = 1;
  int <lower=1> kFor; // Cohort Index for Forecast
  
  // get mu for insample fit
  for(t in 1:T) for(r in 1:R) for (a in 1:A){
    lambdahat[pos]=exp(alpha_age[a]+beta_age1[a]*kappa_time[t]+
                       beta_age2[a]*gamma_cohort[k[pos]]+
                        phi_region[r]+
                        log_E[pos]+normal_rng(0,1)*sigma_eps);
                        
                        
    MHat[pos] = exp(alpha_age[a]+beta_age1[a]*kappa_time[t]+
                       beta_age2[a]*gamma_cohort[k[pos]]+
                        phi_region[r]+
                        eps[pos]*sigma_eps);

    pos += 1; //pos = pos + 1 
   }
  
  kappa_time_pred[1] = drift+kappa_time[T]+sigma_time * normal_rng(0,1);
  
  // Check if Forecast Period is greater than 1
  if(TFor > 1){
  //RW (1) - Drift Cohort Index
   for (t in 2:TFor) kappa_time_pred[t] = drift+kappa_time_pred[t - 1] + sigma_time * normal_rng(0,1);
   }
   // FC of New Cohort Index
 for(t in 1:TFor){
     gamma_cohort_pred[t] = sigma_cohort*normal_rng(0,1);
   }
 //  
    // See that it has the same structure as loop in model block (else predictions are not comparable)
  for (t in 1:TFor) for (r in 1:R) for (a in 1:A)  {
     if(a==1){
      kFor = M*((A-1)-a)+(T+t);
    } else { 
      kFor = M*((A-1)-(a-1))+(T+t);
    }
    mufor[pos_f] = alpha_age[a] + 
                   beta_age1[a] * kappa_time_pred[t]+
                   beta_age2[a] * gamma_cohort_final[kFor]+
                   phi_region[r]+
                   sigma_eps*normal_rng(0,1);
    pos_f += 1;
  }
}


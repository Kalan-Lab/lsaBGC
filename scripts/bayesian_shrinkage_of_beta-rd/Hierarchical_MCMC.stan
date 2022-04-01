// STAN code for pooled hierarchical modeling of beta-RD statistic across GCFs.

data{
  int<lower = 0> n; // number of beta-rds
  int<lower = 0> J; // number of different GCFs

  vector[n] y; // beta-rd values
  int<lower = 1, upper = J> gcf_id[n]; // indicates to which GCF each beta-rd corresponds to 
  
  real<lower = 0> R1; // scale for prior on sigma
  real<lower = 0> R2; // scale for prior on tau
  real<lower = 0> sigma0; // prior for sigma of overline-mu
  real mu0;  
}
parameters{
  real<lower = 0> tau;
  real mu_bar;
  real<lower = 0> sigma;
  vector[n] aux;
}

transformed parameters{
  vector[J] mu; // vector of GCF-specific intercepts (on standardized scale)

  for (j in 1:J) {
      mu[j] = mu_bar + tau * aux[j];
  }
}

model{
  mu_bar ~ normal(mu0, sigma0);
  tau ~ student_t(7, 0.0, R2);
  sigma ~ student_t(7, 0.0, R1);

  for(j in 1:J){
    aux[j] ~ normal(0,1);
  }

  for(i in 1:n){
    y[i] ~ normal(mu[gcf_id[i]], sigma);
  }
}

generated quantities{
  real post_y[J];

  // sample the posterior
  for(j in 1:J){
      post_y[j] = normal_rng(mu[j], sigma);
  }
}

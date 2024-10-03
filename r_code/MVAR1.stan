data {
  int T;         //length of time series
  int p;         //number of variables
  matrix[p,T] Y; //matrix of observations; variables are rows; time is columns
//  real phi_mu;   //PHI prior distribution parameter
//  real phi_sigma; //PHI prior distribution parameter
}
parameters{
  matrix[p,p] PHI;     //dynamics matrix
  vector<lower=1E-15>[p] sigma;     //variances of stochastic forcing
  vector[p] init;      //mean of initial conditions as parameter vector
}
model{
  //for (icol in 1:p){
  //  PHI[,icol] ~ normal(phi_mu,phi_sigma);
  //}
  Y[,1] ~ normal(init, sigma);            //distribution of the initial conditions
  for(i in 2:T){
    Y[,i] ~ normal(PHI*Y[,i-1],sigma);  //conditional predictive distribution
  }
}

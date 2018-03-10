data {
  int Nobs; /* Number of observations */

  real logLobs[Nobs]; /* natural log of luminosity observed. */
  real sigma_logL[Nobs]; /* Uncertainty on each measurement. */

  real Tobs[Nobs]; /* Observed temperature. */
  real sigma_Tobs[Nobs]; /* Temperature uncertainty */
}

parameters {
  real<lower=0,upper=15> beta; /* Power-law slope L = T^beta */
  real<lower=0.8, upper=3>Ttrue[Nobs];  /* True temperature */
}

transformed parameters {
  real Ltrue[Nobs];
  real logLtrue[Nobs];

  for (i in 1:Nobs) {
    Ltrue[i] = Ttrue[i]^beta;
    logLtrue[i] = log(Ltrue[i]);
  }
}

model {
  /* Uniform prior on beta between limits. */
  /* Uniform prior on temperature between limits. */

  Tobs ~ normal(Ttrue, sigma_Tobs);
  logLobs ~ normal(logLtrue, sigma_logL);
}

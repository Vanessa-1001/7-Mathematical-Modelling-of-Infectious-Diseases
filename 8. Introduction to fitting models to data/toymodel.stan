data {
  int<lower=0> N;   // number of coin flips
  int<lower=0, upper=1> y[N]; // 1=heads, 0=tails
}

parameters {
  real<lower=0, upper=1> theta; //prob. heads
}

model {
  theta ~ uniform(0, 1);  // prior
  y ~ bernoulli(theta); //likelihood
}


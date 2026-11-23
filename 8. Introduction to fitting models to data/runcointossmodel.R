N<-20
y<-c(0,0,1,0,1,1,1,0,1,0,0,0,1,1,0,1,1,0,0,1)
coinflipdata<-list(N=N, y=y)
library(rstan)

# This translates model to C++ and compiles and runs it

fit1 <- stan(
  file = "toymodel.stan",  # Stan program
  data = coinflipdata,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  refresh = 0             # no progress shown
)


print(fit1, pars=c("theta"),
      probs=c(0.1, 0.5, 0.9), digits = 3)
plot(fit1)
traceplot(fit1)
check_hmc_diagnostics(fit1)


require(bayesplot)


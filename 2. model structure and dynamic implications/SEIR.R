library(deSolve)

SEIR_model<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+E+I+R)
    rand <- max(rnorm(1, I, 2), 0)
    P_aux <- S + E + rand + R
    lam<- max(rnorm(1,parameters['beta'], 5e-1), 0) * I/P

    
    dS <- mu*P-mu*S-lam*S   
    dE <- lam*S-gamma*E-mu*E
    dI <- gamma*E-tau*I-mu*I 
    dR <- tau*I-mu*R
    
    list(c(dS, dE, dI, dR))}) 
}


initP<-999
initE<-0
initI<-1 
initR<-0 
initS<-initP-initI-initR-initE
istate <- c(S = initS, E = initE, I = initI, R = initR)

parameters <- c(
  mu=(1/(50*52*7)),    
  beta=.35,              
  gamma=1/5,          
  tau=1/10            
)

R0<-parameters["beta"]*parameters["gamma"]/
  ((parameters["gamma"]+parameters["mu"])*(parameters["tau"]+parameters["mu"]))
R0
time_start <- 0
time_stop <- 365 
deltat<- 0.01
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SEIR_model, parms = parameters) 

Flu2009<-as.data.frame(out)
plot(out)
# plot(Flu2009)

N<-rowSums(Flu2009[,2:5])
plot(N)

library(EpiEstim)
Flu2009$incidence<-parameters["gamma"]*Flu2009$E
ww<-which(Flu2009$time %% 1 ==0.0)

res_parametric_si <- estimate_R(abs(Flu2009$incidence[ww]), 
    method="parametric_si", 
    config = make_config(list( mean_si = 2.5, std_si = 1.5)) ) 
plot(res_parametric_si, legend = FALSE)
# epiestime estimate
res_parametric_si$R$`Mean(R)`[1:100]
# compare with  Rt in the simulation
Rt <- R0*Flu2009$S/(Flu2009$S+Flu2009$E+Flu2009$I+Flu2009$R)

plot(Rt[ww], type="l", col="red", lwd=2, ylim=c(0,3), ylab="Rt", xlab="Time (days)")
lines(res_parametric_si$R$`Mean(R)`[1:100], col="blue", lwd=2)



SEIR_model_4n<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+E+I1+I2+I3+I4+R)
    n<-4
    lam<- rnorm(1,parameters['beta'],5e-3)*(I1+I2+I3+I4)/P
    dS <- mu*P-mu*S-lam*S   
    dE <- lam*S-gamma*E-mu*E
    dI1 <- gamma*E-tau*n*I1-mu*I1 
    dI2 <- tau*n*I1-mu*I2-tau*n*I2
    dI3 <- tau*n*I2-mu*I3-tau*n*I3
    dI4 <- tau*n*I3-mu*I4-tau*n*I4
    dR <- n*tau*I4-mu*R
    list(c(dS, dE, dI1, dI2, dI3, dI4, dR))}) 
}

initP<-999
initE<-0
initI1<-1 
initI2<-1 
initI3<-1 
initI4<-1 
initR<-0 
initS<-initP-initI1-initI2-initI3-initI4-initR-initE
istate4n <- c(S = initS, E = initE, I1 = initI1, I2 = initI2, I3 = initI3, I4 = initI4, R = initR)

out4n <- deSolve::ode(y = istate4n, times = tps, func = SEIR_model_4n, parms = parameters) 

Flu2009_4n<-as.data.frame(out4n)
plot(out4n)
plot(Flu2009_4n)


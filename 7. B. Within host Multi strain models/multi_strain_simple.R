library(deSolve)
setwd("C:/MMID2627/day_07/age-dependent models")

SIR_2strains<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),
  {
    P <- S+I1+I2+R1+R2+I12+I21+R12
    lam1 <- beta*I1/P
    lam2 <- beta*I2/P
    
    dS <- mu*P - mu*S- lam1*S - lam2*S + alpha1*R1 + alpha2*R2 
    dI1 <- lam1*S - tau1*I1 - mu*I1 
    dR1 <- tau1*I1 - mu*R1 - sigma*lam2*R1 - alpha1*R1 + alpha2*R12
    dI2 <- lam2*S - tau2*I2 - mu*I2 
    dR2 <- tau2*I2 - mu*R2 - sigma*lam1*R2 - alpha2*R2 + alpha1*R12
    dI12<- sigma*lam2*R1 - tau2*I12 - mu*I12
    dI21<- sigma*lam1*R2 - tau1*I21 - mu*I21
    dR12<- tau2*I12 + tau1*I21 - mu*R12 - alpha1*R12 - alpha2*R12
    
    list(c(dS, dI1, dR1, dI2, dR2, dI12, dI21, R12))
  }
  )
}


initP<-999 
initI1<-1
initR1<-0
initR2<-0
initI2<-1
initI12<-0
initI21<-0
initR12<-0
initS<-initP-initI1-initR1-initR2-initI2-initI12-initI21-initR12
istate <- c(S = initS, I1 = initI1, R1 = initR1, 
            I2 = initI2, R2 = initR2, I12 = initI12, I21=initI21,
            R12 = initR12)

parameters <- c(
  mu = (1/(50*52*7)),
  beta = 0.5,
  tau1 = 1/10,
  tau2 = 1/10,
  sigma= 0.2,
  alpha1= 1/300,
  alpha2= 1/365
)

time_start <- 0
time_stop <- 1000
deltat<- 0.01
tps <- seq(time_start , time_stop , by = deltat)


out <- deSolve::ode(y = istate, times = tps, func = SIR_2strains, parms = parameters) 
plot(out[,'I1']+out[,'I21'])
lines(out[,'I2']+out[,'I12'], col = "cyan")


istate['I1']<-10
out2 <- deSolve::ode(y = istate, times = tps, func = SIR_2strains, parms = parameters) 
plot(out2[,'I1']+out2[,'I21'],ylim=c(0,350))
lines(out2[,'I2']+out2[,'I12'], col = "cyan")


parameters['tau2']<-1/20
out3 <- deSolve::ode(y = istate, times = tps, func = SIR_2strains, parms = parameters) 
lines(out3[,'I2']+out3[,'I12'], col = "magenta")


parameters['sigma']<-0
out4 <- deSolve::ode(y = istate, times = tps, func = SIR_2strains, parms = parameters) 
plot(out4[,'I1']+out4[,'I21'])
lines(out4[,'I2']+out4[,'I12'], col = "purple")

plot(out2[,'I1']+out2[,'I21'])
lines(out4[,'I1']+out4[,'I21'],col="red")



########################################
# Generalised SIR model with X strains and cross-immunity
#########################################

A <- 5 # number of strains

# characteristics of each strain
beta <- seq(0.15, 2.5, length.out = A) # transmission rate
p <- rep(0.2, A) # probability of infection given a contact
omega <- rep(0/(2*365), A) # rate of loss of immunity = 1/(average duration of immunity)
rho <- rep(0.0, A) # relative infectiousness of incubation phase
nui <- rep(1/5, A) # rate of recovery = 1/(average duration of symptomatic infection phase)
nus <- rep(1/7, A) # rate if recovery for non-fatal severe
ratem <- rep(1/7, A) # 1/ time to death for fatal
rhos <- rep(0.1, A) # relative infectiousness*contacts of severe ill patients
gamma <- rep(1/4.5, A) # rate of movement from incubation to infectious stage = 1/(average incubation period)
isr <- rep(0.05, A) # prob of severe disease if case for each strain
sfr <- rep(0.005, A) # prob of death if severe for each strain
mort <- 1/(75*365) # natural mortality rate

# cross immunity matrix between strains as a diagonal matrix with 1s on the diagonal and values between 0 and 1 off-diagonal
for(i in 1:A){
  for(j in 1:A){
    if(i==j){
      cross_immunity[i,j] <- 1
    } else {
      cross_immunity[i,j] <- 2*exp(-abs(i-j))
    }
  }
}
image(cross_immunity, main = "Cross-immunity matrix between strains", xlab = "Strain", ylab = "Strain", col = heat.colors(10), axes = FALSE)

startdate <- as.Date("2020-01-31")
stopdate <- as.Date("2026-12-31")
day_start <- as.numeric(startdate - startdate)
day_stop <- as.numeric(stopdate - startdate)
times <- seq(day_start, day_stop)

##  MODEL PARAMETERS
parameters <- list(
  beta = beta,              # transmission rate
  p = p,                    # probabilty of infection given a contact 
  rho = rho,                # relative infectiousness of incubation phase
  omega = omega,            # rate of loss of immunity = 1/(average duration of immunity)
  gamma = gamma,            # rate of movement from incubation to infectious stage = 1/(average incubation period)
  nui = nui,                # rate of recovery = 1/(average duration of symptomatic infection phase)
  report = 1 / 8,           # proportion of all infections that are reported
  ratem = ratem,            # 1/ time to death for fatal
  mort = mort,              # natural mortality rate
  nus = nus,                # rate if recovery for non-fatal severe
  rhos = rhos,              # relative infectiousness*contacts of severe ill patients 
  isr = isr,                # prob of severe disease if case for each strain
  sfr = sfr,                # prob of death if severe for each strain
  cross_immunity = cross_immunity  # cross-immunity between strains
)


###########################################################################
# Define the indices for each variable
Sindex <- 1 : A
Eindex <- (A + 1) : (2 * A)
Iindex <- (2 * A + 1) : (3 * A)
Rindex <- (3 * A + 1) : (4 * A)
Hindex <- (4 * A + 1) : (5 * A)
Mindex <- (5 * A + 1) : (6 * A)
Cindex <- (6 * A + 1) : (7 * A)
CMindex <- (7 * A + 1) : (8 * A)

aci <- 5 # infecting strain index

###########################################################################
# MODEL INITIAL CONDITIONS
initI <- rep(2, A) # Infected and symptomatic
initE <- rep(1, A) # Incubating
# initE[aci] <- 1 # place random index case in E compartment
initR <- rep(0, A) # Immune
initH <- rep(0, A) # hospitalised 
initM <- rep(0, A) # died 
initC <- rep(0, A) # Cumulative cases (true)
initCM <- rep(0, A) # Cumulative deaths (true)
initP <- rep(1000, A)  # population size
initS <- initP - initE - initI -
          initR - initH - initM - initC - initCM # Susceptible (non-immune)

# initial conditions for the main solution vector
Y <- c(initS, initE, initI, initR,
       initH, initM, initC, initCM)


# set up a function to solve the equations
covid_strains <- function(t, Y, parameters) {
  with(as.list(c(Y, parameters)),
    {
      S <- Y[Sindex]
      E <- Y[Eindex]
      I <- Y[Iindex]
      R <- Y[Rindex]
      H <- Y[Hindex]
      M <- Y[Mindex]
      C <- Y[Cindex]
      CM <- Y[CMindex]
      # print(t)

      P <- (S + E + I + R + H + M)
      # print(dim(P))
      # print(S)
      birth <- mort * P
         
      lam <- p * beta * ((rho * E + I + rhos * H) / P)
      # print(dim(lam))
      dSdt <- - lam * S + omega * R  - mort * S + birth
      dEdt <- lam * S - gamma * E  - mort * E + lam %*% cross_immunity * R
      dIdt <- gamma * (1 - isr) * E - nui * I - mort * I
      dRdt <- nui * I - omega * R  - mort * R + nus * H -
            lam %*% cross_immunity * R
      dHdt <- gamma * isr * (1 - sfr) * E -
               mort * H - nus * H
      dMdt <- gamma * isr * sfr * E -  ratem * M - mort * M
         
      dCdt <- report * gamma * E
      dCMdt <- ratem * M + mort * M
         
      # return the rate of change
      list(c(dSdt, dEdt, dIdt, dRdt, dHdt, dMdt, dCdt, dCMdt))
    }
  )
}

# run the model
out <- ode(y = Y, times = times, func = covid_strains, parms = parameters)

# total population
pop <- out[, (Sindex + 1)] + out[, (Eindex + 1)] + out[, (Iindex + 1)] +
       out[, (Rindex + 1)] + out[, (Hindex + 1)] + out[, (Mindex + 1)]
tpop <- rowSums(pop)
time <- as.Date(out[, 1] + startdate)
plot(tpop)

# model dynamics - all compartments
plot(out[ ,"time"], rowSums(out[, (Sindex + 1)]), type = "l", col = "blue", ylim = c(0, 5000), xlab = "Time", ylab = "Population")
lines(out[ ,"time"], rowSums(out[, (Eindex + 1)]), col = "orange")
lines(out[ ,"time"], rowSums(out[, (Iindex + 1)]), col = "red")
lines(out[ ,"time"], rowSums(out[, (Rindex + 1)]), col = "green")
lines(out[ ,"time"], rowSums(out[, (Hindex + 1)]), col = "purple")
lines(out[ ,"time"], rowSums(out[, (Mindex + 1)]), col = "brown")
legend("topright", legend = c("Susceptible", "Exposed", "Infected", "Recovered", "Hospitalised", "Dead"), 
                      col = c("blue", "orange", "red", "green", "purple", "brown"), lty = 1, cex = 0.8)

# daily incidence
inc <- unlist(parameters["report"]) * unlist(parameters["gamma"]) * out[, (Eindex + 1)]
inc <- as.data.frame(inc)
inc_strains <- cbind(time, inc)
colnames(inc_strains) <- c("time", paste("strain", 1:A))
inc_strains_long <- reshape2::melt(inc_strains, id.vars = colnames(inc_strains)[1], variable.name = "variable", value.name = "value")
colnames(inc_strains_long) <- c("time", "strain", "incidence")

ggplot(inc_strains_long, aes(x = time, y = incidence, color = strain)) +
  geom_line() +
  labs(title = "Daily Incidence of Each Strain Over Time",
       x = "Time",
       y = "Incidence")
       
cum_inc <- colSums(unlist(parameters["report"]) * out[, (Cindex + 1)])  

# prevalence of each strain
prevalence_strains <- (out[, (Iindex + 1)]) # prevalence of infection
prevalence_strains_long <- reshape2::melt(prevalence_strains, id.vars = colnames(prevalence_strains)[1], variable.name = "variable", value.name = "value")
colnames(prevalence_strains_long) <- c("time", "strain", "prevalence")
prevalence_strains_long$strain <- as.factor(prevalence_strains_long$strain)
ggplot(prevalence_strains_long, aes(x = time, y = prevalence, color = strain)) +
  geom_line() +
  labs(title = "Prevalence of Each Strain Over Time",
       x = "Time",
       y = "Prevalence")

# morbidity and mortality
hospreq <- rowSums(out[, (Hindex + 1)]) # requirement for hospital beds
cmortality <- rowSums(out[, (CMindex + 1)]) # cumulative mortality
# plot hospital occupancy against the number of hospital beds available
plot(hospreq)

#### age groups
deaths <- (out[, CMindex + 1])
colnames(deaths) <- c("0-4", "5-9", "10-14", "15-19", "20-24",
            "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
            "55-59", "60-64", "65-69", "70-74", "75-79", "80-84",
            "85-89", "90-94", "95-99", "100+")
deaths_prop <- round(tail(deaths, 1), 0) / sum(round(tail(deaths, 1), 0))
barplot(deaths_prop)
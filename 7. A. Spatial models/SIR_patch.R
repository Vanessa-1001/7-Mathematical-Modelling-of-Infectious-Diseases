library(deSolve)

SIR_patch<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P1 <- S1 + I1 + R1
    P2 <- S2 + I2 + R2

    if (t >= tbit){beta1 <- beta1 * (1-tbeff)}

    lam1 <- beta1 * I1 / P1
    lam2 <- beta2 * I2 / P2

    f12 <- amp*cos(2*pi * t / 365) * f12 + f12
    f21 <- amp*cos(2*pi * t / 365) * f21 + f21
        
    dS1 <- mu*P1 - mu*S1 - lam1*S1 - f12*S1 + f21*S2
    dI1 <- lam1*S1 - tau*I1 - mu*I1 - f12*I1 + f21*I2
    dR1 <- tau*I1 - mu*R1 - f12*R1 + f21*R2
    dS2 <- mu*P2 - mu*S2 - lam2*S2 + f12*S1 - f21*S2
    dI2 <- lam2*S2 - tau*I2 - mu*I2 + f12*I1 - f21*I2
    dR2 <- tau*I2 - mu*R2 + f12*R1 - f21*R2

    N <- P1 + P2
    count <- 1
    list(c(dS1, dI1, dR1, dS2, dI2, dR2, N, count))})
}


initP1 <- 1000
initI1 <- 1
initR1 <- 0
initS1 <- initP1 - initI1 - initR1

initP2 <- 1000
initI2 <- 0
initR2 <- 0
initS2 <- initP2 - initI2 - initR2

istate <- c(S1 = initS1, I1 = initI1, R1 = initR1, 
            S2 = initS2, I2 = initI2, R2 = initR2, 
            N = initS1 + initI1 + initR1 + initS2 + initI2 + initR2,
            cout = 0)

parameters <- c(
  mu = (1 / (50 * 52 * 7)),
  beta1 = 0.75,
  beta2 = 0.45,
  tau = 1 / 10,
  f12 = 0.001,
  f21 = 0.001,
  tbit = 1e100,
  amp = 0.5
)

time_start <- 0
time_stop <- 1000
deltat <- 0.01
tps <- seq(time_start , time_stop , by = deltat)
tt <- seq(1, length(tps), by = 1/deltat)

out <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
plot(out)
ww <- which(out[, 1] %% 1 == 0)
sum(out[ww, 2:4]) + sum(out[ww, 5:7])
tail(out[, 8], 1)
tail(out[, 9], 1)


parameters["f12"] <- 0.0
out_no_migration <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
# plot(out_no_migration)

parameters["f12"] <- 0.001
parameters["tbit"] <- 10
parameters["tbeff"] <- 0.2
out_intervention <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
par(mfrow=c(2,1))
plot(out[,3], type = "l", col = "orange3", main = "Pop1")
lines(out_intervention[,3], type = "l", col = "blue")
legend("topright", legend = c("Control", "Intervention"), lwd = 2, col = c("orange3", "blue"))
plot(out[,6], type = "l", col = "orange3", main = "Pop2")
lines(out_intervention[,6], type = "l", col = "blue")


parameters["tbit"] <- 1e100
parameters["f12"] <- 0.003
parameters["f21"] <- 0.001
parameters["beta2"] <- 0.25
out_high_migr <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
parameters["tbit"] <- 10
out_intervention_high_migr <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
par(mfrow=c(2,1))
plot(out_high_migr[,3], type = "l", col = "orange3", main = "Pop1")
lines(out_intervention_high_migr[,3], type = "l", col = "blue")
legend("topright", legend = c("Control", "Intervention"), lwd = 2, col = c("orange3", "blue"))
plot(out_high_migr[,6], type = "l", col = "orange3", main = "Pop2")
lines(out_intervention_high_migr[,6], type = "l", col = "blue")




########################################################################################

##    Seasonal migration effect on SIR patch model

########################################################################################

SIR_patch_seasonal <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P1 <- S1 + I1 + R1
    P2 <- S2 + I2 + R2

    if (t >= tbit){beta1 <- beta1 * (1-tbeff)}

    lam1 <- beta1 * I1 / P1
    lam2 <- beta2 * I2 / P2

    f12s <- f12 + amp*cos(2*pi*(t-0)/365)
    f21s <- f21 + amp*cos(2*pi*(t-0)/365)
        
    dS1 <- mu*P1 - mu*S1 - lam1*S1 - f12s*S1 + f21s*S2
    dI1 <- lam1*S1 - tau*I1 - mu*I1 - f12s*I1 + f21s*I2
    dR1 <- tau*I1 - mu*R1 - f12s*R1 + f21s*R2
    dS2 <- mu*P2 - mu*S2 - lam2*S2 + f12s*S1 - f21s*S2
    dI2 <- lam2*S2 - tau*I2 - mu*I2 + f12s*I1 - f21s*I2
    dR2 <- tau*I2 - mu*R2 + f12s*R1 - f21s*R2

    N <- P1 + P2
    count <- 1
    list(c(dS1, dI1, dR1, dS2, dI2, dR2, N, count))})
}

parameters["tbit"] <- 1e100
parameters["f12"] <- 0.001
parameters["f21"] <- 0.001
parameters["beta2"] <- 0.0
parameters["amp"] <- 0.00095
out_base_migr_seas <- deSolve::ode(y = istate, times = tps, func = SIR_patch_seasonal, parms = parameters) 
par(mfrow=c(2,1))
plot(out_base_migr_seas[,2], type = "l", col = "#cd0029", main = "Seasonal Migration Effect on S1 and S2", ylab = "Number Susceptible", xlab = "Time")
lines(out_base_migr_seas[,5], type = "l", col = "#00ff7b")
legend("topright", legend = c("S1", "S2"), lwd = 2, col = c("#cd0029",  "#00ff44"))
parameters["amp"] <- 0.0
out_base_migr_no_seas <- deSolve::ode(y = istate, times = tps, func = SIR_patch_seasonal, parms = parameters) 
plot(out_base_migr_no_seas[,2], type = "l", col = "#cd0029", main = "Non-Seasonal Migration", ylab = "Number Susceptible", xlab = "Time")
lines(out_base_migr_no_seas[,5], type = "l", col = "#00ff7b")
legend("topright", legend = c("S1", "S2"), lwd = 2, col = c("#cd0029",  "#00ff44"))



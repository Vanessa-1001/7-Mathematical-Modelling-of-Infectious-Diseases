##################################################################
###       Spatial SEIRS MODEL WITH population movement         ###
##################################################################
require("plot.matrix")
require("deSolve")
require("igraph")
library("sf")
library("spdep")
setwd("C:/MMID2627/day_07/spatial models")

##########################################################################
# MOVEMENT DATA
movement <- as.matrix(read.csv("movement_spatial.csv", header = FALSE))
movement <- movement * 0.0002
# movement <- 0
plot(movement)
##########################################################################
# Population structure
popstruc <- as.data.frame(read.csv("pop_spatial.csv", header = TRUE))
pop_size <- popstruc$pop
pop_x <- popstruc$longitude
pop_y <- popstruc$latitude
areas <- length(pop_size)

###########################################################################
# SPATIAL RISK  DATA
risk <- popstruc$incidence*100

###########################################################################
# VISUALISE NETWORK
X <- as.matrix(movement)
g2 <- graph.adjacency(X, mode = "undirected",
    weighted = TRUE, diag = FALSE)
# Assign attributes to the graph
g2$name <- "SEIR spatial model"
# Assign attributes to the graph's vertices
V(g2)$color <- ifelse(risk > 5, "red", "green")
V(g2)$x <- pop_x
V(g2)$y <- pop_y
V(g2)$size <- pop_size / 5e4
plot(g2, edge.width = E(g2)$weight * 8)

startdate <- as.Date("2022-01-01")
stopdate <- as.Date("2022-12-31")
day_start <- as.numeric(startdate-startdate)
day_stop <- as.numeric(stopdate-startdate)
times <- seq(day_start, day_stop, by = 0.1)

#MODEL PARAMETERS
parameters <- c(
  p = 0.047,                # probability of infection given a contact 
  rho = 0.2,                # relative infectiousness of incubation phase
  omega = (1 / (100 * 365)),# rate of loss of immunity = 1/(average duration of immunity) 
  gamma = 1 / 4.5,          # rate of movement from incubation to infectious stage = 1/(average incubation period)
  nui = 1 / 5,              # rate of recovery = 1/(average duration of symptomatic infection phase)
  report = 1 / 8,           # proportion of all infections that are reported
  ratem = 1 / 7,            # 1/ time to death for fatal
  nus = 1 / 7,              # rate if recovery for non-fatal severe
  rhos = 0.1,               # relative infectiousness*contacts of severely ill patients 
  csr = 0.1,                # case severity ratio
  sfr = 0.2,                # severe case fatality ratio
  mort = (1 / (80 * 365))   # mortality rate
)


###########################################################################
# Define the indices for each variable
Sindex <- 1 : areas
Eindex <- (areas + 1) : (2 * areas)
Iindex <- (2 * areas + 1) : (3 * areas)
Rindex <- (3 * areas + 1) : (4 * areas)
Hindex <- (4 * areas + 1) : (5 * areas)
Mindex <- (5 * areas + 1) : (6 * areas)
Cindex <- (6 * areas + 1) : (7 * areas)
CMindex <- (7 * areas + 1) : (8 * areas)

###########################################################################
# MODEL INITIAL CONDITIONS
initI <- round(0.01 * popstruc$pop) # Infected and symptomatic
initE <- 0 * popstruc$pop # Incubating
initR <- 0 * popstruc$pop # Immune
initH <- 0 * popstruc$pop # hospitalised
initM <- 0 * popstruc$pop # died
initC <- 0 * popstruc$pop # Cumulative cases (true)
initCM <- 0 * popstruc$pop # Cumulative deaths (true)
initS <- popstruc$pop - initE - initI -
          initR - initH - initM # Susceptible (non-immune)

# initial conditions for the main solution vector
Y <- c(initS, initI, initE, initR,
       initH, initM, initC, initCM)

# set up a function to solve the equations
spatial_model <- function(t, Y, parameters) {
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
      contacts <- risk 
    #   print(contacts)
      P <- (S + E + I + R + H + M)
      lam <- p * contacts * ((rho * E + I + rhos * H) / P)
    #   print(lam)
    #   print(movement)

      dSdt <- - S * lam + omega * R - mort * S  + mort * P -
        movement %*% S + t(movement) %*% S
      dEdt <- S * lam - gamma * E - mort * E -
        movement %*% E + t(movement) %*% E
      dIdt <- gamma * (1 - csr) * E - nui * I - mort * I -
        movement %*% I + t(movement) %*% I
      dRdt <- nui * I - omega * R - mort * R + nus * H -
        movement %*% R + t(movement) %*% R

      dHdt <- - mort * H + gamma * csr * (1 - sfr) * E - nus * H
      dMdt <-  gamma * csr * sfr * E - ratem * M - mort * M
         
      dCdt <- report * gamma * E
      dCMdt <- ratem * M + mort * M
         
      # return the rate of change
      list(c(dSdt, dEdt, dIdt, dRdt, dHdt, dMdt, dCdt, dCMdt))
    }
  )
}

# run the model
out <- ode(y = Y, times = times, func = spatial_model, parms = parameters)
# total population
pop <- out[, (Sindex + 1)] + out[, (Eindex + 1)] + out[, (Iindex + 1)] +
       out[, (Rindex + 1)] + out[, (Hindex + 1)] + out[, (Mindex + 1)]
tpop <- rowSums(pop)
plot(tpop)


# Calculate Global Moran's I Over Time
# Define symmetric binary spatial weights matrix W (neighbors = 1, self = 0)
W <- movement
W <- matrix(movement, nrow = areas, byrow = TRUE)

# Extract only the Infected columns 
df <- as.data.frame(out)
I_mat <- as.matrix(df[, Iindex + 1])

nbrs <- neighborhood(g2, order = 1)
class(nbrs) <- "nb"
lw <- nb2listw(nbrs, style = "W", zero.policy = TRUE)
# ---------------------------------------------------------
#  Calculate Global Moran's I & P-Values Over Time
# ---------------------------------------------------------
moran_results <- apply(I_mat, 1, function(x) {
  # Guard against near-zero variance (all patches identical)
  if (max(abs(x - mean(x))) < 1e-6) {
    return(c(moran_I = NA, p_value = NA))
  }
  # Analytical Moran's I test using spdep library
  test <- moran.test(x, listw = lw, zero.policy = TRUE)
  return(c(
    moran_I = as.numeric(test$estimate["Moran I statistic"]),
    p_value = test$p.value
  ))
})
# Combine results with time
moran_df <- cbind(time = df$time, as.data.frame(t(moran_results)))

# ---------------------------------------------------------
#  Plot Spatial Autocorrelation & Significance Over Time
# ---------------------------------------------------------
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
# Plot global Moran's I trajectory
plot(moran_df$time, moran_df$moran_I, type = "l", col = "darkblue", lwd = 2,
     ylim = c(-0.2, 1), xlab = "Time (days)", ylab = "Global Moran's I",
     main = "Spatial Autocorrelation Across Network")
abline(h = 0, lty = 3, col = "gray")

# Plot p-value trajectory showing when spatial clustering is statistically significant
plot(moran_df$time, moran_df$p_value, type = "l", col = "firebrick", lwd = 2,
     ylim = c(0, 1), xlab = "Time (days)", ylab = "p-value",
     main = "Statistical Significance of Spatial Clustering")
abline(h = 0.05, lty = 2, col = "darkgreen")
legend("topright", legend = "p = 0.05 threshold", lty = 2, col = "darkgreen")


# daily incidence
inc <- parameters["report"] * parameters["gamma"] * out[, (Eindex + 1)]
dailyinc <- rowSums(inc)
hospreq <- rowSums(out[, (Hindex + 1)]) # requirement for hospital beds
cmortality <- rowSums(out[, (CMindex + 1)]) # cumulative mortality
cum_inc_cases <- colSums(unlist(parameters["report"]) * out[, (Cindex + 1)])  
# plot hospital occupancy against the number of hospital beds available
plot(hospreq)

#### spatial groups
deaths <- (out[, CMindex + 1])
colnames(deaths) <- c("Area1", "Area2", "Area3", "Area4", "Area5",
            "Area6", "Area7", "Area8", "Area9", "Area10", "Area11",
            "Area12", "Area13", "Area14", "Area15", "Area16")
deaths_prop <- round(tail(deaths, 1), 0) / sum(round(tail(deaths, 1), 0))
barplot(deaths_prop)

movement <- movement * 0
# run the model with no movement
out_no_move <- ode(y = Y, times = times, func = spatial_model, parms = parameters)
inc_no_move <- parameters["report"] * parameters["gamma"] * out_no_move[, (Eindex + 1)]
# total incidence by area
barplot(colSums(inc_no_move), col = "firebrick")
barplot(colSums(inc), col = "dodgerblue2")
## plot incidence over time
plot(rowSums(inc))   ## aggregate incidence in all areas
lines(rowSums(inc_no_move), col = "green")


df2 <- as.data.frame(out_no_move)
# Combine results with time
I_mat2 <- as.matrix(df2[, Iindex + 1])
# ---------------------------------------------------------
#  Calculate Global Moran's I & P-Values Over Time
# ---------------------------------------------------------
moran_results2 <- apply(I_mat2, 1, function(x) {
  # Guard against near-zero variance 
  if (max(abs(x - mean(x))) < 1e-6) {
    return(c(moran_I = NA, p_value = NA))
  }
  
  # Analytical Moran's I test using spdep
  test <- moran.test(x, listw = lw, zero.policy = TRUE)
  return(c(
    moran_I = as.numeric(test$estimate["Moran I statistic"]),
    p_value = test$p.value
  ))
})
# Combine results with time
moran_df2 <- cbind(time = df2$time, as.data.frame(t(moran_results)))
# ---------------------------------------------------------
#  Plot Spatial Autocorrelation & Significance Over Time
# ---------------------------------------------------------
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
# Plot Global Moran's I trajectory
plot(moran_df2$time, moran_df2$moran_I, type = "l", col = "darkblue", lwd = 2,
     ylim = c(-0.2, 1), xlab = "Time (days)", ylab = "Global Moran's I",
     main = "Spatial Autocorrelation Across Network")
abline(h = 0, lty = 3, col = "gray")

# Plot p-value trajectory (showing when spatial clustering is statistically significant)
plot(moran_df2$time, moran_df2$p_value, type = "l", col = "firebrick", lwd = 2,
     ylim = c(0, 1), xlab = "Time (days)", ylab = "p-value",
     main = "Statistical Significance of Spatial Clustering")
abline(h = 0.05, lty = 2, col = "darkgreen")
legend("topright", legend = "p = 0.05 threshold", lty = 2, col = "darkgreen")


inc <- as.data.frame(inc)
inc_strains <- cbind(out[,"time"], inc)
colnames(inc_strains) <- c("time", paste("area", 1:areas))
inc_strains_long <- reshape2::melt(inc_strains, id.vars = colnames(inc_strains)[1], variable.name = "variable", value.name = "value")
colnames(inc_strains_long) <- c("time", "area", "incidence")

library(ggplot2)
ggplot(inc_strains_long, aes(x = time, y = incidence, color = area)) +
  geom_line() +
  labs(title = "Daily Incidence of Each Strain Over Time",
       x = "Time",
       y = "Incidence")


inc_no_move <- as.data.frame(inc_no_move)
inc_strains2 <- cbind(out[,"time"], inc)
colnames(inc_strains2) <- c("time", paste("area", 1:areas))
inc_strains_long2 <- reshape2::melt(inc_strains2, id.vars = colnames(inc_strains2)[1], variable.name = "variable", value.name = "value")
colnames(inc_strains_long2) <- c("time", "area", "incidence")

ggplot(inc_strains_long2, aes(x = time, y = incidence, color = area)) +
  geom_line() +
  labs(title = "Daily Incidence of Each Strain Over Time",
       x = "Time",
       y = "Incidence")




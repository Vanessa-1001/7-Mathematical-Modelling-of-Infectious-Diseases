install.packages(c("deSolve","ggplot2","tidyr","patchwork"))

library(deSolve)
library(ggplot2)
library(tidyr)
library(patchwork)

# -------------------------
# Model equations
# -------------------------
two_strain_tiv <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dT  <- - beta1 * (T/N) * V1 - beta2 * (T/N) * V2 + alpha
    dI1 <-   beta1 * (T/N) * V1 - delta * I1
    dI2 <-   beta2 * (T/N) * V2 - delta * I2
    dV1 <- (1 - m) * p * I1 + m * p * I2 - c * V1 - kappa * E * V1
    dV2 <- (1 - m) * p * I2 + m * p * I1 - c * V2 - kappa*kappa_esc * E * V2
    dE  <-  rho * (I1 + I2) - gamma * E
    
    list(c(dT, dI1, dI2, dV1, dV2, dE))
  })
}

# -------------------------
# Parameters (arbitrary, not based on a particular pathogen though roughly similar to influenza)
# -------------------------
parms <- c(
  alpha = 0.001, ## Generation rate of new target cells
  beta1 = 0.25, ## Infection rate of Virus 1 per virion per day
  beta2 = 0.25, ## Infection rate of Virus 2 per virion per day
  delta = 5.0, ## Death rate of infected cells
  p = 1000, ## virions produced per infected cell per day
  c = 2.0, ## Non-immune-mediated clearance rate of free virus per day
  m = 0, ## Probability of mutation between virus strains per virion produced
  kappa = 0.01, ## Immune-mediated clearance rate of Virus 1 per immune effector per day
  kappa_esc = 1, ## Factor reduction in immune clearance rate for escape variant (Virus 2)
  rho = 0.01, ## Rate of immune effector generation per infected cell per day
  gamma = 0.1, ## Immunity waning rate per day
  N = 1e8 ## Total number of target cells
)

# -------------------------
# Initial conditions & times
# -------------------------
state0 <- c(
  T  = 1e8, ## Initial number of target cells
  I1 = 0.0, ## Initial number of infected cells with Virus 1
  I2 = 0.0, ## Initial number of infected cells with Virus 2
  V1 = 1, ## Initial number of Virus 1 virions
  V2 = 0, ## Initial number of Virus 2 virions
  E  = 0 ## Initial number of immune effectors e.g., cells, antibodies
)

times <- seq(0, 15, length.out = 1000)  # days

# -------------------------
# Solve
# -------------------------
out <- ode(y = state0, times = times, func = two_strain_tiv, parms = parms, atol = 1e-8, rtol = 1e-8)
out <- as.data.frame(out)

# -------------------------
# Tidy for plotting
# -------------------------
df <- pivot_longer(out, cols = -time, names_to = "state", values_to = "value")

# add tiny floor to avoid -Inf on log scale
tiny <- 1
df$value_safe <- pmax(df$value, tiny)

# -------------------------
# Plot (log y)
# -------------------------
p_log <- ggplot(df %>% dplyr::filter(state %in% c("V1","V2","E"))) +
  geom_line(aes(x = time, y = value_safe, color = state),size = 1) +
  geom_hline(yintercept=1000,linetype="dashed") +
  scale_y_log10() +
  labs(x = "time (days)", y = "counts (log scale, arbitrary units)",
       color = "state") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p_natural <- ggplot(df %>% dplyr::filter(state %in% c("V1","V2","E")), aes(x = time, y = value_safe, color = state)) +
  geom_line(size = 1) +
  geom_hline(yintercept=1000,linetype="dashed") +
  labs(x = "time (days)", y = "counts (natural scale, arbitrary units)",
       color = "state") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_log/p_natural

## 1. Baseline, single virus type. What do you think is causing the infection to end?
## 2. What if there is existing immunity? Set E0 to a higher value, try 1e4
## 3. Introduce second virus, change the initial seeding, who wins? What about varying the transmission rates? Think about the seeding scale -- log10. Try doubling the transmission rate of one of the viruses.
## 4. Second virus arises through mutations -- set m to a small value like 1e-4. What happens now?
## 5. Give second virus a fitness advantage -- try setting it to 90% immune escape (kappa_escape=0.1)
## 6. Return the initial immune seeding to the higher value and seed with the second virus with substantial immune escape. (E0=1e4, kappa_escape=0.1) What happens now?

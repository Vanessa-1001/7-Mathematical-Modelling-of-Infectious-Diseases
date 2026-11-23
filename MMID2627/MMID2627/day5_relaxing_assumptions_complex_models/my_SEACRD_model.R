

# Load the packages required
library(deSolve,
        ggplot2,
        tidyr,
        dplyr)

# ==========================================================================
# STEP 1: BASELINE COMPARTMENTS (INITIAL CONDITIONS)
# Define starting population sizes for all compartments
initP  <- 1000   # Total population size (N)
initE  <- 0      # Exposed (first-time, unprotected)
initA  <- 0      # Asymptomatic infectious
initC  <- 1      # Clinical infectious (index case)
initH  <- 0      # Hospitalised (added for Step 2 severity tracking)
initR  <- 0      # Recovered
initER <- 0      # Exposed reinfected (added for Step 4 reinfection tracking)
initD  <- 0      # Deceased

# Step 7: Prophylaxis Intervention compartments
initPr <- 0      # Susceptible individuals currently protected by prophylaxis
initEP <- 0      # Exposed individuals who got infected while on prophylaxis

# Initial Susceptible population (remaining individuals after subtracting all other states)
initS <- initP - initE - initA - initC - initH - initR - initER - initD - initPr - initEP

# Combine into the initial state vector for deSolve
istate <- c(
  S  = initS,
  E  = initE,
  A  = initA,
  C  = initC,
  H  = initH,
  R  = initR,
  ER = initER,
  D  = initD,
  Pr = initPr,
  EP = initEP
)

# ==========================================================================
# PARAMETERS (LABELLED BY STEP)

parameters <- c(
  # --- Step 1: Baseline Compartments (Core transition rates) ---
  beta  = 0.5,        # Baseline transmission rate (per day)
  nu    = 1/5,        # Incubation rate (1 / 5 days)
  tau_A = 1/7,        # Recovery rate for Asymptomatic cases (1 / 7 days)
  tau_C = 1/5,        # Exit rate for Clinical cases (1 / 5 days)
  tau_h = 1/10,       # Exit rate for Hospitalised cases (1 / 10 days)
  
  # --- Step 2: Disease Severity Parameters ---
  gamma = 0.4,        # ICR: Infection-to-Clinical Ratio (40% of first infections become clinical)
  theta = 0.3,        # CHR: Clinical Hospitalisation Ratio (30% of clinical cases are hospitalised)
  omega = 0.02,       # HFR: Hospitalisation Fatality Ratio (2% of hospitalised cases die)
  
  # --- Step 3: Differential Infectivity ---
  rel_inf_A = 0.2,    # Asymptomatic cases have 20% infectivity of clinical cases (E has 0%)
  
  # --- Step 4: Reinfection ---
  sigma = 1.0,        # Relative susceptibility to reinfection (1 = same force of infection as S)
  # Note: Reinfected clinical probability (gamma / 2 = 0.2) is calculated in the ODE
  
  # --- Step 5: Temporary Immunity ---
  alpha = 1/180,      # Rate of waning immunity from R back to S (1 / 180 days)
  
  # --- Step 7: Prophylaxis Intervention ---
  phi_p = 0.01,       # Rate at which susceptible individuals receive prophylaxis (set to 0 for baseline)
  rho   = 1/30        # Waning rate of prophylaxis protection (1 / 30 days = 1 month)
)


# ==========================================================================
# ODE MODEL FUNCTION (STEP 1 TO STEP 7)
SEACRD <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Step 1: Baseline Compartments - Total living population (excludes Deceased D)
    N <- S + E + A + C + H + R + ER + Pr + EP
    
    # Step 3: Differential Infectivity - E, ER, and EP do not transmit; A transmits at 20% of C
    lambda <- beta * (C + rel_inf_A * A) / N
    
    # --- Susceptible (S) ---
    # Step 1: Loses infected individuals (-lambda * S)
    # Step 5: Temporary Immunity - Gains individuals whose natural immunity wanes (+alpha * R)
    # Step 7: Prophylaxis Intervention - Loses people starting prophylaxis (-phi_p * S) and regains them when the 1-month drug wears off (+rho * Pr)
    dS <- -lambda * S + alpha * R - phi_p * S + rho * Pr
    
    # --- Exposed First-Time (E) ---
    # Step 1: Gains newly infected susceptibles (lambda * S) and progresses at incubation rate (-nu * E)
    dE <- lambda * S - nu * E
    
    # --- Step 7: Prophylaxis Intervention Compartments (Pr & EP) ---
    # Pr: Susceptibles protected by prophylaxis; drug wanes after 1 month (-rho * Pr), and they can still be exposed (-lambda * Pr)
    dPr <- phi_p * S - rho * Pr - lambda * Pr
    
    # EP: Exposed while on prophylaxis; progresses at incubation rate (-nu * EP)
    dEP <- lambda * Pr - nu * EP
    
    # --- Asymptomatic Infectious (A) ---
    # Step 1 & Step 2: Gains non-clinical fraction (1 - gamma) from E, and recovers at rate tau_A
    # Step 4: Reinfection - Gains higher asymptomatic fraction (1 - gamma / 2) from ER
    # Step 7: Prophylaxis Intervention - Gains 100% of EP because the drug prevents clinical symptoms
    dA <- nu * (1 - gamma) * E + nu * (1 - gamma / 2) * ER + nu * EP - tau_A * A
    
    # --- Clinical Infectious (C) ---
    # Step 1 & Step 2: Gains clinical fraction (gamma) from E, and exits at rate tau_C
    # Step 4: Reinfection - Gains halved clinical fraction (gamma / 2) from ER
    # Step 7: Prophylaxis Intervention - Receives 0 inflow from EP because prophylaxis prevents symptoms
    dC <- nu * gamma * E + nu * (gamma / 2) * ER - tau_C * C
    
    # --- Hospitalised (H) ---
    # Step 2: Disease Severity Parameters - Gains CHR fraction (theta) of exiting clinical cases, exits at rate tau_h
    dH <- tau_C * theta * C - tau_h * H
    
    # --- Recovered (R) ---
    # Step 1 & Step 2: Gains recoveries from A, non-hospitalised C (1 - theta), and surviving H (1 - omega)
    # Step 4: Reinfection - Loses reinfected individuals (-sigma * lambda * R)
    # Step 5: Temporary Immunity - Loses individuals whose immunity wanes back to S (-alpha * R)
    dR <- tau_A * A + tau_C * (1 - theta) * C + tau_h * (1 - omega) * H - sigma * lambda * R - alpha * R
    
    # --- Exposed Reinfected (ER) ---
    # Step 4: Reinfection - Gains reinfected individuals from R (sigma * lambda * R) and incubates at rate nu
    dER <- sigma * lambda * R - nu * ER
    
    # --- Deceased (D) ---
    # Step 2: Disease Severity Parameters - Gains HFR fraction (omega) of exiting hospitalised cases
    dD <- tau_h * omega * H
    
    
    # ==========================================================================
    # STEP 6: EPIDEMIOLOGICAL OUTPUTS (DAILY INCIDENCE TRACKING)
    case_incidence   <- nu * E + nu * ER + nu * EP   # All daily new infections (A + C)
    hospitalisations <- tau_C * theta * C            # Daily new hospital admissions
    deaths           <- tau_h * omega * H            # Daily new deaths
    
    # Return derivatives and tracked Step 6 outputs
    list(
      c(dS, dE, dA, dC, dH, dR, dER, dD, dPr, dEP),
      Case_Incidence   = case_incidence,
      Hospitalisations = hospitalisations,
      Deaths           = deaths
    )
  })
}

# ==============================================================================
# STEP 6 & STEP 7: RUN SIMULATIONS AND PLOT EPIDEMIOLOGICAL OUTPUTS
times <- seq(0, 365, by = 1)

# Run model WITH Step 7 Prophylaxis Intervention (phi_p = 0.01)
out_with_drug <- as.data.frame(ode(y = istate, times = times, func = SEACRD, parms = parameters))

# Run model WITHOUT Step 7 Prophylaxis Intervention (phi_p = 0) for baseline comparison
params_no_drug <- parameters
params_no_drug["phi_p"] <- 0
out_no_drug <- as.data.frame(ode(y = istate, times = times, func = SEACRD, parms = params_no_drug))

# Label and combine both scenarios
out_with_drug$Scenario <- "With Step 7 Prophylaxis"
out_no_drug$Scenario   <- "No Prophylaxis (Steps 1-6)"
combined_out <- bind_rows(out_no_drug, out_with_drug)

# Reshape data to plot Step 6 Epidemiological Outputs: Case Incidence, Hospitalisations, Deaths
plot_df <- combined_out %>%
  select(time, Scenario, Case_Incidence, Hospitalisations, Deaths) %>%
  pivot_longer(cols = c(Case_Incidence, Hospitalisations, Deaths),
               names_to = "Metric", values_to = "Daily_Count")

ggplot(plot_df, aes(x = time, y = Daily_Count, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 13) +
  labs(
    title = "SEACRD Model: Step 6 Epidemiological Outputs",
    subtitle = "Comparing Baseline (Steps 1-6) vs. Step 7 Prophylaxis Intervention",
    x = "Time (Days)",
    y = "Daily Count"
  )

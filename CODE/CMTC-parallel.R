run_one_sim = function(params) {
  IC <- c(S = (params$Pop-params$I_0), I = params$I_0)
  params_local <- c(gamma = params$gamma, beta = params$beta)
  reactions <- list(
    # propensity function effects name for reaction
    reaction("beta*S*I", c(S=-1,I=+1), "new_infection"),
    reaction("gamma*I", c(S=+1,I=-1), "recovery")
  )
  set.seed(NULL)
  sol <- ssa(
    initial_state = IC,
    reactions = reactions,
    params = params_local,
    method = ssa_exact(),
    final_time = params$t_f,
    log_firings = TRUE    # This way we keep track of events
  )
  # Interpolate result (just I will do)
  wanted_t = seq(from = 0, to = params$t_f, by = 0.01)
  sol$interp_I = approx(x = sol$time, y = sol$state[,"I"], 
                        xout = wanted_t)
  names(sol$interp_I) = c("time", "I")
  # Return result
  return(sol)
}

# Set parameters, initial conditions, etc.
gamma = 1/3
R0 = 1.5
# R0=beta/gamma*S0, so beta=R0*gamma/S0
beta = as.numeric(R0*gamma/IC["S"])
t_f = 100
params <- list(gamma = gamma, beta = beta,
               Pop = 1000, I_0 = 2, R0 = R0,
               t_f = 100, nb_sims = 50)

library(future.apply)
# Prepare parallel processing environment
if (availableCores() >= 64) {
  # We're on the cluster.. Be reasonable to avoid overheating
  nb_cores_to_use = 50
} else {
  # We're on something smaller, reserve 2 cores to avoid sluggish
  nb_cores_to_use = availableCores()-2
}
# Now choose future plan
plan(multisession, workers = nb_cores_to_use)
# To run sequentially (e.g., to debug), run sequentially
# plan(sequential)

# Run the parallel code
SIMS <- 
  future_lapply(
    1:params$nb_sims, 
    function(x) 
      run_one_sim(params),
    future.seed = TRUE
  )

# Find max y value for plot
y_max = max(unlist(lapply(SIMS, function(x) max(x$interp_I$I))), 
            na.rm = TRUE)
# Now plot
plot(SIMS[[1]]$interp_I$time,
     SIMS[[1]]$interp_I$I,
     type = "l", lwd = 0.5,
     xlab = "Time (days)",
     ylab = "Number infectious",
     ylim = c(0, y_max),
     main = paste("CTMC with R0 =", params$R0))
for (i in 2:length(SIMS)) {
  lines(SIMS[[i]]$interp_I$time,
        SIMS[[i]]$interp_I$I,
        type = "l", lwd = 0.5)
}


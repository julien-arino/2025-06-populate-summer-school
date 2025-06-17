## ----set-options,echo=FALSE,warning=FALSE,message=FALSE-----------------------
# Load required libraries
required_packages = c("deSolve",
                      "dplyr", 
                      "ggplot2", 
                      "knitr", 
                      "latex2exp",
                      "lattice",
                      "magick",
                      "readr", 
                      "tidyr",
                      "viridis")

for (p in required_packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, dependencies = TRUE)
    require(p, character.only = TRUE)
  }
}
# Knitr options
opts_chunk$set(echo = TRUE, 
               warning = FALSE, 
               message = FALSE, 
               dev = c("pdf", "png"),
               fig.width = 6, 
               fig.height = 4, 
               fig.path = "FIGS/course-01-",
               fig.keep = "high",
               fig.show = "hide")
knitr::knit_hooks$set(crop = knitr::hook_pdfcrop)
options(knitr.table.format = "latex")
# Date for front title page (if needed)
yyyy = strsplit(as.character(Sys.Date()), "-")[[1]][1]


## ----set-slide-background,echo=FALSE,results='asis'---------------------------
# Are we plotting for a dark background?
plot_blackBG = FALSE
if (plot_blackBG) {
  bg_color = "black"
  fg_color = "white"
  input_setup = "\\input{slides-setup-blackBG.tex}"
} else {
  bg_color = "white"
  fg_color = "black"
  input_setup = "\\input{slides-setup-whiteBG.tex}"
}
cat(input_setup)


## ----KMK_SI_plane,echo=FALSE--------------------------------------------------
# Plot the trajectories of the KMK SIR system in the SI-plane
# In this plot, the total population N is normalised to 1
gamma = 1/5
beta = 0.5
# Store the values of S and I in lists
S = list()
I = list()
i = 1
# Solutions starting on the S axis
for (S0 in seq(0.6, 1, 0.1)) {
  I0 = 0
  S[[i]] = seq(0.001, S0, 0.001)
  I[[i]] = S0+I0-S[[i]]+gamma/beta*log(S[[i]]/S0)
  i = i+1
}
# Solutions starting on the S+I=1 line
for (S0 in seq(1, 0.05, -0.1)) {
  I0 = 1-S0
  S[[i]] = seq(0.001, S0, 0.001)
  I[[i]] = S0+I0-S[[i]]+gamma/beta*log(S[[i]]/S0)
  i = i+1
}
# S+I=1 line
S[[i]] = seq(0, 1, 0.001)
I[[i]] = 1-S[[i]]
if (plot_blackBG) {
  par(bg = 'black', fg = 'white') # set background to black, foreground white
  colour = "white"
} else {
  colour = "black"
}
for (i in 1:length(S)) {
  if (i == 1) {
    plot(S[[i]], I[[i]],
         type = "l", lwd = 3,
         col = ifelse((I[[i]][length(I[[i]])] < max(I[[i]])), "red", colour),
         xlim = c(0,1), ylim = c(0,1),
         xaxs = "i", yaxs = "i",
         xlab = "S", ylab = "I",
         col.axis = colour, cex.axis = 1,
         col.lab = colour, cex.lab = 1,
         bty = "n")
    points(S[[i]][length(S[[i]])], I[[i]][length(I[[i]])],
          pch = 19, cex = 2, 
          col = ifelse((I[[i]][length(I[[i]])] < max(I[[i]])), "red", colour))
  } else if (i<length(S)) {
    lines(S[[i]], I[[i]], 
          col = ifelse((I[[i]][length(I[[i]])] < max(I[[i]])), "red", colour),
          lwd = 3)
    points(S[[i]][length(S[[i]])], I[[i]][length(I[[i]])],
           pch = 19, cex = 2, 
           col = ifelse((I[[i]][length(I[[i]])] < max(I[[i]])), "red", colour))
  } else {
    lines(S[[i]], I[[i]])
  }
}


## ----eval=TRUE----------------------------------------------------------------
rhs_SIR_KMK <- function(t, x, p) {
  with(as.list(c(x, p)), {
    dS = - beta * S * I
    dI = beta * S * I - gamma * I
    dR = gamma * I
    return(list(c(dS, dI, dR)))
  })
}
# Initial condition for S (to compute R_0)
S0 = 1000
gamma = 1/14
# Set beta so that R_0 = 1.5
beta = 1.5 * gamma / S0 
params = list(gamma = gamma, beta = beta)
IC = c(S = S0, I = 1, R = 0)
times = seq(0, 365, 1)
sol_KMK <- ode(IC, times, rhs_SIR_KMK, params)  


## ----KMK_R0eq1dot5------------------------------------------------------------
plot(sol_KMK[, "time"], sol_KMK[, "I"], 
     type = "l", lwd = 2,
     main = TeX("Kermack-McKendrick SIR, $R_0=1.5$"),
     xlab = "Time (days)", ylab = "Prevalence")


## ----eval=TRUE----------------------------------------------------------------
final_size_eq = function(S_inf, S0 = 999, I0 = 1, R_0 = 2.5) {
  OUT = S0*(log(S0)-log(S_inf)) - (S0+I0-S_inf)*R_0
  return(OUT)
}


## -----------------------------------------------------------------------------
uniroot(f = final_size_eq, interval = c(0.05, 999))


## -----------------------------------------------------------------------------
final_size = function(L) {
  with(as.list(L), {
  S_inf = uniroot(f = function(x) 
    final_size_eq(S_inf = x, 
                  S0 = S0, I0 = I0, 
                  R_0 = R_0),
    interval = c(0.05, S0))
  return(S_inf$root)
  })
}


## ----KMK_final_size_0p8-------------------------------------------------------
N0 = 1000
I0 = 1
S0 = N0-I0
R_0 = 0.8
S = seq(0.1, S0, by = 0.1)
fs = final_size_eq(S, S0 = S0, I0 = I0, R_0 = R_0)
S_inf = uniroot(f = function(x) final_size_eq(S_inf = x, 
                                              S0 = S0, I0 = I0, 
                                              R_0 = R_0),
                interval = c(0.05, S0))
plot(S, fs, type = "l", ylab = "Value of equation (10)")
abline(h = 0)
points(x = S_inf$root, y = 0, pch = 19)
text(x = S_inf$root, y = 0, labels = "S_inf", adj = c(-0.25,-1))


## ----KMK_final_size_2p5,echo=FALSE--------------------------------------------
N0 = 1000
I0 = 1
S0 = N0-I0
R_0 = 2.5
S = seq(0.1, S0, by = 0.1)
fs = final_size_eq(S, S0 = S0, I0 = I0, R_0 = R_0)
S_inf = uniroot(f = function(x) final_size_eq(S_inf = x, 
                                              S0 = S0, I0 = I0, 
                                              R_0 = R_0),
                interval = c(0.05, S0))
plot(S, fs, type = "l", ylab = "Value of equation (10)")
abline(h = 0)
points(x = S_inf$root, y = 0, pch = 19)
text(x = S_inf$root, y = 0, labels = "S_inf", adj = c(-0.25,-1))    


## ----KMK_attack_rate----------------------------------------------------------
values = expand.grid(
  R_0 = seq(0.01, 3, by = 0.01),
  I0 = seq(1, 100, 1)
)
values$S0 = N0-values$I0
L = split(values, 1:nrow(values))
values$S_inf = sapply(X = L, FUN = final_size)
values$final_size = values$S0-values$S_inf+values$I0
values$attack_rate = (values$final_size / N0)*100

p = levelplot(attack_rate ~ R_0*I0, data = values, 
              xlab = TeX("$R_0$"), ylab = "I(0)",
              col.regions = viridis(100))
print(p)


## ----SIRS_one_sim_prevalence--------------------------------------------------
library(deSolve)
rhs_SIRS <- function(t, x, p) {
  with(as.list(c(x, p)), {
    dS = b + nu * R - d * S - beta * S * I
    dI = beta * S * I - (d + gamma) * I
    dR = gamma * I - (d + nu) * R
    return(list(c(dS, dI, dR)))
  })
}
# Initial conditions
N0 = 1000
I0 = 1
R0 = 0
IC = c(S = N0-(I0+R0), I = I0, R = R0)
# "Known" parametres
d = 1/(80*365.25)
b = N0 * d
gamma = 1/14
nu = 1/365.25
# Set beta s.t. R_0 = 1.5
R_0 = 1.5
beta = R_0 * (d + gamma) / (N0-I0-R0)
params = list(b = b, d = d, gamma = gamma, beta = beta, nu = nu)
times = seq(0, 500, 1)
# Call the numerical integrator
sol_SIRS <- ode(y = IC, times = times, func = rhs_SIRS, 
                parms = params, method = "ode45")
# Plot the result
plot(sol_SIRS[,"time"], sol_SIRS[,"I"], 
     type = "l", lwd = 2,
     xlab = "Time (days)", ylab = "Prevalence")


## ----SIRS_3_sims_prevalence,echo=FALSE----------------------------------------
# Compute the EPs
valeur_PE = function(params) {
  with(as.list(c(params)), {
    OUT = list()
    if (R_0<1) {
      OUT$S_EP = Pop
      OUT$I_EP = 0
      OUT$col = "dodgerblue4"
    } else {
      OUT$S_EP = 1/R_0*Pop
      OUT$I_EP = (1-1/R_0)*(d+nu)/(d+nu+gamma)*Pop
      OUT$col = "darkorange4"
    }
    return(OUT)
  })
}
# RHS function set in previous chunk

# Put the parameters in a list
# "Known" parametres
params = list()
params$Pop = N0
params$d = 1/(80 * 365.25)
params$b = params$Pop * params$d
params$gamma = 1/14
params$nu = 1/365.25
params$t_f = 1200
params$I_0 = I0
# Note that we did not set R_0 or beta. This is done in a loop

# IC. "Static" part (N0, I0, R0) of IC are set in previous chunk
IC = c(S = N0-(I0+R0), I = I0, R = R0)

# Times at which the solution will be returned.
tspan = seq(from = 0, to = params$t_f, by = 0.1)

# Now simulate the ODE. Loop on several values of R_0
R_0 = c(0.8, 1.5, 2.5)
# Save results in a list together with EP values
sol_ODE = list()
EP = list()
# Now loop on R_0
for (r_0 in R_0) {
  # Name for list entry
  entry_name = sprintf("$R_0$=%1.1f",r_0)
  # Keep the current value of R_0 to compute EPs
  params$R_0 = r_0
  # R0=(beta/(d+gamma)) => beta=R0*(d+gamma)
  params$beta = r_0 * (params$d+params$gamma) / (N0-I0-R0)
  # Call numerical integrator
  sol_ODE[[entry_name]] = ode(y = IC,
                              func = rhs_SIRS,
                              times = tspan,
                              parms = params)
  EP[[entry_name]] = valeur_PE(params)
  EP[[entry_name]]$lty = which(r_0 == R_0)
}

# Get maximum value of I across all simulations for plot. Note the use of lapply.
max_I = max(unlist(lapply(sol_ODE, function(x) max(x[,"I"]))))

# Plot
plot(sol_ODE[[1]][,"time"], sol_ODE[[1]][,"I"],
     ylim = c(0, max_I),
     type = "l", lwd = 5, col = EP[[1]]$col, lty = EP[[1]]$lty,
     xlab = "Time (days)", ylab = "Prevalence")
points(x = params$t_f, y = EP[[1]]$I_EP, 
       col = EP[[1]]$col, pch = 19, cex = 2)
for (i in 2:length(sol_ODE)) {
  lines(sol_ODE[[i]][,"time"], sol_ODE[[i]][,"I"],
        type = "l", lwd = 5, col = EP[[i]]$col, lty = EP[[i]]$lty)
  points(x = params$t_f, y = EP[[i]]$I_EP, 
         col = EP[[i]]$col, pch = 19, cex = 2)
}
legend("topright", legend = TeX(names(EP)), cex = 0.8,
       col = unlist(lapply(EP, function(x) x$col)),
       lty = unlist(lapply(EP, function(x) x$lty)),
       lwd = c(3,3,3))


## ----SIRS_bifurcation_R0,echo=FALSE-------------------------------------------
# Values of the EPs
value_EPs = function(R_0, N) {
  EP_I = ifelse(R_0 < 1, 0, (1-1/R_0)*N)
  return(EP_I)
}

R_0 = seq(0.5, 5, by = 0.01)
EP_I = value_EPs(R_0, N = 1000)
# We also show the DFE when R_0>1, so prepare this
R_0_geq_1 = R_0[which(R_0>=1)]
DFE = rep(0, length(R_0_geq_1))

plot(R_0, EP_I,
     type = "l", lwd = 3,
     xlab = TeX("$R_0$"),
     las = 1,
     ylab = "Prevalence at equilibrium")
lines(R_0_geq_1, DFE,
      type = "l", lwd = 3,
      lty = 2)
legend("topleft", legend = c("LAS EP", "Unstable EP"),
       lty = c(1, 2), lwd = c(2,2),
       bty = "n")


## ----convert-Rnw-to-R,warning=FALSE,message=FALSE-----------------------------
# From https://stackoverflow.com/questions/36868287/purl-within-knit-duplicate-label-error
rmd_chunks_to_r_temp <- function(file){
  callr::r(function(file, temp){
    out_file = sprintf("../CODE/%s", gsub(".Rnw", ".R", file))
    knitr::purl(file, output = out_file, documentation = 1)
  }, args = list(file))
}
rmd_chunks_to_r_temp("course-01-introduction-math-epi.Rnw")


## ----eval=FALSE---------------------------------------------------------------
# pp = ggplot(...)
# print(pp)


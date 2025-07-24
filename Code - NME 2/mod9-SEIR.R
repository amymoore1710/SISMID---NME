##
## NME @ SIDMID Extending EpiModel Module: An SEIR Epidemic
##

library(EpiModel)

rm(list = ls())


# EpiModel Model Extensions -----------------------------------------------


  #Running a simple SIR model 
# Network & epidemic models
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
param <- param.net(inf.prob = 0.3, rec.rate = 0.1)
init <- init.net(i.num = 10, r.num = 0)
control <- control.net(type = "SIR", nsteps = 25, nsims = 1, verbose = FALSE)
mod1 <- netsim(est1, param, init, control)
mod1

# Print function arguments
args(control.net)

# Inspect function contents
View(infection.net)

View(recovery.net)

# General function design
fx <- function(dat, at) {

  ## function processes that update dat

  return(dat)
}

# See infection module

# See progress module


# Network Model -----------------------------------------------------------

# Initialize the network
nw <- network_initialize(500)

# Define the formation model: edges + degree terms
formation = ~edges + degree(1)

# Input the appropriate target statistics for each term
target.stats <- c(500, 180)

# Parameterize the dissolution model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 25)
coef.diss

# Fit the model
est <- netest(nw, formation, target.stats, coef.diss)

# Model diagnostics
dx <- netdx(est, nsims = 10, ncores = 5, nsteps = 500,
            nwstats.formula = ~edges + degree(0:7),
            keep.tedgelist = TRUE)
print(dx)
plot(dx)


# Epidemic Model Parameterization  ----------------------------------------

# Parameterize model
param <- param.net(inf.prob = 0.5, act.rate = 2,
                   ei.rate = 0.01, ir.rate = 0.01)

# Initial conditions
init <- init.net(i.num = 10)

# Set control settings & source modules from SEIR-fx script
source("mod9-SEIR-fx.R")
control <- control.net(type = NULL,
                       nsteps = 500,
                       nsims = 1,
                       ncores = 1,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       resimulate.network = FALSE)

# Epidemic simulation
sim <- netsim(est, param, init, control)

# Simulate model again with 10 sims
source("mod9-SEIR-fx.R")
control <- control.net(type = NULL,
                       nsteps = 500,
                       nsims = 10,
                       ncores = 5,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       resimulate.network = FALSE)
sim <- netsim(est, param, init, control)

# Print sim
print(sim)

# Default plots
par(mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(sim)

plot(sim, y = c("se.flow", "ei.flow", "ir.flow"), legend = TRUE)

# Save sim to data frame
df <- as.data.frame(sim)
df[which(df$time == 100), ]

# Plot transmission matrices as phylograms
tm1 <- get_transmat(sim)
plot(tm1)


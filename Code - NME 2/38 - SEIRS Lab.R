
  # 7.23.2025
  # Amy Moore
  # SISMID - NME II
  # 38 - SEIRS Lab

  #Load Packages
library(EpiModel)


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
  ###Adding rs.rate -> the rate at which recovered individuals transition back to being susceptible
  ###rs.rate = 1/100 ---> it takes 100 time steps for immunity to wane
param <- param.net(inf.prob = 0.5, act.rate = 2,
                   ei.rate = 0.01, ir.rate = 0.01, 
                   rs.rate = 0.01)

# Initial conditions
init <- init.net(i.num = 10)

# Set control settings & source modules from SEIR-fx script
source("38 - SEIRS Lab Helper Functions.R")
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
source("38 - SEIRS Lab Helper Functions.R")
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

plot(sim, y = c("se.flow", "ei.flow", "ir.flow", "rs.flow"), legend = TRUE)



    ##### Model #2 - Shorter Immunity After Recovery #####

# Parameterize model
###Adding rs.rate -> the rate at which recovered individuals transition back to being susceptible
###rs.rate = 1/10 ---> it takes 10 time steps for immunity to wane
param2 <- param.net(inf.prob = 0.5, act.rate = 2,
                   ei.rate = 0.01, ir.rate = 0.01, 
                   rs.rate = 0.1)

# Simulate model again with 10 sims
sim2 <- netsim(est, param2, init, control)

# Print sim
print(sim2)

# Default plots
par(mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(sim2)

plot(sim2, y = c("se.flow", "ei.flow", "ir.flow", "rs.flow"), legend = TRUE)




##### Model #3 - Moderate Immunity After Recovery #####

# Parameterize model
###Adding rs.rate -> the rate at which recovered individuals transition back to being susceptible
###rs.rate = 1/50 ---> it takes 10 time steps for immunity to wane
param3 <- param.net(inf.prob = 0.5, act.rate = 2,
                    ei.rate = 0.01, ir.rate = 0.01, 
                    rs.rate = 0.02)

# Simulate model again with 10 sims
sim3 <- netsim(est, param3, init, control)

# Print sim
print(sim3)

# Default plots
par(mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(sim3)

plot(sim3, y = c("se.flow", "ei.flow", "ir.flow", "rs.flow"), legend = TRUE)



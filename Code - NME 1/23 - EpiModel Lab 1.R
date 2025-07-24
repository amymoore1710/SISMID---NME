
  # 7.22.2025
  # Amy Moore
  # SISMID - NME
  # 23 - EpiModel Lab 1

#Load Packages
library(EpiModel)


#Initialize Network
nw <- network_initialize(n = 500)



#Model Parameterization
formation <- ~edges + degrange(from = 2)
# edges - counts the number of edges
# degrange - hard-codes that no nodes have degree higher than 2 (implied from 2 to inf)
  ### this should prevent concurrency since nodes can no longer have degree 2 or higher
  ### could also keep concurrent with a target stat of 0

#List of target stats (# edges, # nodes degree 2+)
target.stats <- c(175, 0)



#Dissolution Model
#offset for the # of edges plus a term for the duration of the edge
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
#coefficient is log(Odds) of dissolution
coef.diss



#Estimate Network Model with tergm
#nw - the initialized empty network
#formation - formation formula (model parameterization)
#target.stats - data (i.e. sufficient stats)
#coef.diss - dissolution coefs object 
est <- netest(nw, formation, target.stats, coef.diss)



#Model Diagnostics

# Non-parallel version
# dx <- netdx(est, nsims = 10, nsteps = 1000,
#             nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent,
#             keep.tedgelist = TRUE)

#Parallelize this diagnostic procedure
numCores <- parallel::detectCores() #check how many cores
#run the same code in parallel (leave one core open)
dx <- netdx(est, nsims = 10, nsteps = 1000, ncores = numCores - 1,
            nwstats.formula = ~edges + meandeg + degree(0:2) + concurrent,
            keep.tedgelist = TRUE)

#summary table
print(dx)
#summary plot
plot(dx)

#Extract the simulated raw network summary stats
# nwstats1 <- get_nwstats(dx, sim = 1)
# head(nwstats1, 20)
# 
# par(mfrow = c(1, 2))
# plot(dx, type = "duration")
# plot(dx, type = "dissolution")


#timed edge list - tracking how long each edge lasts and who it is between
# tel <- as.data.frame(dx, sim = 1)
# head(tel, 20)

#MCMC Chain diagnostics
# dx.static <- netdx(est, nsims = 10000, dynamic = FALSE)
# print(dx.static)
# par(mfrow = c(1,1))
# plot(dx.static, sim.lines = TRUE, sim.lwd = 0.1)
# nwstats2 <- get_nwstats(dx.static)
# head(nwstats2, 20)




##### Now to the Epidemic Model #####

#Parameter Definitions 
# act.rate -> # of acts that occur within a partnership per time unit
# inf.prob -> risk of transmission given an act with an infected person
# rec.rate -> speed at which infected people return to susceptibility
param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.1)

#Initialize Network
#i.num -> number of infected individuals (randomly selected in this case but you can specify who if you want)
init <- init.net(i.num = 10)

#Control Settings
#type -> SIS mean Susceptible, Infected, Susceptible
#nsims -> number of simulations
#nsteps -> number of time steps
#ncores -> for parallelization
control <- control.net(type = "SIS", nsims = 5, nsteps = 500, ncores = numCores)



#Simulating Epidemic
#est -> network from above
#param -> parameters defined above
#init -> intial conditions defined above
#control -> control settings defined above
sim <- netsim(est, param, init, control)
print(sim)


#Model Analysis
#Compartment Averages
par(mfrow = c(1, 1))
plot(sim, ylim = c(0,600), xlim = c(0,500))


# par(mfrow = c(1, 1))
# #plot each sim
# plot(sim, sim.lines = TRUE, mean.line = FALSE, qnts = FALSE, popfrac = TRUE)
#plot smoothed summary prevalence instead of incidence on y-axis
plot(sim, mean.smooth = FALSE, qnts = 1, qnts.smooth = FALSE, popfrac = TRUE, ylim = c(0,1.05), xlim = c(0,500))

par(mfrow = c(1,1))
#plotting flow into the compartment at each time step
plot(sim, y = c("si.flow", "is.flow"), qnts = FALSE, 
     ylim = c(0, 2), legend = TRUE, main = "Flow Sizes")

#Network Plots
par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
simNum <- 5
#Simulation 1 - Time Step 1
plot(sim, type = "network", col.status = TRUE, at = 1, sims = simNum)
#Simulation 2 - Time Step 500
plot(sim, type = "network", col.status = TRUE, at = 500, sims = simNum)

# #Time specific summaries
# summary(sim, at = 500)
# 
# 
#Extract Simulated Data
df <- as.data.frame(sim)
summary(df)
# 
# #Average the simulations together
# df <- as.data.frame(sim, out = "mean")
# head(df, 10)
# tail(df, 10)
# 
# 
# #Save Network Dynamics
# nw1 <- get_network(sim, sim = 1)
# nw1
# #Network dataframe
# nwdf <- as.data.frame(nw1)
# head(nwdf, 25)

# transmission events
tm1 <- get_transmat(sim, sim = 1)
tm1


# df <- as.data.frame(sim)
# df.mean <- as.data.frame(sim, out = "mean")
# 
# library(ggplot2)
# ggplot() +
#   geom_line(data = df, mapping = aes(time, i.num, group = sim), alpha = 0.25,
#             lwd = 0.25, color = "firebrick") +
#   geom_bands(data = df, mapping = aes(time, i.num),
#              lower = 0.1, upper = 0.9, fill = "firebrick") +
#   geom_line(data = df.mean, mapping = aes(time, i.num)) +
#   theme_minimal()

  # 7.23.2025
  # Amy Moore
  # SISMID - NME
  # 34 - Open Population Lab



  #Load Packages
library(EpiModel)


    ##### Part 1: Explore a model with increased birth rate compared to death rate


#Initialize Network
nw <- network_initialize(n = 500)  
nw <- set_vertex_attribute(nw, attrname = "risk", value = rep(0:1, each = 250))



##### Growing Population Size #####
# Assume Entry rate is higher than Exit rates
# Death Rate is the same among all groups

  #Formation Process
formation <- ~edges + nodefactor("risk") + nodematch("risk")
#Values chosen to match a mean degree of 0.75 in group 1 and 90% of group 1 edges are between two group 1 members.
target.stats <- c(125, 187.5, 112.5)

  #Dissolution Process
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 
                               duration = 40, d.rate = 0.001)
coef.diss

#estimate the network
est1 <- netest(nw, formation, target.stats, coef.diss)
summary(est1)

#diagnostics
dx1 <- netdx(est1, nsims = 10, nsteps = 1000, ncores = 4,
             nwstats.formula = ~edges + 
               nodefactor("risk", levels = NULL) + 
               nodematch("risk", diff = TRUE))
dx1
plot(dx1)


#Epidemic Simulation
#a.rate -> Arrival rate
#ds.rate -> departure rate from susceptible
#di.rate -> departure rate from infected
param <- param.net(inf.prob = 0.1, act.rate = 5,
                   a.rate = 0.0015, ds.rate = 0.001, di.rate = 0.001)

init <- init.net(i.num = 50)

control <- control.net(type = "SI", nsteps = 300, nsims = 20, ncores = 8, 
                       resimulate.network = TRUE, epi.by = "risk", tergmLite = TRUE,
                       nwstats.formula = ~edges + 
                         nodefactor("risk", levels = NULL) + 
                         nodematch("risk", diff = TRUE) + 
                         meandeg)

sim1 <- netsim(est1, param, init, control)

#Post-Sim Network Diagnostics
sim1
plot(sim1, type = "formation", plots.joined = FALSE)

#individual terms
par(mfrow = c(1, 2))
plot(sim1, type = "formation", stats = "edges", qnts = FALSE, sim.lines = TRUE)
plot(sim1, type = "formation", stats = "meandeg", qnts = FALSE, sim.lines = TRUE)


#Demographic outcomes

#Overall Population Size
par(mfrow = c(1, 2))
plot(sim1, y = "num", legend = TRUE, sim.lines = TRUE)
plot(sim1, y = c("num.risk0", "num.risk1"), legend = TRUE, sim.lines = TRUE)

#Compartment Flows
par(mfrow = c(1, 1))
plot(sim1, y = c("a.flow", "ds.flow", "di.flow"), mean.lwd = 2.5,
     qnts = FALSE, legend = TRUE, ylim = c(0, 1))

#Standardized Flows
sim1 <- mutate_epi(sim1, ds.flow.st = ds.flow / s.num,
                   di.flow.st = di.flow / i.num)
plot(sim1, y = c("ds.flow.st", "di.flow.st"), 
     qnts = FALSE, legend = TRUE, mean.lwd = 2.5)


#Epidemic Outcomes

#Prevalence
plot(sim1, y = "i.num", qnts = 1, main = "Total Prevalence", ylim = c(0, 500))

#Prevalence by group
plot(sim1, y = c("i.num.risk0", "i.num.risk1"),  legend = TRUE, qnts = 1,
     ylim = c(0, 500), main = "Prevalence by Group")


# Do you need to adjust the dissolution_coefs parameterization here?
    ### Depends, do we care that number of edges is increasing?
    ### Mean degree is staying the same, but population is growing
par(mfrow = c(1, 2))
plot(sim1, type = "formation", stats = "edges", qnts = FALSE, sim.lines = TRUE)
plot(sim1, type = "formation", stats = "meandeg", qnts = FALSE, sim.lines = TRUE)


# What happens to the overall population structure (population size)?
    ### Total Population is increasing throughout
    ### Arrival Flow is steady overtime
    ### Flow into Susceptible goes down because more people get infected
    ### Flow into Infected goes up because more people get infected
par(mfrow = c(1, 2))
plot(sim1, y = "num", legend = TRUE, sim.lines = TRUE)
plot(sim1, y = c("a.flow", "ds.flow", "di.flow"), mean.lwd = 2.5,
     qnts = FALSE, legend = TRUE, ylim = c(0, 1))


# What happens to network structure (edges and meandeg post-simulation diagnostics)?
    ### As population increases, number of edges increases
    ### As population increases, mean degree stays constant

par(mfrow = c(1, 2))
plot(sim1, type = "formation", stats = "edges", qnts = FALSE, sim.lines = TRUE)
plot(sim1, type = "formation", stats = "meandeg", qnts = FALSE, sim.lines = TRUE)




    ##### Part 2: assume a higher overall death rates with disease-related mortality #####
# Similar to above but changing back to original birth rate
# and changing death rate to vary by disease status

  #try 1 - simple average of ds.rate and di.rate
coef.diss2 <- dissolution_coefs(dissolution = ~offset(edges), 
                                duration = 40, d.rate = 0.0019)
coef.diss2

  #try 2 - increase effective d.rate closer to di.rate
coef.diss2 <- dissolution_coefs(dissolution = ~offset(edges), 
                                duration = 40, d.rate = 0.0023)
coef.diss2

param2 <- param.net(inf.prob = 0.1, act.rate = 5,
                   a.rate = 0.001, ds.rate = 0.0012, di.rate = 0.0025)



est2 <- netest(nw, formation, target.stats, coef.diss2)
sim2 <- netsim(est2, param2, init, control)

par(mfrow = c(1, 2))
plot(sim2, y = "num", legend = TRUE, sim.lines = TRUE)
plot(sim2, y = c("num.risk0", "num.risk1"), legend = TRUE, sim.lines = TRUE)

#Number of edges
par(mfrow = c(1, 1))
plot(sim2, type = "formation", stats = "edges", qnts = FALSE, sim.lines = TRUE)

###  # edges is slightly lower than expected because deaths take away edges and we didn't adjust for it

# mean degree
par(mfrow = c(1, 1))
plot(sim2, type = "formation", stats = "meandeg", qnts = FALSE, sim.lines = TRUE)
abline(h = 0.5, lty = 2, lwd = 2)

print(sim2)






##### Model 3: Varying Population Size #####
#higher rate of disease-specific mortality

#new parameters - adjusting death rate
param3 <- param.net(inf.prob = 0.1, act.rate = 5,
                    a.rate = 0.001, ds.rate = 0.001, di.rate = 0.002)

#adjust for differing death rate
### 0.0018 is a trial and error approximation of a weighted average 
### of person time spent in each disease state by the death rates for each state
coef.diss3 <- dissolution_coefs(dissolution = ~offset(edges), 
                                duration = 40, d.rate = 0.0018)
coef.diss3


#Fit Network and Epidemic
est3 <- netest(nw, formation, target.stats, coef.diss3)
sim3 <- netsim(est3, param3, init, control)



par(mfrow = c(1, 2))
#population size - declines over time
plot(sim3, y = "num", sim.lines = TRUE, 
     main = "Population Size", legend = FALSE)
#disease prevalence - increases over time
plot(sim3, y = "i.num", sim.lines = TRUE, 
     main = "Disease Prevalence", legend = FALSE)

#number of edges - decreases as population decreases
par(mfrow = c(1, 1))
plot(sim3, type = "formation", stats = "edges", qnts = FALSE, sim.lines = TRUE)


#mean degree - preserved over time
par(mfrow = c(1, 1))
plot(sim3, type = "formation", stats = "meandeg", qnts = FALSE, sim.lines = TRUE)
abline(h = 0.5, lty = 2, lwd = 2)
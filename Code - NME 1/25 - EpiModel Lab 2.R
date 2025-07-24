
  # 7.22.2025
  # Amy Moore
  # SISMID - NME
  # 25 - EpiModel Lab 2

#Load Packages
library(EpiModel)

#Set Seed
set.seed(12345)

#Initialization of the Network

#Set sample size for each group
num.g1 <- num.g2 <- 250
#Set total sample size
n <- num.g1 + num.g2
#initialize the network with n nodes
nw <- network_initialize(n)
nw

#define a list of the nodal attribute "group"
group <- rep(1:2, times = c(num.g1, num.g2))
#set the attribute on the network
nw <- set_vertex_attribute(nw, attrname = "group", value = group)
nw

#proof of group assignment
get_vertex_attribute(nw, "group")



#Network Model Specification

#overall mean degree (equal for both groups since all ties are cross-group)
meandeg <- 0.66

#set concurrency values for each group 
pconc.g1 <- c(0.05)
pconc.g2 <- c(0.15)



#Formation of the Network Model
#edges -> # of total edges
#nodematch("group") -> # of edges between members of the same group --> set to 0 to prevent within-group edges
#concurrent(...) -> # of concurrent nodes in group 1 and group 2 reported separately
formation <- ~edges + nodematch("group") + concurrent(by = "group", levels = NULL)

#set edges to match expected number based on mean degree
#set nodematch to zero (no edges between members of the same group)
#set concurrent grp1 to expected number of grp1 members with concurrent edges
#set concurrent grp2 to expected number of grp2 members with concurrent edges 
target.stats <- c(n*meandeg/2, 0, num.g1*pconc.g1, num.g2*pconc.g2)
target.stats
### Totally fine for the target stats to not be whole numbers


#Dissolution Model
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
coef.diss


#Network Model Estimation
est <- netest(nw, formation, target.stats, coef.diss)

#Diagnostics
dx <- netdx(est, nsims = 25, nsteps = 500, ncores = 5,
            nwstats.formula = ~edges + concurrent(by = "group") + degree(0:3, by = "group"))
dx
plot(dx)





##### Epidemic model #####

  # Parameterization
#Rec.rate --> implies that this is an SIR or SIS model
#inf.rate --> rate of infection
  ### set inf.prob = 0.6 to increase infection probability in group 1
param <- param.net(inf.prob = 0.6, inf.prob.g2 = 0.3,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)

param.og <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.3,
                   rec.rate = 0.02, rec.rate.g2 = 0.02)

param.v3 <- param.net(inf.prob = 0.6, inf.prob.g2 = 0.1,
                      rec.rate = 0.02, rec.rate.g2 = 0.02)

  #Initial Conditions
#r.num -> implies "recovered" category
init <- init.net(i.num = 10, i.num.g2 = 10, 
                 r.num = 0, r.num.g2 = 0)


  #Control Settings
control <- control.net(type = "SIR", nsims = 5, nsteps = 500, ncores = 5)



  # Simulations
sim <- netsim(est, param, init, control)
sim

sim.og <- netsim(est, param.og, init, control)
sim.og

sim.v3 <- netsim(est, param.v3, init, control)
sim.v3



  #Analysis

#Compartment Flow 
par(mar = c(3,3,2,1))
par(mfrow = c(1, 1))
plot(sim.og, main = "Disease State Sizes")
plot(sim, main = "Disease State Sizes")
plot(sim.v3, main = "Disease State Sizes")
###The group with less concurrency has higher incidence

#Incidence
par(mfrow = c(1,2))
plot(sim.og, y = c("i.num", "i.num.g2"), popfrac = TRUE,
     qnts = FALSE, ylim = c(0, 0.4), legend = TRUE)
plot(sim.og, y = c("si.flow", "si.flow.g2"), 
     qnts = FALSE, ylim = c(0, 2.5), legend = TRUE)
par(mfrow = c(1,2))
plot(sim, y = c("i.num", "i.num.g2"), popfrac = TRUE,
     qnts = FALSE, ylim = c(0, 0.4), legend = TRUE)
plot(sim, y = c("si.flow", "si.flow.g2"), 
     qnts = FALSE, ylim = c(0, 2.5), legend = TRUE)
par(mfrow = c(1,2))
plot(sim.v3, y = c("i.num", "i.num.g2"), popfrac = TRUE,
     qnts = FALSE, ylim = c(0, 0.4), legend = TRUE)
plot(sim.v3, y = c("si.flow", "si.flow.g2"), 
     qnts = FALSE, ylim = c(0, 2.5), legend = TRUE)


#Network Plots
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
plot(sim, type = "network", col.status = TRUE, at = 50, 
     sims = "mean", shp.g2 = "square")
plot(sim, type = "network", col.status = TRUE, at = 100, 
     sims = "mean", shp.g2 = "square")


#Computing Cumulative Incidence

#Extract Average Simulation Data
df.og <- as.data.frame(sim.og, out = "mean")
#cumulative incidence
c(sum(df.og$si.flow, na.rm = TRUE),sum(df.og$si.flow.g2, na.rm = TRUE))


#Extract Average Simulation Data
df <- as.data.frame(sim, out = "mean")
#cumulative incidence
c(sum(df$si.flow, na.rm = TRUE),sum(df$si.flow.g2, na.rm = TRUE))

#Extract Average Simulation Data
df.v3 <- as.data.frame(sim.v3, out = "mean")
#cumulative incidence
c(sum(df.v3$si.flow, na.rm = TRUE),sum(df.v3$si.flow.g2, na.rm = TRUE))



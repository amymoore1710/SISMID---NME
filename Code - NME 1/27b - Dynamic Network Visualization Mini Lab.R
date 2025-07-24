
  # 7.22.2025
  # Amy Moore
  # SISMID - NME
  # 27b - Dynamic Network Visualization Mini Lab

#Load Packages
library(EpiModel)
library("ndtv")
library("htmlwidgets")

#Set Seed
set.seed(12345)

#Initialization of the Network

#Set sample size for each group
num.g1 <- num.g2 <- 50
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
meandeg <- 2

#set concurrency values for each group 
pconc.g1 <- c(0.50)
pconc.g2 <- c(0.50)



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

# #Diagnostics
# dx <- netdx(est, nsims = 25, nsteps = 500, ncores = 5,
#             nwstats.formula = ~edges + concurrent(by = "group") + degree(0:3, by = "group"))
# dx
# plot(dx)





##### Epidemic model #####

# Parameterization
#Rec.rate --> implies that this is an SIR or SIS model
#inf.rate --> rate of infection
param <- param.net(inf.prob = 0.6, inf.prob.g2 = 0.3,
                      rec.rate = 0.02, rec.rate.g2 = 0.02)

#Initial Conditions
#r.num -> implies "recovered" category
init <- init.net(i.num = 10, i.num.g2 = 10, 
                 r.num = 0, r.num.g2 = 0)


#Control Settings
control <- control.net(type = "SIR", nsims = 1, nsteps = 25, ncores = 5)



# Simulations
sim <- netsim(est, param, init, control)
sim

#Fit Network and Visualize
nw <- get_network(sim)
nw <- color_tea(nw)
slice.par <- list(start = 1, end = 25, interval = 1,
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))
compute.animation(nw, slice.par = slice.par, verbose = TRUE)
community <- get_vertex_attribute(nw, "group")
community.shape <- ifelse(community == 1, 4, 50) #defines shape based on community

render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.sides = community.shape,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))







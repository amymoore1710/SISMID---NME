
  # 7.23.2025
  # Amy Moore
  # SISMID - NME
  # 33 - Feedback through Nodal Attribute Changes

  #initialize network
n <- 500
nw <- network_initialize(n)

#randomly assign 100 nodes to be infected initially
prev <- 0.2
infIds <- sample(1:n, n*prev)
nw <- set_vertex_attribute(nw, "status", "s")
nw <- set_vertex_attribute(nw, "status", "i", infIds)
get_vertex_attribute(nw, "status")

#expected number of edges between infected
mean.deg.inf <- 0.3
inedges.inf <- mean.deg.inf * n * prev
inedges.inf

#expected number of edges between susceptible
mean.deg.sus <- 0.8
inedges.sus <- mean.deg.sus * n * (1 - prev)
inedges.sus

#expected total number of edges
edges <- (inedges.inf + inedges.sus)/2
edges



p <- inedges.sus/(edges*2)
q <- 1 - p
nn <- p^2
np <- 2*p*q
pp <- q^2
round(nn + pp, 3)


fit <- netest(nw,
              formation = ~edges + nodefactor("status"),
              target.stats = c(edges, inedges.sus),
              coef.diss = dissolution_coefs(~offset(edges), duration = 1))
sim <- netdx(fit, dynamic = FALSE, nsims = 1e4,
             nwstats.formula = ~edges + nodematch("status"))

stats <- get_nwstats(sim)
head(stats)

round(mean(stats$nodematch.status/stats$edges), 3)

nmatch <- edges * 0.91


  #Model 1
formation <- ~edges + nodefactor("status") + nodematch("status")
target.stats <- c(edges, inedges.sus, nmatch)

coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 50)

est <- netest(nw, formation, target.stats, coef.diss)

dx <- netdx(est, nsims = 10, nsteps = 500, ncores = 5,
            nwstats.formula = ~edges + 
              meandeg + 
              nodefactor("status", levels = NULL) + 
              nodematch("status"), verbose = FALSE)
dx

plot(dx)

plot(dx, type = "dissolution")

  #Model 2
est2 <- netest(nw, formation = ~edges, target.stats = edges, coef.diss)

dx2 <- netdx(est2, nsims = 10, nsteps = 1000, ncores = 5,
             nwstats.formula = ~edges + 
               meandeg + 
               nodefactor("status", levels = NULL) + 
               nodematch("status"), verbose = FALSE)
dx2

plot(dx2)



    ##### Epidemic Model #####

param <- param.net(inf.prob = 0.03)
init <- init.net()
control <- control.net(type = "SI", nsteps = 500, nsims = 5, ncores = 5,
                       resimulate.network = TRUE, tergmLite = TRUE,
                       nwstats.formula = ~edges + 
                         meandeg + 
                         nodefactor("status", levels = NULL) + 
                         nodematch("status"))
  #SIM 1 - heterogeniety in mixing
sim <- netsim(est, param, init, control)
  #Classical Mixing
sim2 <- netsim(est2, param, init, control)


  #Results
par(mfrow = c(1,2))
plot(sim, main = "Seroadaptive Behavior")
plot(sim2, main = "No Seroadaptive Behavior")

#Compartment - Infected overtime
par(mfrow = c(1,1))
plot(sim, y = "i.num", popfrac = TRUE, sim.lines = FALSE, qnts = 1)
plot(sim2, y = "i.num", popfrac = TRUE, sim.lines = FALSE, qnts = 1, 
     mean.col = 2, qnts.col = 2, add = TRUE)
legend("topleft", c("Serosort", "Non-Serosort"), lty = 1, lwd = 3,
       col = c(4, 2), cex = 0.9, bty = "n")


#Network Stats
plot(sim, type = "formation", qnts = FALSE, sim.lines = TRUE)

#Comparing Node Factor Status between S and I
plot(sim, type = "formation", qnts = FALSE, sim.lines = TRUE,
     stats = c("nodefactor.status.s", 
               "nodefactor.status.i"))

plot(sim, type = "formation", 
     stats = c("edges", "nodematch.status"), 
     qnts = FALSE, sims = 1, sim.lwd = 2)

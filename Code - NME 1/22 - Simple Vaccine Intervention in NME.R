
  # 7.22.2025
  # Amy Moore
  # SISMID - NME
  # 22 - Simple Vaccine Intervention in NME

  #Start with simple tergm
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 50
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)

  #Simple SI model example
#inter.eff -> Effectiveness of intervention
#inter.start -> When intervention begins
param <- param.net(inf.prob = 0.5, inter.eff = 0.96, inter.start = 25)
init <- init.net(i.num = 5)
control <- control.net(type = "SI", nsteps = 100, nsims = 10, ncores = 5)
sim <- netsim(est, param, init, control)
plot(sim)

  #Simple SIS model example with less effective intervention
param <- param.net(inf.prob = 0.5, inter.eff = 0.8, inter.start = 100, 
                   rec.rate = 0.07)
init <- init.net(i.num = 10)
control <- control.net(type = "SIS", nsteps = 250, nsims = 10, ncores = 5)
sim <- netsim(est, param, init, control)
plot(sim)

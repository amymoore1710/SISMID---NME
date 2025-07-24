
  # 7.22.2025
  # Amy Moore
  # SISMID - NME
  # 27 - Dynamic Network Visualization

  #Load Packages
library("EpiModel")
library("ndtv")
library("htmlwidgets")

  #Set Seed
set.seed(1234)

  #Build a simple Network Model
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 40
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss)

  #Model Epidemic on the Network
param <- param.net(inf.prob = 1)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 25, nsims = 1)
sim <- netsim(est, param, init, control)



  #Extract Network
nw <- get_network(sim)
nw

  #Color Diseased Nodes
nw <- color_tea(nw, verbose = FALSE)


  #Revisit static visuals
tm <- get_transmat(sim)
tm

par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))
plot(tm, style = "phylo")
plot(tm, style = "network", displaylabels = TRUE)
plot(tm, style = "transmissionTimeline")

  #Proximity Timeline
par(mfrow = c(1, 1), mar = c(3, 3, 2, 1))
proximity.timeline(nw, 
                   vertex.col = "ndtvcol",
                   spline.style = "color.attribute",
                   mode = "sammon",
                   default.dist = 10,
                   chain.direction = "reverse")


  #Dynamic Network Movie
#control settings
slice.par <- list(start = 1, end = 25, interval = 1, 
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))

  #determine where to place nodes in space for each time step
compute.animation(nw, slice.par = slice.par, verbose = TRUE)

  #render animation
render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))



  ##### What if we want to model without Concurrency? #####

  #How much concurrency is occuring?
ppois(1, lambda = 0.8, lower.tail = FALSE) * 100

  #epidemic model
nw <- network_initialize(n = 100)
formation <- ~edges + concurrent
target.stats <- c(40, 0)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss)
set.seed(12345)
param <- param.net(inf.prob = 1)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 25, nsims = 1)
sim <- netsim(est, param, init, control)

  #Static Transition Visuals
par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))
tm <- get_transmat(sim)
plot(tm, style = "phylo")
plot(tm, style = "network", displaylabels = TRUE)
plot(tm, style = "transmissionTimeline")


nw <- get_network(sim)
nw <- color_tea(nw, verbose = FALSE)
compute.animation(nw, slice.par = slice.par)
render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))




  ##### Very Long Durations #####

  #Network Building
nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 40
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1e5)
est <- netest(nw, formation, target.stats, coef.diss)


  #Epidemic Model
param <- param.net(inf.prob = 1)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 25, nsims = 1)
sim <- netsim(est, param, init, control)

  #Visual Movie
nw <- get_network(sim)
nw <- color_tea(nw, verbose = FALSE)
compute.animation(nw, slice.par = slice.par)
render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))



  ##### Very Short Duration #####
#Essentially what a typical SIR model assumes

nw <- network_initialize(n = 100)
formation <- ~edges
target.stats <- 40
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 2)
est <- netest(nw, formation, target.stats, coef.diss)

param <- param.net(inf.prob = 1)
init <- init.net(i.num = 1)
control <- control.net(type = "SI", nsteps = 25, nsims = 1)
sim <- netsim(est, param, init, control)

nw <- get_network(sim)
nw <- color_tea(nw, verbose = FALSE)
compute.animation(nw, slice.par = slice.par)
render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))



    ##### Triangles #####

  #Build Network
nw <- network_initialize(n = 100)
formation <- ~edges + gwesp(0, TRUE)
target.stats <- c(40, 25)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
est <- netest(nw, formation, target.stats, coef.diss)
sim <- netsim(est, param, init, control)

  #Create Epidemic
nw <- get_network(sim)
nw <- color_tea(nw, verbose = FALSE)
compute.animation(nw, slice.par = slice.par)
render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))



    ##### Extra Nodal Attributes #####

  #Create Network with attributes
nw <- network_initialize(n = 100)
nw <- set_vertex_attribute(nw, "community", rbinom(100, 1, 0.5))
nw <- set_vertex_attribute(nw, "age", sample(18:50, 100, TRUE))
nw

  #Forming the network
formation <- ~edges + nodematch("community") + nodefactor("community") +
  absdiff("age") + concurrent
target.stats <- c(50, 40, 70, 100, 30)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges) + offset(nodematch("community")),
                               duration = c(20, 10))
est <- netest(nw, formation, target.stats, coef.diss)

  #Add Epidemic
param <- param.net(inf.prob = 1)
init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsteps = 25, nsims = 1)
sim <- netsim(est, param, init, control)

  #Fit Network and Visualize
nw <- get_network(sim)
nw <- color_tea(nw)
slice.par <- list(start = 1, end = 25, interval = 1,
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))
compute.animation(nw, slice.par = slice.par, verbose = TRUE)
community <- get_vertex_attribute(nw, "community")
community.shape <- ifelse(community == 1, 4, 50) #defines shape based on community

age <- get_vertex_attribute(nw, "age")
age.size <- age/25 #defines size of node based on age
render.d3movie(
  nw,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = age.size,
  vertex.sides = community.shape,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))






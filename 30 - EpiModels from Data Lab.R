
  # 7.22.2025
  # Amy Moore
  # SISMID - NME
  # 30 - EpiModels from Data Lab

  #Load Packages
library(EpiModel)
library(ergm.ego)
library(ndtv)
#install.packages("Rglpk")

data("faux.mesa.high")
mesa <- faux.mesa.high

# model 5 from yesterday with degree 0 term

formation <- ~ edges + 
  nodefactor("Grade") + nodematch("Grade", diff=T) +
  nodefactor("Race") + nodematch("Race", diff=T) +
  nodefactor("Sex") +   nodematch("Sex") +
  degree(0)


targets <- summary(mesa ~ edges + 
                     nodefactor("Grade") + nodematch("Grade", diff=T) +
                     nodefactor("Race") + nodematch("Race", diff=T) +
                     nodefactor("Sex") +   nodematch("Sex") +
                     degree(0))

targets


set.seed(1)
myfit <- netest(mesa,
                formation = formation,
                target.stats = targets,
                coef.diss = dissolution_coefs(~offset(edges), 60))


mydx <- netdx(myfit, nsims = 25, nsteps = 500, ncores = 5)
plot(mydx)



  #Epidemic Models
myparam <- param.net(inf.prob = 0.2,
                     act.rate = 1.8,
                     rec.rate = 0.1)
myinit <- init.net(i.num = 10)

mycontrol <- control.net("SIS", nsteps = 100, nsims = 20,   
                         epi.by = "Grade",
                         nwstats.formula = ~edges + nodematch("Grade") + 
                           degree(0:5) + triangles, verbose = TRUE)

mySIS <- netsim(myfit, param = myparam, init = myinit, control = mycontrol)

mySIS <- mutate_epi(mySIS, 
                    prev.7 = i.num.Grade7 / num.Grade7,
                    prev.8 = i.num.Grade8 / num.Grade8,
                    prev.9 = i.num.Grade9 / num.Grade9,
                    prev.10 = i.num.Grade10 / num.Grade10,
                    prev.11 = i.num.Grade11 / num.Grade11,
                    prev.12 = i.num.Grade12 / num.Grade12)

plot(mySIS, y = c("prev.7", "prev.8", "prev.9", "prev.10",
                  "prev.11", "prev.12"), 
     qnts.alpha = 0.2, qnts = 0.5, 
     main="Prevalence Rate",
     legend = TRUE,
     xlim = c(0,mycontrol$nsteps+25))



  #Dynamic Visualization
myND <- get_network(mySIS)
myND <- color_tea(myND, verbose = FALSE)

slice.par <- list(start = 1, end = 20, interval = 1, 
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))

compute.animation(myND, slice.par = slice.par, verbose = TRUE)
render.d3movie(
  myND,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  output.mode = "htmlWidget")


  # 7.22.2025
  # Amy Moore
  # SISMID - NME
  # 29 - Network Data Inputs

  #Load Packages
library(EpiModel)
library(Rglpk)
library(ndtv)

  #Load Data
data("faux.mesa.high")
mesa <- faux.mesa.high
class(mesa)

  #Looking at the network
mesa

  #Looking at attribute effects
fill.col <- RColorBrewer::brewer.pal(8, "Set2") # 8 colors

# Overall
barplot(degreedist(mesa)/network.size(mesa),  xaxt = "n",
        main = "Overall degree distribution", 
        xlab = "Degree", ylab = "Fraction of nodes",
        col = "lightblue")

degs <- gsub("degree", "", names(degreedist(mesa)))
pos <- barplot(degreedist(mesa), plot=F)

axis(1, at = pos, labels = degs)
points(x=pos, y=dbinom(c(0:10, 13), size=205, prob=1.98/204), col="red", pch=19)
text(7, 0.2, paste("Mean degree: ", round(summary(mesa ~ meandeg), 2)))
legend(10, 0.27, pch=19, col="red", "RBG mean", cex=.5)


par(mfrow=c(1,2))
# by Grade

barplot(summary(mesa ~ nodefactor("Grade", levels=T)) / table(mesa %v% "Grade"), 
        horiz = T, 
        yaxt = "n", col = fill.col,
        main = "Mean degree by Grade", xlab = "Mean degree", ylab = "Grade")
grades <- c(7:12)
pos <- barplot(summary(mesa ~ meandeg:nodefactor("Grade", levels=TRUE)), 
               horiz = T, plot=F)
axis(2, at = pos, labels = grades)

# by sex

barplot(summary(mesa ~ nodefactor("Sex", levels=T)) / table(mesa %v% "Sex"), horiz = T, 
        yaxt = "n", col=fill.col[7:8],
        main = "Mean degree by Sex", xlab = "Mean degree", ylab = "Sex")
sex <- c("F", "M")
pos <- barplot(summary(mesa ~ meandeg:nodefactor("Sex", levels=TRUE)), 
               horiz = T, plot=F)
axis(2, at = pos, labels=sex)


# Isolates by sex
round(100*(summary(mesa ~ degree(0, by="Sex")) / table(mesa %v% "Sex")), 1)


set.seed(0)
coords <- plot(mesa)

mixingmatrix(mesa, "Grade")
mixingmatrix(mesa, "Sex")

par(mfrow = c(1,2), oma=c(0,0,3,0))

plot(mesa, coord=coords,
     vertex.col='Grade', main="By Grade", cex.main = 1)
legend('bottomleft',fill=7:12,
       legend=paste('Grade',7:12),cex=0.5)

plot(mesa, coord=coords,
     vertex.col='Sex', main="By Sex", cex.main = 1)
legend('bottomleft',fill=1:2,
       legend=paste('Sex',c("F","M")),cex=0.5)

mtext("Mesa friendship network", side=3, outer=T, cex=1.5)






  ##### Specify the Model #####

  #Model Formation
#Node factor is like a first-order term effect and Node match is a higher order term
formation <- ~ edges + 
  nodefactor("Grade") + nodematch("Grade") +
  nodefactor("Sex") +   nodematch("Sex") +
  degree(0, by = "Sex")

  #Take targets from the sample
#rather than trying to come up with numbers, just use data observed
targets <- summary(mesa ~ edges + 
                     nodefactor("Grade") + nodematch("Grade") +
                     nodefactor("Sex") +   nodematch("Sex") +
                     degree(0, by = "Sex"))

targets

  #Fit the tergm
set.seed(0)
myfit <- netest(mesa,
                formation = formation,
                target.stats = targets,
                coef.diss = dissolution_coefs(~offset(edges), 60))

  #Diagnostics
mydx <- netdx(myfit, nsims = 25, nsteps = 100, ncores = 5)

plot(mydx)

  #Run diagnostics on stats not in the model
summary(mesa ~ meandeg + gwdegree(.5, fixed=T) + kstar(2) + triangles)
set.seed(1)
mydx <- netdx(myfit, nsims = 25, nsteps = 100, ncores = 5,
              nwstats.formula = ~ meandeg + gwdegree(.5, fixed=T) + 
                kstar(2) + triangles,
              keep.tedgelist = TRUE)
plot(mydx)


  ##### Epidemic Simulation #####
myparam <- param.net(inf.prob = 0.2,
                     act.rate = 1.8,
                     rec.rate = 0.1)
myinit <- init.net(i.num = 10)

# Here we set which attribute to use for epi stats breakdowns 
# (only one can be specified in base EpiModel) 
# and also which stats we want to monitor

mycontrol <- control.net("SIS", nsteps = 100, nsims = 10,   
                         epi.by = "Grade",
                         nwstats.formula = ~edges + nodematch("Grade") + 
                           degree(0:2, by = "Sex"),
                         verbose = TRUE)

  #Simulate Epidemic
set.seed(0)
mySIS <- netsim(myfit, param = myparam, init = myinit, control = mycontrol)

par(mar = c(3,3,2,1))
par(mfrow = c(1,1))
plot(mySIS, formation=T)

  #Prevalence by Grade
mySIS <- mutate_epi(mySIS, 
                    prev.7 = i.num.Grade7 / num.Grade7,
                    prev.8 = i.num.Grade8 / num.Grade8,
                    prev.9 = i.num.Grade9 / num.Grade9,
                    prev.10 = i.num.Grade10 / num.Grade10,
                    prev.11 = i.num.Grade11 / num.Grade11,
                    prev.12 = i.num.Grade12 / num.Grade12)

plot(mySIS, y = c("prev.7", "prev.8", "prev.9", "prev.10",
                  "prev.11", "prev.12"), 
     qnts.alpha = 0.2, qnts = 0.5, legend = TRUE)


  #Term Monitoring throughout Epidemic simulation
plot(mySIS, type = "formation", qnts = FALSE, sim.lines = TRUE)



  ##### Dynamic Visualization #####

myEpiSim <- get_network(mySIS)
myEpiSim <- color_tea(myEpiSim, verbose = FALSE)

slice.par <- list(start = 1, end = 20, interval = 1, 
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))

compute.animation(myEpiSim, slice.par = slice.par, verbose = TRUE)
render.d3movie(
  myEpiSim,
  render.par = render.par,
  plot.par = plot.par,
  vertex.cex = 0.9,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  output.mode = "htmlWidget")















    ##### Ego Centric Data #####

  #Load Package
library(ergm.ego)

  #Load Data
mesa.ego <- ergm.ego:::as.egor.network(mesa)
class(mesa.ego)

mesa.ego
summary(mesa.ego)
table(mesa %v% "Grade", mesa %v% "Sex")
table(mesa.ego$ego$Grade, mesa.ego$ego$Sex)

degreedist(mesa.ego, by="Sex", brgmod=T,
           main="Mesa degree distribution by Sex")
degreedist(mesa.ego, by="Grade", brgmod=T,
           main="Mesa degree distribution by Grade")


cbind(
  network = summary(mesa ~ edges + 
                      nodefactor("Grade") + nodematch("Grade") +
                      nodefactor("Sex") +   nodematch("Sex") +
                      degree(0, by = "Sex")),
  egocensus = summary(mesa.ego ~ edges + 
                        nodefactor("Grade") + nodematch("Grade") +
                        nodefactor("Sex") + nodematch("Sex") +
                        degree(0, by = "Sex")))

  #Fit the Network
set.seed(0)
myfit.ego <- netest(mesa.ego,
                    formation = formation,
                    coef.diss = dissolution_coefs(~offset(edges), 60))

round(cbind(
  network = myfit$coef.form.crude,
  egocensus = myfit.ego$coef.form.crude,
  diff = myfit$coef.form.crude - myfit.ego$coef.form.crude), 3)

  #Diagnostics
set.seed(0)
mydx.ego <- netdx(myfit, nsims = 25, nsteps = 500, ncores = 5)
plot(mydx.ego)






  ##### Scale Network to a different size #####

  #fit to n=1000
set.seed(0)
myfit.ego.pp1000 <- netest(mesa.ego,
                           formation = formation,
                           coef.diss = dissolution_coefs(~offset(edges), 60),
                           set.control.ergm.ego = ergm.ego::control.ergm.ego(
                             ppopsize = 1000))

  #fit to n=2000
set.seed(0)
myfit.ego.pp2000 <- netest(mesa.ego,
                           formation = formation,
                           coef.diss = dissolution_coefs(~offset(edges), 60),
                           set.control.ergm.ego = ergm.ego::control.ergm.ego(
                             ppopsize = 2000))

# network sizes:

cat("Network sizes: \n\n")
data.frame(
  network = network.size(myfit$newnetwork),
  egocensus = network.size(myfit.ego$newnetwork),
  egopp1000 = network.size(myfit.ego.pp1000$newnetwork),
  egopp2000 = network.size(myfit.ego.pp2000$newnetwork))

# scaling factor for network size:

cat("Network size scaling factor: \n\n")
data.frame(
  network = network.size(myfit$newnetwork),
  egocensus = network.size(myfit.ego$newnetwork),
  egopp1000 = network.size(myfit.ego.pp1000$newnetwork),
  egopp2000 = network.size(myfit.ego.pp2000$newnetwork)) / network.size(mesa)

# target stats:

cat("Target stats: \n\n")
data.frame(
  network = myfit$target.stats,
  egocensus = myfit.ego$target.stats,
  egopp1000 = myfit.ego.pp1000$target.stats,
  egopp2000 = myfit.ego.pp2000$target.stats)


# scaling factor for target stats:

cat("Target stat scaling: \n\n")
data.frame(
  network = myfit$target.stats,
  egocensus = myfit.ego$target.stats,
  egopp1000 = myfit.ego.pp1000$target.stats,
  egopp2000 = myfit.ego.pp2000$target.stats) / myfit$target.stats

# mean degree:

cat("Mean degree: \n\n")
data.frame(
  network = summary(myfit$newnetwork ~ meandeg),
  egocensus = summary(myfit.ego$newnetwork ~ meandeg),
  egopp1000 = summary(myfit.ego.pp1000$newnetwork ~ meandeg),
  egopp2000 = summary(myfit.ego.pp2000$newnetwork ~ meandeg))



round(data.frame(
  network = myfit$coef.form.crude,
  egocensus = myfit.ego$coef.form.crude,
  egopp1000 = myfit.ego.pp1000$coef.form.crude,
  egopp2000 = myfit.ego.pp2000$coef.form.crude,
  pct.diff2000 = (100*(myfit$coef.form.crude - myfit.ego.pp2000$coef.form.crude)/myfit$coef.form.crude)), 3)
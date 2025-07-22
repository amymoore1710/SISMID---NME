library(knitr)
opts_chunk$set(comment = NA)
par(mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))

library(EpiModel)

nw <- network_initialize(n = 500)

formation <- ~edges + concurrent + degrange(from = 4)

target.stats <- c(175, 110, 0)

coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 50)
coef.diss

args(netest)

est <- netest(nw, formation, target.stats, coef.diss)

# dx <- netdx(est, nsims = 10, nsteps = 1000,
#             nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent,
#             keep.tedgelist = TRUE)

parallel::detectCores()

dx <- netdx(est, nsims = 10, nsteps = 1000, ncores = 4,
            nwstats.formula = ~edges + meandeg + degree(0:4) + concurrent,
            keep.tedgelist = TRUE)

print(dx)

plot(dx)

nwstats1 <- get_nwstats(dx, sim = 1)
head(nwstats1, 20)

par(mfrow = c(1, 2))
plot(dx, type = "duration")
plot(dx, type = "dissolution")

tel <- as.data.frame(dx, sim = 1)
head(tel, 20)

dx.static <- netdx(est, nsims = 10000, dynamic = FALSE)
print(dx.static)

par(mfrow = c(1,1))
plot(dx.static, sim.lines = TRUE, sim.lwd = 0.1)

nwstats2 <- get_nwstats(dx.static)
head(nwstats2, 20)

param <- param.net(inf.prob = 0.4, act.rate = 2, rec.rate = 0.1)

init <- init.net(i.num = 10)

# control <- control.net(type = "SIS", nsims = 5, nsteps = 500)

control <- control.net(type = "SIS", nsims = 5, nsteps = 500, ncores = 5)

sim <- netsim(est, param, init, control)

print(sim)

par(mfrow = c(1, 1))
plot(sim)

par(mfrow = c(1, 2))
plot(sim, sim.lines = TRUE, mean.line = FALSE, qnts = FALSE, popfrac = TRUE)
plot(sim, mean.smooth = FALSE, qnts = 1, qnts.smooth = FALSE, popfrac = TRUE)

par(mfrow = c(1,1))
plot(sim, y = c("si.flow", "is.flow"), qnts = FALSE, 
     ylim = c(0, 25), legend = TRUE, main = "Flow Sizes")

par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
plot(sim, type = "network", col.status = TRUE, at = 1, sims = 1)
plot(sim, type = "network", col.status = TRUE, at = 500, sims = 1)

summary(sim, at = 500)

df <- as.data.frame(sim)
head(df, 10)
tail(df, 10)

df <- as.data.frame(sim, out = "mean")
head(df, 10)
tail(df, 10)

nw1 <- get_network(sim, sim = 1)
nw1

nwdf <- as.data.frame(nw1)
head(nwdf, 25)

tm1 <- get_transmat(sim, sim = 1)
head(tm1, 10)

df <- as.data.frame(sim)
df.mean <- as.data.frame(sim, out = "mean")

library(ggplot2)
ggplot() +
  geom_line(data = df, mapping = aes(time, i.num, group = sim), alpha = 0.25,
            lwd = 0.25, color = "firebrick") +
  geom_bands(data = df, mapping = aes(time, i.num),
             lower = 0.1, upper = 0.9, fill = "firebrick") +
  geom_line(data = df.mean, mapping = aes(time, i.num)) +
  theme_minimal()

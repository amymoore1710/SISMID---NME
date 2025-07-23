
  # 7.23.2025
  # Amy Moore
  # SISMID - NME
  # 32 - Epidemic Models with Demography


  #Load Packages
library(EpiModel)



  #Initialize Network
nw <- network_initialize(n = 500)  
nw <- set_vertex_attribute(nw, attrname = "risk", value = rep(0:1, each = 250))



    ##### Model 1: Stable Population Size #####

  #Formation Process
formation <- ~edges + nodefactor("risk") + nodematch("risk")
target.stats <- c(125, 187.5, 112.5)

  #Dissolution Process
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), 
                               duration = 40, d.rate = 0.001)
coef.diss








































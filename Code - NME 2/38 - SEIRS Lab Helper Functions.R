
  # 7.23.2025
  # Amy Moore
  # SISMID - NME II
  # 38 - SEIRS Lab Helper Functions


infect <- function(dat, at) {
  
  ## Uncomment this to run environment interactively
  # browser()
  
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  infTime <- get_attr(dat, "infTime")
  
  ## Parameters ##
  inf.prob <- get_param(dat, "inf.prob")
  act.rate <- get_param(dat, "act.rate")
  
  ## Find infected nodes ##
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  
  ## Initialize default incidence at 0 ##
  nInf <- 0
  
  ## If any infected nodes, proceed with transmission ##
  if (nElig > 0 && nElig < nActive) {
    
    ## Look up discordant edgelist ##
    del <- discord_edgelist(dat, at)
    
    ## If any discordant pairs, proceed ##
    if (!(is.null(del))) {
      
      # Set parameters on discordant edgelist data frame
      del$transProb <- inf.prob
      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      
      # Stochastic transmission process
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      
      # Keep rows where transmission occurred
      del <- del[which(transmit == 1), ]
      
      # Look up new ids if any transmissions occurred
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)
      
      # Set new attributes and transmission matrix
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        infTime[idsNewInf] <- at
        dat <- set_attr(dat, "status", status)
        dat <- set_attr(dat, "infTime", infTime)
        
        dat <- set_transmat(dat, del, at)
      }
    }
  }
  
  ## Save summary statistic for S->E flow
  dat <- set_epi(dat, "se.flow", at, nInf)
  
  return(dat)
}


progress <- function(dat, at) {
  
  ## Uncomment this to function environment interactively
  # browser()
  
  ## Attributes ##
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  
  ## Parameters ##
  ei.rate <- get_param(dat, "ei.rate")
  ir.rate <- get_param(dat, "ir.rate")
  rs.rate <- get_param(dat, "rs.rate")
  
  ## E to I progression process ##
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "e")
  nEligInf <- length(idsEligInf)
  
  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }
  
  ## I to R progression process ##
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }
  
  ## I to R progression process ##
  nSus <- 0
  idsEligSus <- which(active == 1 & status == "r")
  nEligSus <- length(idsEligSus)
  
  if (nEligSus > 0) {
    vecSus <- which(rbinom(nEligSus, 1, rs.rate) == 1)
    if (length(vecSus) > 0) {
      idsSus <- idsEligSus[vecSus]
      nSus <- length(idsSus)
      status[idsSus] <- "s"
    }
  }
  
  ## Write out updated status attribute ##
  dat <- set_attr(dat, "status", status)
  
  ## Save summary statistics ##
  dat <- set_epi(dat, "ei.flow", at, nInf)
  dat <- set_epi(dat, "ir.flow", at, nRec)
  dat <- set_epi(dat, "rs.flow", at, nSus)
  
  dat <- set_epi(dat, "e.num", at,
                 sum(active == 1 & status == "e"))
  dat <- set_epi(dat, "r.num", at,
                 sum(active == 1 & status == "r"))
  
  return(dat)
}
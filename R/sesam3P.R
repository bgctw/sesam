#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam3P <- function(
  ### Soil Enzyme Steady Allocation model including phosporous dynamics
  t,x,parms
){
  ##details<<
  ## Simplified version of sesam3s
  ## model corresponding to Seam3 with enzyme levels computed by quasi steady state
  ## Alpha as an explicit state variable that changes with turnover
  ## Simplified computation of target based on current revenue based on current alpha
  x <- pmax(unlist(x),1e-16)      # no negative masses
  # compute steady state enzyme levels for N and for C limitation
  dRPot <- parms$kR * x["R"]
  dLPot <- parms$kL * x["L"]
  immoPot <- parms$iB * x["I"]
  immoPPot <- parms$iBP * x["IP"]
  cnR <- x["R"]/x["RN"]
  cnL <- x["L"]/x["LN"]
  cnE <- parms$cnE
  cnB <- parms$cnB
  cpR <- x["R"]/x["RP"]
  cpL <- x["L"]/x["LP"]
  cpE <- parms$cpE
  cpB <- parms$cpB
  alpha <- x["alpha"]
  B <- x["B"]
  aeB <- parms$aE*B        # aeB without associanted growth respiration
  kmN <- parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  tvrB <- parms$tau*B      # microbial turnover
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * (1 - alpha)*aeB/(kmN + (1 - alpha)*aeB)
  decR <- dRPot * alpha*aeB/(kmN + alpha*aeB)
  #
  tvrERecycling <- parms$kNB*synE
  uNOrg <- parms$nu*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  uPOrg <- parms$nuP*(decL/cpL + decR/cpR + tvrERecycling/cpE)
  uC <- decL + decR + tvrERecycling
  CsynBC <- uC - rM - synE/parms$eps
  CsynBN <- cnB/parms$eps*(uNOrg + immoPot - synE/cnE)
  CsynBP <- cpB/parms$eps*(uPOrg + immoPPot - synE/cpE)
  CsynB <- min(CsynBC, CsynBN, CsynBP)
  if (CsynB > 0) {
    synB <- parms$eps*CsynB
    rG <- (1 - parms$eps)*CsynB
  } else {
    synB <- CsynB # with negative biomass change, do growth respiration
    rG <- 0
  }
  PhiB <- uNOrg - synB/cnB - synE/cnE
  PhiPB <- uPOrg - synB/cpB - synE/cpE
  alphaC <- computeSesam3sAllocationPartitioning(
    dR = dRPot, dL = dLPot, B = B
    , kmkN = kmN, aE = parms$aE
    , alpha = alpha
  )
  alphaN <- computeSesam3sAllocationPartitioning(
    dR = dRPot/cnR, dL = dLPot/cnL, B = B
    , kmkN = kmN, aE = parms$aE
    , alpha = alpha
  )
  alphaP <- computeSesam3sAllocationPartitioning(
    dR = dRPot/cpR, dL = dLPot/cpL, B = B
    , kmkN = kmN, aE = parms$aE
    , alpha = alpha
  )
  resBalance <- balanceAlphaBetweenElementLimitations(
    structure(c(alphaC, alphaN, alphaP), names = c("C","N","P"))
    , c(CsynBC, CsynBN, CsynBP)
    , tauB = parms$tau*B
  )
  alphaTarget <- resBalance$alpha
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) *  (parms$tau + abs(synB)/B)
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC - (synE/parms$eps + synB + rG + rM)
  respTvr <- (1 - parms$epsTvr) * tvrB
  # assuming same cnRatio of predators to equal cn ratio of microbes
  PhiTvr <- respTvr/parms$cnB
  PhiPTvr <- respTvr/parms$cpB
  #
  # tvr that feeds R pool, assume that N in SOM for resp (by epsTvr) is mineralized
  tvrC <-  +parms$epsTvr*tvrB   + (1 - parms$kNB)*synE
  tvrN <-  +parms$epsTvr*tvrB/parms$cnB   + (1 - parms$kNB)*synE/parms$cnE
  tvrP <-  +parms$epsTvr*tvrB/parms$cpB   + (1 - parms$kNB)*synE/parms$cpE
  # fluxes leaving the system (will be set in scen where trv does not feed back)
  tvrExC <- tvrExN <- tvrExP <- 0
  #
  leach <- parms$l*x["I"]
  PhiU <- (1 - parms$nu)*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  leachP <- parms$lP*x["IP"]
  PhiPU <- (1 - parms$nuP)*(decL/cpL + decR/cpR + tvrERecycling/cpE)
  #
  dB <- synB - tvrB
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dLP <- -decL/cpL   + parms$iL/parms$cpIL
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  dRP <- -decR/cpR  + parms$iR/parms$cpIR  + tvrP
  # here plant uptake as absolute parameter
  dI <-  +parms$iI  - parms$kIPlant  - leach  + PhiU  + PhiB  + PhiTvr
  dIP <-  +parms$iIP  - parms$kIPPlant  - leachP  + PhiPU  + PhiPB  + PhiPTvr
  #
  if (isTRUE(parms$isFixedS)) {
    # scenario of fixed substrate
    dR <- dL <- dRN <- dLN <- dRP <- dLP <- dI <- dIP <- 0
  } else if (isTRUE(parms$isTvrNil)) {
    # scenario of enzymes and biomass not feeding back to R
    dR <- +parms$iR - decR
    dRN <- +parms$iR/parms$cnIR - decR/cnR
    dRP <- +parms$iR/parms$cpIR - decR/cpR
    tvrExC <- tvrC
    tvrExN <- tvrN
    tvrExP <- tvrP
  }
  #
  resDeriv <- structure(as.numeric(
    c( dB, dR, dRN, dRP, dL, dLN, dLP, dI, dIP, dAlpha))
    ,names = c("dB","dR","dRN","dRP","dL","dLN","dLP","dI","dIP","dAlpha"))
  if (any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
  sqrEps <- sqrt(.Machine$double.eps)
  # parms$iL - (decL + dL)
  # parms$iR + tvrC -(decR + dR)
  #
  # checking the mass balance of fluxes
  plantNUp <- plantPUp <- 0 # checked in mass balance but is not (any more) in model
  respB <- (synE)/parms$eps*(1 - parms$eps)  + rG + rM + respO
  resp <- respB + respTvr
  if (diff( unlist(c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps )  stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(uN = uNOrg, usage = synE/parms$cnE + synB/parms$cnB + PhiB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (diff( unlist(
    c(uP = uPOrg, usage = synE/parms$cpE + synB/parms$cpB + PhiPB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dR + dL + tvrExC + resp,    parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + dRN + dLN + dI + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL - plantNUp  + parms$iI -
         parms$kIPlant - parms$l*x["I"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
    if (diff(unlist(
      c( dB/parms$cpB  + dRP + dLP + dIP + tvrExP
         , parms$iR/parms$cpIR  + parms$iL/parms$cpIL - plantPUp  + parms$iIP -
         parms$kIPPlant - parms$lP*x["IP"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) { resDeriv["dR"] <- resDeriv["dRN"] <- resDeriv["dRP"] <- 0 }
  if (isTRUE(parms$isFixedL)) { resDeriv["dL"] <- resDeriv["dLN"] <- resDeriv["dLP"] <- 0 }
  if (isTRUE(parms$isFixedI)) { resDeriv["dI"] <-  resDeriv["dIP"] <- 0   }
  #
  # further computations just for output for tacking the system
  ER <- alpha * parms$aE * x["B"] / parms$kN
  EL <- (1 - alpha) * parms$aE * x["B"] / parms$kN
  limER <- ER / (parms$kmR + ER)
  limEL <- EL / (parms$kmL + EL)
  revRC <- dRPot / (parms$km*parms$kN + alphaC*aeB)
  revLC <- dLPot / (parms$km*parms$kN + (1 - alphaC)*aeB)
  revRN <- dRPot/cnR / (parms$km*parms$kN + alphaN*aeB)
  revLN <- dLPot/cnL / (parms$km*parms$kN + alphaN*aeB)
  # net mic mineralization/immobilization when accounting uptake mineralization
  PhiBU <- PhiB + PhiU
  # total mineralization flux including microbial turnover
  PhiTotal <- PhiBU + PhiTvr
  # do not match in other limitation
  # c(alphaC, revRC/(revRC + revLC)); c(alphaN, revRN/(revRN + revLN))
  # compute C available for biomass, this time without accounting immobilization flux
  NsynBNSubstrate <- uNOrg - synE/cnE
  CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
  # N limited based on substrate uptake (without accounting for immobilization)
  isLimNSubstrate <-  CsynBNSubstrate < CsynBC
  #
  if (isTRUE(parms$isRecover) ) recover()
  list( resDeriv,  c(
     resp = as.numeric(resp)
     , respO = as.numeric(respO)
     , respB = as.numeric(respB)
     , respTvr = as.numeric(respTvr)
     #, ER = as.numeric(ER), EL = as.numeric(EL)
    #, MmB = as.numeric(MmB)
    , PhiTotal = as.numeric(PhiTotal)
    , PhiB = as.numeric(PhiB), PhiU = as.numeric(PhiU)
    , PhiTvr = as.numeric(PhiTvr)
    , PhiBU = as.numeric(PhiBU)
    , immoPot = as.numeric(immoPot)
    , PhiPTotal = as.numeric(PhiPB + PhiPU + PhiTvr)
    , PhiPB = as.numeric(PhiPB), PhiPU = as.numeric(PhiPU)
    , PhiPTvr = as.numeric(PhiPTvr)
    , PhiPBU = as.numeric(PhiPB + PhiPU)
    , immoPPot = as.numeric(immoPPot)
    , resBalance$wELim
    , alphaTarget = as.numeric(alphaTarget)
    , alphaC = as.numeric(alphaC), alphaN = as.numeric(alphaN), alphaP = as.numeric(alphaP)
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , cpR = as.numeric(cpR), cpL = as.numeric(cpL)
    , limER = as.numeric(limER), limEL = as.numeric(limEL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , tvrB = as.numeric(tvrB)
    , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    , revRN = as.numeric(revRN)
    , revLN = as.numeric(revLN)
    , CsynB = as.numeric(CsynB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    #, pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    #, NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    #, dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    #, dI = as.numeric(dI)
    #, uC = as.numeric(uC), synB = as.numeric(synB)
    #, decN = as.numeric(decN)
  ))
}

.depr.balanceAlphaBetweenElementLimitations <- function(
  ### compute balance between alphas of different element limitations
  alpha    ##<< numeric vector of allocation coefficients for different elements
  , CsynBE ##<< numeric vector of carbon availale for biomass synthesis
  ## for each based on element limitaiton, first is carbon
  , ce     ##<< numeric vector for C to E microbial biomass ratios
  ## The first entry C to C is not considered.
  , eps    ##<< numeric scalar: intrinsic carbon use eff. for biomass synthesis
  , delta = 200  ##<< scalar smoothing factor
  #, alphaC, alphaN, CsynBN, CsynBC, NsynBC, NsynBN
){
  ##details<< Select the alpha corresponding to the smallest CsynBE.
  ## However, if elemental limitations are close,
  ## do a smooth transition between corresponding alpha values
  wCLim = min( .Machine$double.xmax, (min(CsynBE[-1]/CsynBE[1]))^delta )
  NSynBE <- eps*CsynBE[2]/ce[2]
  wNLim = min( .Machine$double.xmax, (min(CsynBE[-2]/CsynBE[2]))^delta )
  ce[1] <- 1/eps  # so that EsynBE equals CSynBE for first component
  wELim <- sapply( seq_along(CsynBE), function(iE){
    EsynBE <- eps*CsynBE/ce[iE]  # available mass of E for biomass synthesis
    # min to avoid  + Inf
    wELimi <- min( .Machine$double.xmax, (min(EsynBE[-iE]/EsynBE[iE]))^delta )
  })
  #alpha <- (wCLim*alphaC + wNLim*alphaN) / (wCLim + wNLim)
  alphaBalanced <- sum(wELim * alpha)/sum(wELim)
  alphaBalanced
}

.depr.balanceAlphaBetweenElementLimitations <- function(
  ### compute balance between alphas of different element limitations
  alpha    ##<< numeric vector of allocation coefficients for different elements
  , CsynBE ##<< numeric vector of carbon availale for biomass synthesis
  , delta = 200  ##<< scalar smoothing factor
){
  ##details<< Select the alpha corresponding to the smallest CsynBE.
  ## However, if elemental limitations are close,
  ## do a smooth transition between corresponding alpha values
  wELim <- sapply( seq_along(CsynBE), function(iE){
    # min to avoid  + Inf
    wELimi <- min( .Machine$double.xmax, (min(CsynBE[-iE]/CsynBE[iE]))^delta )
  })
  alphaBalanced <- sum(wELim * alpha)/sum(wELim)
  alphaBalanced
}

balanceAlphaBetweenElementLimitations <- function(
  ### compute balance between alphas of different element limitations
  alpha    ##<< numeric vector of allocation coefficients for different elements
  , CsynBE ##<< numeric vector of carbon availale for biomass synthesis
  , tauB   ##<< scakar typical microbial turnover flux for scaling
  , delta = 10  ##<< scalar smoothing factor, the higher, the steeper the transition
){
  ##details<< Select the alpha corresponding to the smallest CsynBE.
  ## However, if elemental limitations are close,
  ## do a smooth transition between corresponding smallest alpha values
  wELim <- sapply( seq_along(CsynBE), function(iE){
    wELimi <- min( .Machine$double.xmax
                   , exp( delta/tauB*(min(CsynBE[-iE]) - CsynBE[iE])))
  })
  names(wELim) <- names(alpha)
  wELimNorm <- wELim/sum(wELim)
  alphaBalanced <- sum(wELimNorm * alpha)
  ##value<< a list with entries
  list(
    alpha = alphaBalanced  ##<< numeric scalar: balanced alpha
    , wELim = wELimNorm    ##<< numeric vector: fractional element limitations
  )
}


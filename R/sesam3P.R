#library(deSolve)

# gC/m2 and gN/m2, /yr

#' @export
derivSesam3P <- function(
  ### Soil Enzyme Steady Allocation model including phosporous dynamics
  t,x,parms
){
  ##details<<
  ## Sesam3P (Sesam3a extended by Phosphorous) extended by biomineralization
  x <- pmax(unlist(x),1e-16)      # no negative masses
  # compute steady state enzyme levels for N and for C limitation
  dRPot <- parms$kR * x["R"]
  dLPot <- parms$kL * x["L"]
  dLPPot <- parms$kLP * x["LP"] # potential biomineralization
  dRPPot <- parms$kRP * x["RP"]
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
  # one state variable less, here LP, because alphas sum to one
  alpha <- cbind(L = x["alphaL"], R = x["alphaR"], LP = NA, RP = x["alphaRP"])[1,]
  alpha["LP"] <- 1 - sum(alpha, na.rm = TRUE)
  B <- x["B"]
  aeB <- parms$aE*B        # aeB without associated growth respiration
  kmN <- parms$kmN #parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  tvrB <- parms$tau*B      # microbial turnover
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * alpha["L"]*aeB/(kmN + alpha["L"]*aeB)
  decR <- dRPot * alpha["R"]*aeB/(kmN + alpha["R"]*aeB)
  decLP <- dLPPot * alpha["LP"]*aeB/(kmN + alpha["LP"]*aeB)
  decRP <- dRPPot * alpha["RP"]*aeB/(kmN + alpha["RP"]*aeB)
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
  # community composition and enzyme allocation
  limE <- computeElementLimitations(
    cbind(C = CsynBC, N = CsynBN, P = CsynBP)[1,]
    , tauB = parms$tau*B)
  alphaTarget <- computeSesam3PAllocationPartitioning(
    dS = cbind(R = dRPot, L = dLPot)[1,]
    ,dSP = cbind(R = dRPPot, L = dLPPot)[1,]
    , B = B
    ,kmkN = kmN, aE =  parms$aE
    ,alpha = alpha
    ,limE = limE
    ,beta = cbind(L = cnL, R = cpR, E = parms$cnE)[1,]
    ,gamma = cbind(L = cpL, R = cpR, E = parms$cpE)[1,]
  )
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
  plantNUpPot <- parms$kIPlant*x["I"]
  plantNUp <- if (!is.null(parms$plantNUpAbs)) {
    min(parms$plantNUpAbs, plantNUpPot)
  } else {
    plantNUpPot
  }
  plantPUpPot <- parms$kIPPlant*x["IP"]
  plantPUp <- if (!is.null(parms$plantPUpAbs)) {
    min(parms$plantPUpAbs, plantPUpPot)
  } else {
    plantPUpPot
  }
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
  dI <-  +parms$iI  - plantNUp  - leach  + PhiU  + PhiB  + PhiTvr
  dIP <-  +parms$iIP  - plantPUp  - leachP  + PhiPU  + PhiPB  + PhiPTvr
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
    c( dB, dR, dRN, dRP, dL, dLN, dLP, dI, dIP
       , dAlpha["L"], dAlpha["R"], dAlpha["RP"]))
    ,names = c("dB","dR","dRN","dRP","dL","dLN","dLP","dI","dIP"
               ,"dAlphaL","dAlphaR","dAlphaRP"))
  if (any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
  sqrEps <- sqrt(.Machine$double.eps)
  # parms$iL - (decL + dL)
  # parms$iR + tvrC -(decR + dR)
  #
  # checking the mass balance of fluxes
  respB <- (synE)/parms$eps*(1 - parms$eps)  + rG + rM + respO
  resp <- respB + respTvr
  if (diff( unlist(c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps )  stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(uN = uNOrg, usage = synE/parms$cnE + synB/parms$cnB + PhiB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (diff( unlist(
    c(uP = uPOrg, usage = synE/parms$cpE + synB/parms$cpB + PhiPB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance P error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dR + dL + tvrExC + resp,    parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + dRN + dLN + dI + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL + parms$iI -
         plantNUp - parms$l*x["I"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
    if (diff(unlist(
      c( dB/parms$cpB  + dRP + dLP + dIP + tvrExP
         , parms$iR/parms$cpIR  + parms$iL/parms$cpIL + parms$iIP -
         plantPUp - parms$lP*x["IP"])))^2 >
      .Machine$double.eps )  stop("mass balance dP error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) { resDeriv["dR"] <- resDeriv["dRN"] <- resDeriv["dRP"] <- 0 }
  if (isTRUE(parms$isFixedL)) { resDeriv["dL"] <- resDeriv["dLN"] <- resDeriv["dLP"] <- 0 }
  if (isTRUE(parms$isFixedI)) { resDeriv["dI"] <-  resDeriv["dIP"] <- 0   }
  #
  # further computations just for output for tacking the system
  limZ <- alpha * aeB / (parms$kmN + alpha * aeB)
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
    , structure(limE, names = paste0("lim",names(limE)))
    , structure(alphaTarget, names = paste0("alpha",names(alphaTarget)))
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , cpR = as.numeric(cpR), cpL = as.numeric(cpL)
    , structure(limZ, names = paste0("limZ",names(limZ)))
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , tvrB = as.numeric(tvrB)
    , synB = as.numeric(synB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    , CsynBP = as.numeric(CsynBP)
    #, pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    #, NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    #, dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    #, dI = as.numeric(dI)
    #, uC = as.numeric(uC), synB = as.numeric(synB)
    #, decN = as.numeric(decN)
  ))
}

computeElementLimitations <- function(
  ### compute element limitations
  CsynBE   ##<< numeric vector of carbon available for biomass synthesis
  , tauB   ##<< scalar typical microbial turnover flux for scaling
  , delta = 10  ##<< scalar smoothing factor, the higher, the steeper the transition
){
  ##details<< Select the alpha corresponding to the smallest CsynBE.
  ## However, if elemental limitations are close,
  ## do a smooth transition between corresponding smallest alpha values
  wELim <- sapply( seq_along(CsynBE), function(iE){
    wELimi <- min( .Machine$double.xmax
                   , exp( delta/tauB*(min(CsynBE[-iE]) - CsynBE[iE])))
  })
  names(wELim) <- names(CsynBE)
  wELimNorm <- wELim/sum(wELim)
  ##value<< numeric vector: fractional element limitations with names
  ## corresponding to names of CsynBE
}

#' @export
computeSesam3PAllocationPartitioning <- function(
  ###  allocation partitioning alpha of four enzymes and biomineralization
  dS    ##<< numeric vector (L,R) of potential depolymerization C-fluxes
  ,dSP  ##<< numeric vector of potential biomineralizationdepolymerization P-fluxes
  ,B		##<< numeric vector of microbial biomass carbon
  ,kmkN	##<< numeric vector of product of (half-saturation constant in
  ## decomposition equation) x (enzyme turnover rate)
  ,aE		##<< numeric vector of proportion of microbial biomass allocated
  ## to enzyme production per time
  ,alpha  ##<< numeric vector (L,R,LP,RP) of current community allocation
  ## coefficients in (0,1)
  ,limE   ##<< numeric vector (nElement) of weights of elemental limitations
  ,beta    ##<< numeric vector of C/N ratios of substrates and enzymes
  ,gamma   ##<< numeric vector of C/P ratios of substrates (L,R) and enzymes (E)
){
  synEnz <- aE*B
  cost <- (kmkN + alpha*synEnz) *
    (limE["C"] + limE["N"]/beta["E"] + limE["P"]/gamma["E"])
  returnL <- dS["L"]*(limE["C"] + limE["N"]/beta["L"] + limE["P"]/gamma["L"])
  returnR <- dS["R"]*(limE["C"] + limE["N"]/beta["R"] + limE["P"]/gamma["R"])
  returnLP <- dSP["L"]*(limE["P"])
  returnRP <- dSP["R"]*(limE["P"])
  returnE <- structure(c(returnL, returnR, returnLP, returnRP),
                       names = c("L","R","LP","RP"))
  revenue <- returnE/cost[names(returnE)]
  alpha <- revenue/sum(revenue)
  alpha
  ##value<< numeric vector of revenue-optimal alpha for enzymes L,R,LP,RP
}



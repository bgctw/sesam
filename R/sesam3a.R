#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam3a <- function(
  ### Soil Enzyme Steady Allocation model
  t,x,parms
){
  ##details<<
  ## Simplified version of sesam3s (kmN instead of km and kN)
  ## model corresponding to Seam3 with enzyme levels computed by quasi steady state
  ## Alpha as an explicit state variable that changes with turnover
  ## Simplified computation of target based on current revenue based on current alpha
  if (!length(names(x))) names(x) <- c(
    "B", "R", "RN", "L", "LN", "IN", "alpha")
  x <- pmax(unlist(x),1e-16)      # no negative masses
  # compute steady state enzyme levels for N and for C limitation
  dRPot <- parms$kR * x["R"]
  dLPot <- parms$kL * x["L"]
  immoNPot <- parms$iBN * x["IN"]
  cnR <- x["R"]/x["RN"]
  cnL <- x["L"]/x["LN"]
  cnE <- parms$cnE
  cnB <- parms$cnB
  alpha <- x["alpha"]
  B <- x["B"]
  aeB <- parms$aE*B        # aeB without associanted growth respiration
  kmN <- parms$kmN #parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  tvrB <- parms$tau*B      # microbial turnover
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * (1 - alpha)*aeB/(kmN + (1 - alpha)*aeB)
  decR <- dRPot * alpha*aeB/(kmN + alpha*aeB)
  #
  tvrERecycling <- parms$kNB*synE
  uNOrg <- parms$nuN*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  uC <- decL + decR + tvrERecycling
  CsynBCt <- uC - rM - synE/parms$eps
  CsynBC <- if(CsynBCt > 0) parms$eps*CsynBCt else CsynBCt
  NsynBN <- (uNOrg + immoNPot - synE/cnE)
  CsynBN <- cnB*NsynBN
  synB <- min(CsynBC, CsynBN)
  rG <- if (synB > 0) (1 - parms$eps)/parms$eps*synB else 0
  PhiNB <- uNOrg - synB/cnB - synE/cnE
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
  alphaTarget <- balanceAlphaBetweenCNLimitationsExp(
    alphaC, alphaN, CsynBN, CsynBC, tauB = parms$tau*B  )
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) *  (parms$tau + abs(synB)/B)
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC - (synE/parms$eps + synB + rG + rM)
  respTvr <- (1 - parms$epsTvr) * tvrB
  # assuming same cnRatio of predators to equal cn ratio of microbes
  PhiNTvr <- respTvr/parms$cnB
  #
  # tvr that feeds R pool, assume that N in SOM for resp (by epsTvr) is mineralized
  tvrC <-  +parms$epsTvr*tvrB   + (1 - parms$kNB)*synE
  tvrN <-  +parms$epsTvr*tvrB/parms$cnB   + (1 - parms$kNB)*synE/parms$cnE
  # fluxes leaving the system (will be set in scen where trv does not feed back)
  tvrExC <- tvrExN <- 0
  #
  leachN <- parms$lN*x["IN"]
  plantNUpPot <- parms$kINPlant*x["IN"]
  plantNUp <- if (!is.null(parms$plantNUpAbs)) {
    min(parms$plantNUpAbs, plantNUpPot)
  } else {
    plantNUpPot
  }
  PhiNU <- (1 - parms$nuN)*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  #
  dB <- synB - tvrB
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  # here plant uptake as absolute parameter
  dIN <-  +parms$iIN  - plantNUp  - leachN  + PhiNU  + PhiNB  + PhiNTvr
  #
  if (isTRUE(parms$isFixedS)) {
    # scenario of fixed substrate
    dR <- dL <- dRN <- dLN <- dIN <- 0
  } else if (isTRUE(parms$isTvrNil)) {
    # scenario of enzymes and biomass not feeding back to R
    dR <- +parms$iR - decR
    dRN <- +parms$iR/parms$cnIR - decR/cnR
    tvrExC <- tvrC
    tvrExN <- tvrN
  }
  #
  resDeriv <- structure(as.numeric(
    c( dB, dR, dRN, dL, dLN, dIN, dAlpha))
    ,names = c("dB","dR","dRN","dL","dLN","dIN","dAlpha"))
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
    c(uN = uNOrg, usage = synE/parms$cnE + synB/parms$cnB + PhiNB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dR + dL + tvrExC + resp,    parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + dRN + dLN + dIN + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL + parms$iIN -
         plantNUp - parms$lN*x["IN"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) { resDeriv["dR"] <- resDeriv["dRN"] <-  0   }
  if (isTRUE(parms$isFixedL)) { resDeriv["dL"] <- resDeriv["dLN"] <-  0   }
  if (isTRUE(parms$isFixedI)) { resDeriv["dIN"] <-  0   }
  if (isTRUE(parms$isFixedAlpha)) { resDeriv["dAlpha"] <-  0   }
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
  PhiNBU <- PhiNB + PhiNU
  # total mineralization flux including microbial turnover
  PhiNTotal <- PhiNBU + PhiNTvr
  # do not match in other limitation
  # c(alphaC, revRC/(revRC + revLC)); c(alphaN, revRN/(revRN + revLN))
  # compute C available for biomass, this time without accounting immobilization flux
  NsynBNSubstrate <- uNOrg - synE/cnE
  CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
  # N limited based on substrate uptake (without accounting for immobilization)
  isLimNSubstrate <-  CsynBNSubstrate < CsynBC
  #
  if (isTRUE(parms$isRecover) ) recover()
  list( resDeriv, c(
     resp = as.numeric(resp)
     , respO = as.numeric(respO)
     , respB = as.numeric(respB)
     , respTvr = as.numeric(respTvr)
     #, ER = as.numeric(ER), EL = as.numeric(EL)
    #, MmB = as.numeric(MmB)
    , PhiNTotal = as.numeric(PhiNTotal)
    , PhiNB = as.numeric(PhiNB)
    , PhiNU = as.numeric(PhiNU)
    , PhiNTvr = as.numeric(PhiNTvr)
    , PhiNBU = as.numeric(PhiNBU)
    , immoNPot = as.numeric(immoNPot)
    , alphaTarget = as.numeric(alphaTarget)
    , alphaC = as.numeric(alphaC), alphaN = as.numeric(alphaN)
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , limER = as.numeric(limER), limEL = as.numeric(limEL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , tvrB = as.numeric(tvrB)
    , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    , revRN = as.numeric(revRN)
    , revLN = as.numeric(revLN)
    , synB = as.numeric(synB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    #, pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    #, NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    #, dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    #, dIN = as.numeric(dIN)
    , uC = as.numeric(uC)
    #, decN = as.numeric(decN)
    , plantNUp = plantNUp
  ))
}

computeOutputsSesam3 <- function(
  ### compute additional outputs based on state
  x      ##<< matrix (nTime x (1+nStateVar) as returned by lsoda
  ## first column is time
  , parms ##<< parameters used
){
  x[] <- pmax(unlist(x),1e-16)      # no negative masses
  # compute steady state enzyme levels for N and for C limitation
  dRPot <- parms$kR * x[,"R"]
  dLPot <- parms$kL * x[,"L"]
  immoNPot <- parms$iBN * x[,"IN"]
  cnR <- x[,"R"]/x[,"RN"]
  cnL <- x[,"L"]/x[,"LN"]
  cnE <- parms$cnE
  cnB <- parms$cnB
  alpha <- x[,"alpha"]
  B <- x[,"B"]
  aeB <- parms$aE*B        # aeB without associanted growth respiration
  kmN <- parms$kmN #parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  tvrB <- parms$tau*B      # microbial turnover
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * (1 - alpha)*aeB/(kmN + (1 - alpha)*aeB)
  decR <- dRPot * alpha*aeB/(kmN + alpha*aeB)
  #
  tvrERecycling <- parms$kNB*synE
  uNOrg <- parms$nuN*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  uC <- decL + decR + tvrERecycling
  CsynBC <- uC - rM - synE/parms$eps
  NsynBN <- (uNOrg + immoNPot - synE/cnE)
  CsynBN <- ifelse(NsynBN > 0, cnB/parms$eps*NsynBN, cnB*NsynBN)
  synB <- pmin(CsynBC, CsynBN)
  synB <- ifelse(synB > 0, parms$eps*synB, synB )
  rG <- ifelse(synB > 0, (1 - parms$eps)*synB, 0)
  PhiNB <- uNOrg - synB/cnB - synE/cnE
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
  alphaTarget <- balanceAlphaBetweenCNLimitationsExp(
    alphaC, alphaN, CsynBN, CsynBC, tauB = parms$tau*B  )
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) *  (parms$tau + abs(synB)/B)
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC - (synE/parms$eps + synB + rG + rM)
  respTvr <- (1 - parms$epsTvr) * tvrB
  # assuming same cnRatio of predators to equal cn ratio of microbes
  PhiNTvr <- respTvr/parms$cnB
  #
  # tvr that feeds R pool, assume that N in SOM for resp (by epsTvr) is mineralized
  tvrC <-  +parms$epsTvr*tvrB   + (1 - parms$kNB)*synE
  tvrN <-  +parms$epsTvr*tvrB/parms$cnB   + (1 - parms$kNB)*synE/parms$cnE
  # fluxes leaving the system (will be set in scen where trv does not feed back)
  tvrExC <- tvrExN <- 0
  #
  leachN <- parms$lN*x[,"IN"]
  plantNUpPot <- parms$kINPlant*x[,"IN"]
  if (!is.null(parms$plantNUpAbs)) {
    plantNUp <- min(plantNUpPot, parms$plantNUpAbs)
  } else {
    plantNUp <- plantNUpPot
  }
  PhiNU <- (1 - parms$nuN)*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  #
  dB <- synB - tvrB
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  # here plant uptake as absolute parameter
  dIN <-  +parms$iIN  - plantNUp  - leachN  + PhiNU  + PhiNB  + PhiNTvr
  #
  if (isTRUE(parms$isFixedS)) {
    # scenario of fixed substrate
    dR[] <- dL[] <- dRN[] <- dLN[] <- dIN[] <- 0
  } else if (isTRUE(parms$isTvrNil)) {
    # scenario of enzymes and biomass not feeding back to R
    dR <- +parms$iR - decR
    dRN <- +parms$iR/parms$cnIR - decR/cnR
    tvrExC <- tvrC
    tvrExN <- tvrN
  }
  #
  resDeriv <- cbind( dB, dR, dRN, dL, dLN, dIN, dAlpha)
  colnames(resDeriv) <- c("dB","dR","dRN","dL","dLN","dIN","dAlpha")
  if (any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
  sqrEps <- sqrt(.Machine$double.eps)
  # parms$iL - (decL + dL)
  # parms$iR + tvrC -(decR + dR)
  #
  # checking the mass balance of fluxes
  plantNUp <- 0 # checked in mass balance but is not (any more) in model
  respB <- (synE)/parms$eps*(1 - parms$eps)  + rG + rM + respO
  resp <- respB + respTvr
  if (any( (uC - (respB + synB + synE ))^2 > sqrEps ))  stop(
    "biomass mass balance C error")
  if (any(
    (uNOrg - (synE/parms$cnE + synB/parms$cnB + PhiNB))^2 >
    .Machine$double.eps))  stop("biomass mass balance N error")
  if (!isTRUE(parms$isFixedS)) {
    if (any(
      ((dB + dR + dL + tvrExC + resp) - (parms$iR + parms$iL))^2 >
      sqrEps ))  stop("mass balance C error")
    if (any(
      ((dB/parms$cnB  + dRN + dLN + dIN + tvrExN) -
         (parms$iR/parms$cnIR  + parms$iL/parms$cnIL - plantNUp  + parms$iIN -
         plantNUp - parms$lN*x[,"IN"]))^2 >
      .Machine$double.eps ))  stop("mass balance dN error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) { resDeriv[,"dR"] <- resDeriv[,"dRN"] <-  0   }
  if (isTRUE(parms$isFixedL)) { resDeriv[,"dL"] <- resDeriv[,"dLN"] <-  0   }
  if (isTRUE(parms$isFixedI)) { resDeriv[,"dIN"] <-  0   }
  if (isTRUE(parms$isFixedAlpha)) {
    resDeriv[,"dAlpha"] <-  0
    }
  #
  # further computations just for output for tacking the system
  ER <- alpha * parms$aE * x[,"B"] / parms$kN
  EL <- (1 - alpha) * parms$aE * x[,"B"] / parms$kN
  limER <- ER / (parms$km + ER)
  limEL <- EL / (parms$km + EL)
  revRC <- dRPot / (parms$km*parms$kN + alphaC*aeB)
  revLC <- dLPot / (parms$km*parms$kN + (1 - alphaC)*aeB)
  revRN <- dRPot/cnR / (parms$kmN + alphaN*aeB)
  revLN <- dLPot/cnL / (parms$kmN + alphaN*aeB)
  # net mic mineralization/immobilization when accounting uptake mineralization
  PhiNBU <- PhiNB + PhiNU
  # total mineralization flux including microbial turnover
  PhiNTotal <- PhiNBU + PhiNTvr
  # do not match in other limitation
  # c(alphaC, revRC/(revRC + revLC)); c(alphaN, revRN/(revRN + revLN))
  # compute C available for biomass, this time without accounting immobilization flux
  NsynBNSubstrate <- uNOrg - synE/cnE
  CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
  # N limited based on substrate uptake (without accounting for immobilization)
  isLimNSubstrate <-  CsynBNSubstrate < CsynBC
  rownames(resDeriv) <- rownames(x) <- NULL
  cbind( x, resDeriv, data.frame(
    resp = as.numeric(resp)
    , respO = as.numeric(respO)
    , respB = as.numeric(respB)
    , respTvr = as.numeric(respTvr)
    #, ER = as.numeric(ER), EL = as.numeric(EL)
    #, MmB = as.numeric(MmB)
    , PhiNTotal = as.numeric(PhiNTotal)
    , PhiNB = as.numeric(PhiNB)
    , PhiNU = as.numeric(PhiNU)
    , PhiNTvr = as.numeric(PhiNTvr)
    , PhiNBU = as.numeric(PhiNBU)
    , immoNPot = as.numeric(immoNPot)
    , alphaTarget = as.numeric(alphaTarget)
    , alphaC = as.numeric(alphaC), alphaN = as.numeric(alphaN)
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , limER = as.numeric(limER), limEL = as.numeric(limEL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , tvrB = as.numeric(tvrB)
    , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    , revRN = as.numeric(revRN)
    , revLN = as.numeric(revLN)
    , synB = as.numeric(synB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    , plantNUp = plantNUp # inorganic N taken up by plant
    #, pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    #, NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    #, dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    #, dIN = as.numeric(dIN)
    #, uC = as.numeric(uC), synB = as.numeric(synB)
    #, decN = as.numeric(decN)
  ))
}

computeSesam3sAllocationPartitioning <- function(
  ###  allocation partitioning alpha assuming enzymes to be in quasi steady state
  dR		##<< numeric vector of potential decomposition flux from
  ## residue pool (gC or gN/time)
  ,dL		##<< numeric vector of potential decomposition flux from labile pool
  ,B		##<< numeric vector of microbial biomass carbon
  ,kmkN	##<< numeric vector of product of (half-saturation constant in
  ## decomposition equation) x (enzyme turnover rate)
  ,aE		##<< numeric vector of proportion of microbial biomass allocated
  ## to enzyme production per time
  ,alpha  ##<< numeric vector of current community allocation coefficient in (0,1)
){
  revL = dL/(kmkN + (1 - alpha)*aE*B)
  revR = dR/(kmkN + alpha*aE*B)
  alphaTarget = revR/(revL + revR)
  structure(alphaTarget, names = rep("alpha",length(alphaTarget)))
  ##value<< sacalar numeric revenue-optimal alpha value
}
attr(computeSesam3sAllocationPartitioning,"ex") <- function(){
  parms <- list(
    cnB = 11
    ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
    #,cnIR = 4.5   ##<< between micr and enzyme signal
    ,cnIR = 7      ##<< between micr and enzyme signal
    ,cnIL = 30     ##<< N poor substrate, here no inorganic N supply,
    # need to have low C/N for CO2 increase scenario
    ,kN = 60       ##<< /yr enzyme turnover 60 times a year, each 6 days ->
    # fast priming pulses
    ,km = 0.05     ##<< enzyme half-saturation constant, in magnitude of enzymes,
    # determined by kN
    ,kR = 1/(10)   ##<< 1/(x years) # to demonstrate changes on short time scale
    ,kL = 5        ##<< 1/(x years) # formerly 1 year
    ,aE = 0.001*365   ##<< C biomass allocated to enzymes 1/day /microbial biomass
  )
  x <- c( #aE = 0.001*365
    B = 17                     ##<< microbial biomass
    ,ER  = 2*parms$km                  ##<< total enzyme pool
    ,EL  = 4*parms$km                   ##<< total enzyme pool
    ,R = 1000                ##<< N rich substrate
    ,RN = 1000/parms$cnIR   ##<< N rich substrate N pool
    ,L = 100                 ##<< N poor substrate
    ,LN = 100/parms$cnIL    ##<< N poor substrate N pool
    ,IN =  0.4                ##<< inorganic pool gN/m2
  )
  B <- x["B"]
  # allocate only few to R (and most to L)
  computeSesam3sAllocationPartitioning(
    dR = parms$kR*x["R"], dL = parms$kL*x["L"], B = B
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    , alpha = 0.5
  )
  # if decL is small, allocate much to R
  computeSesam3sAllocationPartitioning(
    dR = parms$kR*x["R"], dL = parms$kL*x["L"]/500, B = B
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    , alpha = 0.5
  )
  # equal partitioning
  computeSesam3sAllocationPartitioning(
    dR = 1, dL = 1, B = B
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    , alpha = 0.5
  )
  # allocate 2/3 to R
  computeSesam3sAllocationPartitioning(
    dR = 2, dL = 1, B = B
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    , alpha = 0.5
  )
  # with missing dL, allocate all enzymes to R
  computeSesam3sAllocationPartitioning(
    dR = 2, dL = 0, B = B
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    , alpha = 0.5
  )
  # with missing dR, allocate zero enzymes to R
  computeSesam3sAllocationPartitioning(
    dR = 0, dL = 1, B = B
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    , alpha = 0.5
  )
  # with no decomposition, return NaN
  computeSesam3sAllocationPartitioning(
    dR = 0, dL = 0, B = B
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    , alpha = 0.5
  )
}


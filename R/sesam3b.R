#library(deSolve)

# gC/m2 and gN/m2, /yr

# derivSesam3b conflicts with derivSesam3B in inlinedocs,
derivSesam3b_ <- function(
  ### Soil Enzyme Steady Allocation model
  t,x,parms
){
  ##details<<
  ## Modification of derivSesam3a
  ## with limitation-weighted cost and multidimensional alpha
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
  alpha <- cbind(R = x["alphaR"], L = 1 - x["alphaR"])[1,]
  B <- x["B"]
  aeB <- parms$aE*B        # aeB without associanted growth respiration
  kmN <- parms$kmN #parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  tvrB <- parms$tau*B      # microbial turnover
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * (alpha["L"])*aeB/(kmN + (alpha["L"])*aeB)
  decR <- dRPot * alpha["R"]*aeB/(kmN + alpha["R"]*aeB)
  #
  tvrERecycling <- parms$kNB*synE
  uNOrg <- parms$nuN*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  uC <- decL + decR + tvrERecycling
  CsynBC <- uC - rM - synE/parms$eps
  NsynBN <- (uNOrg + immoNPot - synE/cnE)
  CsynBN <- if (NsynBN > 0) cnB/parms$eps*NsynBN else cnB*NsynBN
  CsynB <- min(CsynBC, CsynBN)
  if (CsynB > 0) {
    synB <- parms$eps*CsynB
    rG <- (1 - parms$eps)*CsynB
  } else {
    synB <- CsynB # with negative biomass change, do growth respiration
    rG <- 0
  }
  PhiNB <- uNOrg - synB/cnB - synE/cnE
  # alphaCR <- computeSesam3sAllocationPartitioning(
  #   dR = dRPot, dL = dLPot, B = B
  #   , kmkN = kmN, aE = parms$aE
  #   , alpha = alpha["R"]
  # )
  # alphaN <- computeSesam3sAllocationPartitioning(
  #   dR = dRPot/cnR, dL = dLPot/cnL, B = B
  #   , kmkN = kmN, aE = parms$aE
  #   , alpha = alpha["R"]
  # )
  # alphaTarget <- balanceAlphaBetweenCNLimitationsExp(
  #   alphaCR, alphaN, CsynBN, CsynBC, tauB = parms$tau*B  )
  wELim <- computeElementLimitations(
    cbind(C = CsynBC, N = CsynBN)[1,]
    , tauB = parms$tau*B
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
    )
  alphaTarget <- computeSesam3bAllocationPartitioning(
    dS = cbind(R = dRPot, L = dLPot)[1,]
    , B = B
    ,kmkN = kmN, aE =  parms$aE
    ,alpha = alpha
    ,wELim = wELim
    ,betaN = cbind(L = cnL, R = cnR, E = parms$cnE)[1,]
  )
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
    c( dB, dR, dRN, dL, dLN, dIN, dAlpha["R"]))
    ,names = c("dB","dR","dRN","dL","dLN","dIN","dAlphaR"))
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
  if (isTRUE(parms$isFixedAlpha)) { resDeriv["dAlphaR"] <-  0   }
  #
  # further computations just for output for tacking the system
  ER <- alpha["R"] * parms$aE * x["B"] / parms$kN
  EL <- (alpha["L"]) * parms$aE * x["B"] / parms$kN
  limER <- ER / (parms$kmR + ER)
  limEL <- EL / (parms$kmL + EL)
  # revRC <- dRPot / (parms$km*parms$kN + alphaC["R"]*aeB)
  # revLC <- dLPot / (parms$km*parms$kN + (1 - alphaCR)*aeB)
  # revRN <- dRPot/cnR / (parms$km*parms$kN + alphaN*aeB)
  # revLN <- dLPot/cnL / (parms$km*parms$kN + alphaN*aeB)
  # net mic mineralization/immobilization when accounting uptake mineralization
  PhiNBU <- PhiNB + PhiNU
  # total mineralization flux including microbial turnover
  PhiNTotal <- PhiNBU + PhiNTvr
  # do not match in other limitation
  # c(alphaCR, revRC/(revRC + revLC)); c(alphaN, revRN/(revRN + revLN))
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
    #, alphaCR = as.numeric(alphaCR), alphaN = as.numeric(alphaN)
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , limER = as.numeric(limER), limEL = as.numeric(limEL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , tvrB = as.numeric(tvrB)
    # , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    # , revRN = as.numeric(revRN)
    # , revLN = as.numeric(revLN)
    , CsynB = as.numeric(CsynB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    #, pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    #, NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    #, dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    #, dIN = as.numeric(dIN)
    , uC = as.numeric(uC), synB = as.numeric(synB)
    #, decN = as.numeric(decN)
    , plantNUp = plantNUp
  ))
}

#' @export
computeSesam3bAllocationPartitioning <- function(
  ###  allocation partitioning alpha of four enzymes and biomineralization
  dS    ##<< numeric vector (L,R) of potential depolymerization C-fluxes
  ,B		##<< numeric vector of microbial biomass carbon
  ,kmkN	##<< numeric vector of product of (half-saturation constant in
  ## decomposition equation) x (enzyme turnover rate)
  ,aE		##<< numeric vector of proportion of microbial biomass allocated
  ## to enzyme production per time
  ,alpha  ##<< numeric vector (L,R) of current community allocation
  ## coefficients in (0,1)
  ,wELim   ##<< numeric vector (nElement) of weights of elemental limitations
  ,betaN    ##<< numeric vector of C/N ratios of substrates and enzymes
){
  synEnz <- aE*B
  cost <- (kmkN + alpha*synEnz) *
    (wELim["C"] + wELim["N"]/betaN["E"])
  returnL <- dS["L"]*(wELim["C"] + wELim["N"]/betaN["L"])
  returnR <- dS["R"]*(wELim["C"] + wELim["N"]/betaN["R"])
  returnE <- structure(c(returnL, returnR), names = c("L","R"))
  revenue <- returnE/cost[names(returnE)]
  alpha <- revenue/sum(revenue)
  alpha
  ##value<< numeric vector of revenue-optimal alpha for enzymes L,R
}



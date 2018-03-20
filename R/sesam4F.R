#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam4F <- function(
  ### Soil Enzyme Steady Allocation model with detailed turnover
  t,xvec,parms
){
  # based on Sesam4a with tracking different fractions of pools
  # assure no negative masses
  sqrEps <- sqrt(.Machine$double.eps)
  xOrig <- xvec
  xvec <- pmax(unlist(xvec), 0.0)      # no negative masses
  # sum pool fractions
  x <- parms$multiPoolFractions$setX(parms$multiPoolFractions, xvec)
  # elemental ratios
  cnR <- x$tot["R"]/x$tot["RN"]
  cnL <- x$tot["L"]/x$tot["LN"]
  cpR <- x$tot["R"]/x$tot["RP"]
  cpL <- x$tot["L"]/x$tot["LP"]
  cW <- parms$cW
  cnE <- parms$cnE;  cnB <- parms$cnB;  cnBW <- parms$cnBW
  cpE <- parms$cpE;  cpB <- parms$cpB;  cpBW <- parms$cpBW
  cnBL <- if (cW == 1 ) cnB else (1 - cW)/(1/cnB - cW/cnBW)
  cpBL <- if (cW == 1 ) cpB else (1 - cW)/(1/cpB - cW/cpBW)
  # abbreviations of variables from x and parms
  alpha <- x$tot["alpha"]
  B <- x$tot["B"]
  dRPot <- parms$kR * x$tot["R"]
  dLPot <- parms$kL * x$tot["L"]
  immoPot <- parms$iB * x$tot["I"]
  immoPPot <- parms$iBP * x$tot["IP"]
  aeB <- parms$aE*B
  kmN <- parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * (1 - alpha)*aeB/(kmN + (1 - alpha)*aeB)
  decR <- dRPot * alpha*aeB/(kmN + alpha*aeB)
  #
  # microbial turnover, partly feeds back to uptake
  tvrB <- parms$tau*B      # without predation/mineralization
  tvrBPred <- if (B < parms$B0) 0 else parms$tauP*(B - parms$B0)*B
  tvrBOrg <- tvrB + parms$epsP*tvrBPred # turnover feeding back to organic
  respTvr <- (1 - parms$epsP)*tvrBPred  # turnover mineralized (respired)
  PhiTvr <- respTvr/parms$cnB
  PhiPTvr <- respTvr/parms$cpB
  #
  # organic uptake
  tvrERecycling <- parms$kNB*synE # part of enzyme turnover feeds back to uptake
  decNLR <- decL/cnL + decR/cnR
  decNE <- tvrERecycling/cnE
  decNB <- tvrBOrg*(1 - cW)/cnBL
  uC <- decL + decR + tvrERecycling + tvrBOrg*(1 - cW)
  # relURecyc depends on relUC because composition of enzymes is that of uptake
  #fracUC <- decL*x$rel[["L"]] + decR*x$rel[["R"]] +
  #  tvrERecycling*relSC + tvrBOrg*(1 - cW)*x$rel[["B"]]
  fracUNOrg <- parms$nu*(decL/cnL*x$rel[["LN"]] + decR/cnR*x$rel[["RN"]] +
                           (tvrERecycling/cnE + tvrBOrg*(1 - cW)/cnBL)*x$rel[["BN"]])
  uNOrg <- parms$nu*(decL/cnL + decR/cnR + tvrERecycling/cnE + tvrBOrg*(1 - cW)/cnBL)
  #uNOrg1 <- parms$nu*(decNLR + decNE + decNB)
  #if (abs(uNOrg - uNOrg1) > 1e-14) stop("fracUNOrg error")
  uPOrg <- parms$nuP*(decL/cpL + decR/cpR +
                        tvrERecycling/cpE + tvrBOrg*(1 - cW)/cpBL)
  #uPOrg1 <- parms$nuP*(decL/cpL + decR/cpR + tvrERecycling/cpE + tvrBOrg*(1 - cW)/cpBL)
  #if (abs(uPOrg - uPOrg1) > 1e-14) stop("fracUPOrg error")
  #
  # elemental limitations by potential biomass synthesis and respiration
  CsynBC <- uC - rM - synE/parms$eps
  CsynBN <- cnB/parms$eps*(uNOrg + immoPot - synE/cnE)
  CsynBP <- cpB/parms$eps*(uPOrg + immoPPot - synE/cpE)
  CsynB <- min(CsynBC, CsynBN, CsynBP)
  if (CsynB > 0) {
    recycB <- 0
    synB <- parms$eps*CsynB
  } else {
    recycB <- -CsynB
    synB <- 0
  }
  rG <- (1 - parms$eps)*synB
  #
  # microbial community composition
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
  resBalance <- resBalance0 <- balanceAlphaBetweenElementLimitations(
    structure(c(alphaC, alphaN, alphaP), names = c("limC","limN","limP"))
    , c(CsynBC, CsynBN, CsynBP)
    , tauB = (tvrB + tvrBPred)
  )
  alphaTarget <- resBalance$alpha
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) * (synB + recycB + tvrB + tvrBPred)/B
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC + recycB - (synE/parms$eps + synB + rG + rM)
  PhiB <- uNOrg + recycB/cnB - synB/cnB - synE/cnE
  PhiPB <- uPOrg + recycB/cpB - synB/cpB - synE/cpE
  #
  tvrN <-  +cW*tvrBOrg/parms$cnBW   + (1 - parms$kNB)*synE/parms$cnE
  tvrP <-  +cW*tvrBOrg/parms$cpBW   + (1 - parms$kNB)*synE/parms$cpE
  # fluxes leaving the system (will be set in scen where trv does not feed back)
  tvrExC <- x$units$C;  tvrExN <- x$units$N;  tvrExP <- x$units$P
  tvrExC[] <- tvrExN[] <- tvrExP[] <- 0
  #
  leach <- parms$l*x$tot["I"]
  PhiU <- (1 - parms$nu)*(decL/cnL + decR/cnR + tvrERecycling/cnE + tvrBOrg*(1 - cW)/cnBL)
  leachP <- parms$lP*x$tot["IP"]
  PhiPU <- (1 - parms$nuP)*(decL/cpL + decR/cpR + tvrERecycling/cpE + tvrBOrg*(1 - cW)/cpBL)
  #
  respB <- (synE)/parms$eps*(1 - parms$eps)  + rG + rM + respO
  immoN <- max(0,-PhiB); minN <- max(0,PhiB)
  immoP <- max(0,-PhiPB); minP <- max(0,PhiPB)
  #
  # after sC and immo is known, relSC can be computed (see SteadyEnzyme4a.Rmd)
  # fluxes that feed synthesis, s, and mineralization fluxes
  sC <- uC + recycB
  sN <- uNOrg + recycB/cnB + immoN
  sP <- uPOrg + recycB/cpB + immoP
  cnS <- sC/sN
  cpS <- sC/sP
  .aC <- decL*x$rel[["L"]] + decR*x$rel[["R"]] + tvrBOrg*(1 - cW)*x$rel[["B"]] +
    recycB*x$rel[["B"]]
  .aN <- parms$nu * (decL/cnL*x$rel[["LN"]] + decR/cnR*x$rel[["RN"]] +
    tvrBOrg*(1 - cW)/cnBL*x$rel[["BN"]]) +
    recycB/cnB*x$rel[["BN"]] + immoN*x$rel[["I"]]
  .aP <- parms$nuP * (decL/cpL*x$rel[["LP"]] + decR/cpR*x$rel[["RP"]] +
    tvrBOrg*(1 - cW)/cpBL*x$rel[["BP"]]) +
    recycB/cpB*x$rel[["BP"]] * immoP*x$rel[["IP"]]
  relSC <- .aC/(sC - tvrERecycling)
  relSN <- .aN/(sN - parms$nu*tvrERecycling/cnE)
  relSP <- .aP/(sP - parms$nuP*tvrERecycling/cpE)
  relUC <- (sC*relSC - recycB*x$rel[["B"]])/uC
  #relUNOrg <- (sN*relSN - recycB/cnB*x$rel[["BN"]] - immoN*x$rel[["I"]])/uNOrg
  #relUPOrg <- (sP*relSP - recycB/cpB*x$rel[["BP"]] - immoP*x$rel[["IP"]])/uPOrg
  if (any(abs( uC*relUC - (.aC + tvrERecycling*relSC)) > sqrEps)) stop(
    "composition C mass balance error in synthesis flux")
  if (any(abs( sN*relSN - (.aN + parms$nu*tvrERecycling/cnE*relSN)) > sqrEps)) stop(
    "composition N mass balance error in synthesis flux")
  if (any(abs( sP*relSP - (.aP + parms$nuP*tvrERecycling/cpE*relSP)) > sqrEps)) stop(
    "composition P mass balance error in synthesis flux")
  fracPhiU <- (1 - parms$nu)*(decL/cnL*x$rel[["LN"]] + decR/cnR*x$rel[["RN"]]
                              + tvrERecycling/cnE*relSN +
                                tvrBOrg*(1 - cW)/cnBL*x$rel[["BN"]])
  fracPhiPU <- (1 - parms$nuP)*(decL/cpL*x$rel[["LP"]] + decR/cpR*x$rel[["RP"]]
                                + tvrERecycling/cpE*relSP +
                                  tvrBOrg*(1 - cW)/cpBL*x$rel[["BP"]])
  #
  # abbreviations for tvr, used again below that adds to R pool
  # source of tvrBOrg is B, source synE is synthesis pool, s
  tvrC <- cW*tvrBOrg*x$rel[["B"]] + (1 - parms$kNB)*synE*relSC
  tvrN <- cW*tvrBOrg/parms$cnBW*x$rel[["BN"]] +
    (1 - parms$kNB)*synE/parms$cnE*relSN
  tvrP <- cW*tvrBOrg/parms$cpBW*x$rel[["BP"]] +
    (1 - parms$kNB)*synE/parms$cpE*relSP
  .bookmarkDB <- function(){}
  dB <- uC*relUC - (respB + synE)*relSC -
    (tvrB + tvrBPred)*x$rel[["B"]]
  .dB <- (sC - respB - synE)*relSC -
    (recycB + tvrB + tvrBPred)*x$rel[["B"]]
  dBN <- (sN - synE/cnE - minN)*relSN +
    (-tvrB - tvrBPred - recycB)/cnB*x$rel[["BN"]]
  dBP <- (sP - minP - synE/cpE)*relSP +
    (-tvrB - tvrBPred - recycB)/cpB*x$rel[["BP"]]
  dBTot <- sum(dB*x$units$C)
  #if ((xOrig$frac["B"] <= 1e-16) && (dBTot < 0)) dB[] <- dBN[] <- dBP[] <- dBTot <- 0
  dL <- -decL*x$rel[["L"]]  + parms$iL*parms$relIL$C
  dLN <- -decL/cnL*x$rel[["L"]] + parms$iL/parms$cnIL*parms$relIL$N
  dLP <- -decL/cpL*x$rel[["L"]] + parms$iL/parms$cpIL*parms$relIL$P
  dR <- -decR*x$rel[["R"]]  + parms$iR*parms$relIR$C  + tvrC
  dRN <- -decR/cnR*x$rel[["RN"]]  + parms$iR/parms$cnIR*parms$relIR$N  + tvrN
  dRP <- -decR/cpR*x$rel[["RP"]]  + parms$iR/parms$cpIR*parms$relIR$P  + tvrP
  # here plant uptake as absolute parameter
  dI <-  +parms$iI*parms$relII + fracPhiU  + minN*relSN + PhiTvr*x$rel[["BN"]] -
    (parms$kIPlant + leach + immoN)*x$rel[["I"]]
  dIP <-  +parms$iIP*parms$relIIP + fracPhiPU + minP*relSP + PhiPTvr*x$rel[["BP"]] -
    (parms$kIPPlant + leachP + immoP)*x$rel[["IP"]]
  #
  if (isTRUE(parms$isFixedS)) {
    # scenario of fixed substrate
    dR[] <- dL[] <- dRN[] <- dLN[] <- dRP[] <- dLP[] <- dI[] <- dIP[] <- 0
  } else if (isTRUE(parms$isTvrNil)) {
    # scenario of enzymes and biomass not feeding back to R
    # subtract from R again an regard in mass balance check
    tvrExC <- tvrC; tvrExN <- tvrN; tvrExP <- tvrP
    dR[] <- dR - tvrC
    dRN[] <- dRN - tvrN
    dRP[] <- dRP - tvrP
  }
  #
  # make sure same order as stateNames in x (C, N P, scalars)
  resDeriv <- structure(as.numeric(
    c( dB, dR, dL, dBN, dRN, dLN, dI, dBP, dRP, dLP, dIP, dAlpha))
    , names = names(x$stateVec(x)) )
  if (any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
  #
  # checking the mass balance of fluxes
  plantNUp <- plantPUp <- 0 # checked in mass balance but is not (any more) in model
  resp <- respB + respTvr
  # biomass mass balance
  if (diff( unlist(c(uC = uC + recycB, usage = respB + synB + synE )))^2 > sqrEps )  stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(sN = uNOrg + recycB/cnB, usage = synE/parms$cnE + synB/parms$cnB + PhiB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (diff( unlist(
    c(uP = uPOrg + recycB/cpB, usage = synE/parms$cpE + synB/parms$cpB + PhiPB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance P error")
  if ((sum(dB*x$units$C) - parms$cnB*sum(dBN*x$units$N)) > sqrEps) stop(
    "biomass CN error")
  if ((sum(dB*x$units$C) - parms$cpB*sum(dBP*x$units$P)) > sqrEps) stop(
    "biomass CP error")
  # biomass turnover mass balance
  if (diff( unlist(c(tvrB + tvrBPred, usage = respTvr + tvrBOrg )))^2 > sqrEps )  stop(
    "biomass turnover mass balance C error")
  if (diff( unlist(c(
    tvrB/cnB + tvrBPred/cnB
    , usage = PhiTvr + tvrBOrg/cnB
    #, usage = PhiTvr + cW*tvrBOrg/cnBW + (1 - cW)*tvrBOrg/cnBL
  )))^2 > sqrEps )  stop(
    "biomass turnover mass balance N error")
  if (diff( unlist(c(
    (tvrB + tvrBPred)/cpB
    , usage = PhiPTvr + cW*tvrBOrg/cpBW + (1 - cW)*tvrBOrg/cpBL
     )))^2 > sqrEps )  stop("biomass turnover mass balance P error")
  if (!isTRUE(parms$isFixedS)) {
    .bookmarkDCBalance <- function(){}
    if (any(abs(
      (dB + dR + dL + tvrExC + respB*relSC + respTvr*x$rel[["B"]]) -
      (parms$iR*parms$relIR$C + parms$iL*parms$relIL$C)
      ) > 1e-3))  stop("mass balance dC error")
    if (any(abs(
      (dBN  + dRN + dLN + dI + tvrExN) -
      ((parms$iR/parms$cnIR*parms$relIR$N  + parms$iL/parms$cnIL*parms$relIL$N
       + parms$iI*parms$relII) +
      (parms$kIPlant - parms$l*x$tot["I"])*x$rel[["I"]])
    ) > sqrEps))  stop("mass balance dN error")
    if (any(abs(
      (dBP  + dRP + dLP + dIP + tvrExP) -
      ((parms$iR/parms$cpIR*parms$relIR$P  + parms$iL/parms$cpIL*parms$relIL$P
       + parms$iIP*parms$relIIP) +
      (parms$kIPPlant - parms$lP*x$tot["IP"])*x$rel[["IP"]])
    ) > sqrEps))  stop("mass balance dP error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) {
    resDeriv[names(x$frac[["R"]])] <- resDeriv[names(x$frac[["RN"]])] <-
      resDeriv[names(x$frac[["RP"]])] <- 0 }
  if (isTRUE(parms$isFixedL)) {
    resDeriv[names(x$frac[["L"]])] <- resDeriv[names(x$frac[["LN"]])] <-
      resDeriv[names(x$frac[["LP"]])] <- 0 }
  if (isTRUE(parms$isFixedI)) {
    resDeriv[names(x$frac[["I"]])] <- resDeriv[names(x$frac[["I"]])] <- 0   }
  #
  # further computations just for output for tacking the system
  ER <- alpha * parms$aE * B / parms$kN
  EL <- (1 - alpha) * parms$aE * B / parms$kN
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
  #
  if (isTRUE(parms$isRecover) ) recover()
  list( resDeriv,  c(
     resp = as.numeric(resp)
     , respO = as.numeric(respO)
     , respB = as.numeric(respB)
     , respTvr = as.numeric(respTvr)
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
    , alphaC = as.numeric(alphaC)
    , alphaN = as.numeric(alphaN), alphaP = as.numeric(alphaP)
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , cpR = as.numeric(cpR), cpL = as.numeric(cpL)
    , limER = as.numeric(limER), limEL = as.numeric(limEL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , synB = as.numeric(synB)
    , recycB = as.numeric(recycB)
    , tvrB = as.numeric(tvrB)
    , tvrBPred = as.numeric(tvrBPred)
    , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    , revRN = as.numeric(revRN)
    , revLN = as.numeric(revLN)
    , CsynB = as.numeric(CsynB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    , CsynBP = as.numeric(CsynBP)
    , uptakeC = as.numeric(uC)
    , decNLR = as.numeric(decNLR), decNE = as.numeric(decNE), decNB = as.numeric(decNB)
    #, ER = as.numeric(ER), EL = as.numeric(EL)
  , structure(as.numeric(resp*relSC)
              , names = paste("resp",names(x$units$C), sep = "_"))
  , structure(as.numeric(immoN*x$rel[["I"]])
              , names = paste("immoN",names(x$units$N), sep = "_"))
  , structure(as.numeric(minN*x$rel[["B"]])
              , names = paste("minN",names(x$units$N), sep = "_"))
  , structure(as.numeric(immoP*x$rel[["IP"]])
              , names = paste("immoP",names(x$units$P), sep = "_"))
  , structure(as.numeric(minP*x$rel[["B"]])
              , names = paste("minP",names(x$units$P), sep = "_"))
  ))
}


balanceAlphaBetweenElementLimitations <- function(
  ### compute balance between alphas of different element limitations
  alpha    ##<< numeric vector of allocation coefficients for different elements
  , CsynBE ##<< numeric vector of carbon availale for biomass synthesis
  , tauB   ##<< numeric scalar: typical microbial turnover flux for scaling
  , delta = 20  ##<< scalar smoothing factor, the higher, the steeper the transition
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


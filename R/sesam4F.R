#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam4F <- function(
  ### Soil Enzyme Steady Allocation model with detailed turnover and pool fractions
  t,xvec,parms
){
  # based on Sesam4a with tracking different fractions of pools
  # assure no negative masses
  sqrEps <- sqrt(.Machine$double.eps)
  xOrig <- xvec
  xvec <- pmax(unlist(xvec), 0.0)      # no negative masses
  # sum pool fractions
  x <- parms$multiPoolFractions$setX(parms$multiPoolFractions, xvec)
  if (any(x$tot[c("BC","BN","BP","RC","RN","RP","LC","LN","LP")] <= 0)) stop(
    "encountered zero or negative masses")
  # elemental ratios
  cnR <- x$tot["RC"]/x$tot["RN"]
  cnL <- x$tot["LC"]/x$tot["LN"]
  cpR <- x$tot["RC"]/x$tot["RP"]
  cpL <- x$tot["LC"]/x$tot["LP"]
  cW <- parms$cW
  cnE <- parms$cnE;  cnB <- parms$cnB;  cnBW <- parms$cnBW
  cpE <- parms$cpE;  cpB <- parms$cpB;  cpBW <- parms$cpBW
  if (((cW == 1) || (cW == 0)) && ((cnB != cnBW) || (cpB != cpBW))) stop(
    "If ceB != ceBW, organic turnover must be partitioned between R and DOM.",
    "cW=1 and cW=0 are not allowed.")
  cnBL <- if (cW == 1 ) cnB else (1 - cW)/(1/cnB - cW/cnBW)
  cpBL <- if (cW == 1 ) cpB else (1 - cW)/(1/cpB - cW/cpBW)
  # abbreviations of variables from x and parms
  alpha <- x$tot["alpha"]
  B <- x$tot["BC"]
  dRPot <- parms$kR * x$tot["RC"]
  dLPot <- parms$kL * x$tot["LC"]
  immoNPot <- parms$iBN * x$tot["IN"]
  immoPPot <- parms$iBP * x$tot["IP"]
  aeB <- parms$aE*B
  kmN <- parms$kmN #parms$km*parms$kN
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
  PhiNTvr <- respTvr/parms$cnB
  PhiPTvr <- respTvr/parms$cpB
  #
  # organic uptake
  tvrERecycling <- parms$kNB*synE # part of enzyme turnover feeds back to uptake
  decNLR <- decL/cnL + decR/cnR
  decNE <- tvrERecycling/cnE
  decNB <- tvrBOrg*(1 - cW)/cnBL
  uC <- decL + decR + tvrERecycling + tvrBOrg*(1 - cW)
  # relURecyc depends on relUC because composition of enzymes is that of uptake
  #fracUC <- decL*x$rel[["LC"]] + decR*x$rel[["RC"]] +
  #  tvrERecycling*relSC + tvrBOrg*(1 - cW)*x$rel[["BC"]]
  fracUNOrg <- parms$nuN*(decL/cnL*x$rel[["LN"]] + decR/cnR*x$rel[["RN"]] +
                           (tvrERecycling/cnE + tvrBOrg*(1 - cW)/cnBL)*x$rel[["BN"]])
  uNOrg <- parms$nuN*(decL/cnL + decR/cnR + tvrERecycling/cnE + tvrBOrg*(1 - cW)/cnBL)
  #uNOrg1 <- parms$nuN*(decNLR + decNE + decNB)
  #if (abs(uNOrg - uNOrg1) > 1e-14) stop("fracUNOrg error")
  uPOrg <- parms$nuP*(decL/cpL + decR/cpR +
                        tvrERecycling/cpE + tvrBOrg*(1 - cW)/cpBL)
  #uPOrg1 <- parms$nuP*(decL/cpL + decR/cpR + tvrERecycling/cpE + tvrBOrg*(1 - cW)/cpBL)
  #if (abs(uPOrg - uPOrg1) > 1e-14) stop("fracUPOrg error")
  #
  # elemental limitations by potential biomass synthesis and respiration
  CsynBCt <- uC - rM - synE/parms$eps
  CsynBC <- if(CsynBCt > 0) parms$eps*CsynBCt else CsynBCt
  CsynBN <- cnB*(uNOrg + immoNPot - synE/cnE)
  CsynBP <- cpB*(uPOrg + immoPPot - synE/cpE)
  synB <- min(CsynBC, CsynBN, CsynBP)
  starvB <- if (synB > 0) 0 else -synB
  rG <- if (synB > 0) (1 - parms$eps)/parms$eps*synB else 0
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
  dAlpha <- (alphaTarget - alpha) * (synB + tvrB + tvrBPred)/B
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC - (synE/parms$eps + synB + rG + rM)
  PhiNB <- uNOrg - synB/cnB - synE/cnE
  PhiPB <- uPOrg - synB/cpB - synE/cpE
  #
  tvrN <-  +cW*tvrBOrg/parms$cnBW   + (1 - parms$kNB)*synE/parms$cnE
  tvrP <-  +cW*tvrBOrg/parms$cpBW   + (1 - parms$kNB)*synE/parms$cpE
  # fluxes leaving the system (will be set in scen where trv does not feed back)
  tvrExC <- x$units$C;  tvrExN <- x$units$N;  tvrExP <- x$units$P
  tvrExC[] <- tvrExN[] <- tvrExP[] <- 0
  #
  leachN <- parms$lN*x$tot["IN"]
  PhiNU <- (1 - parms$nuN)*(decL/cnL + decR/cnR + tvrERecycling/cnE + tvrBOrg*(1 - cW)/cnBL)
  leachP <- parms$lP*x$tot["IP"]
  PhiPU <- (1 - parms$nuP)*(decL/cpL + decR/cpR + tvrERecycling/cpE + tvrBOrg*(1 - cW)/cpBL)
  #
  respB <- (synE)/parms$eps*(1 - parms$eps)  + rG + rM + respO
  resp <- respB + respTvr
  immoN <- max(0,-PhiNB); minN <- max(0,PhiNB)
  immoP <- max(0,-PhiPB); minP <- max(0,PhiPB)
  #
  # after sC and immo is known, relSC can be computed (see SteadyEnzyme4a.Rmd)
  # fluxes that feed synthesis, s, and mineralization fluxes
  sC <- uC + starvB
  sN <- uNOrg + starvB/cnB + immoN
  sP <- uPOrg + starvB/cpB + immoP
  cnS <- sC/sN
  cpS <- sC/sP
  .aC <- decL*x$rel[["LC"]] + decR*x$rel[["RC"]] + tvrBOrg*(1 - cW)*x$rel[["BC"]] +
    starvB*x$rel[["BC"]]
  .aN <- parms$nuN * (decL/cnL*x$rel[["LN"]] + decR/cnR*x$rel[["RN"]] +
    tvrBOrg*(1 - cW)/cnBL*x$rel[["BN"]]) +
    starvB/cnB*x$rel[["BN"]] + immoN*x$rel[["IN"]]
  .aP <- parms$nuP * (decL/cpL*x$rel[["LP"]] + decR/cpR*x$rel[["RP"]] +
    tvrBOrg*(1 - cW)/cpBL*x$rel[["BP"]]) +
    starvB/cpB*x$rel[["BP"]] + immoP*x$rel[["IP"]]
  relSC <- .aC/(sC - tvrERecycling)
  relSN <- .aN/(sN - parms$nuN*tvrERecycling/cnE)
  relSP <- .aP/(sP - parms$nuP*tvrERecycling/cpE)
  relUC <- (sC*relSC - starvB*x$rel[["BC"]])/uC
  if (any(abs( sC*relSC - (.aC + tvrERecycling*relSC)) > sqrEps)) stop(
    "composition C mass balance error in synthesis flux")
  if (any(abs( sN*relSN - (.aN + parms$nuN*tvrERecycling/cnE*relSN)) > sqrEps)) stop(
    "composition N mass balance error in synthesis flux")
  if (any(abs( sP*relSP - (.aP + parms$nuP*tvrERecycling/cpE*relSP)) > sqrEps)) stop(
    "composition P mass balance error in synthesis flux")
  fracPhiU <- (1 - parms$nuN)*(decL/cnL*x$rel[["LN"]] + decR/cnR*x$rel[["RN"]]
                              + tvrERecycling/cnE*relSN +
                                tvrBOrg*(1 - cW)/cnBL*x$rel[["BN"]])
  fracPhiPU <- (1 - parms$nuP)*(decL/cpL*x$rel[["LP"]] + decR/cpR*x$rel[["RP"]]
                                + tvrERecycling/cpE*relSP +
                                  tvrBOrg*(1 - cW)/cpBL*x$rel[["BP"]])
  #
  # abbreviations for tvr, used again below that adds to R pool
  # source of tvrBOrg is B, source synE is synthesis pool, s
  tvrC <- cW*tvrBOrg*x$rel[["BC"]] + (1 - parms$kNB)*synE*relSC
  tvrN <- cW*tvrBOrg/parms$cnBW*x$rel[["BN"]] +
    (1 - parms$kNB)*synE/parms$cnE*relSN
  tvrP <- cW*tvrBOrg/parms$cpBW*x$rel[["BP"]] +
    (1 - parms$kNB)*synE/parms$cpE*relSP
  .bookmarkDB <- function(){}
  # dBC <- uC*relUC - (respB + synE)*relSC -
  #   (tvrB + tvrBPred)*x$rel[["BC"]]
  dBC <- (sC - respB - synE)*relSC +
    (-tvrB - tvrBPred - starvB)*x$rel[["BC"]]
  dBN <- (sN - minN - synE/cnE)*relSN +
    (-tvrB - tvrBPred - starvB)/cnB*x$rel[["BN"]]
  dBP <- (sP - minP - synE/cpE)*relSP +
    (-tvrB - tvrBPred - starvB)/cpB*x$rel[["BP"]]
  dBTot <- sum(dBC*x$units$C)
  #if ((xOrig$frac["BC"] <= 1e-16) && (dBTot < 0)) dBC[] <- dBN[] <- dBP[] <- dBTot <- 0
  dLC <- -decL*x$rel[["LC"]]  + parms$iL*parms$relIL$C
  dLN <- -decL/cnL*x$rel[["LC"]] + parms$iL/parms$cnIL*parms$relIL$N
  dLP <- -decL/cpL*x$rel[["LC"]] + parms$iL/parms$cpIL*parms$relIL$P
  dRC <- -decR*x$rel[["RC"]]  + parms$iR*parms$relIR$C  + tvrC
  dRN <- -decR/cnR*x$rel[["RN"]]  + parms$iR/parms$cnIR*parms$relIR$N  + tvrN
  dRP <- -decR/cpR*x$rel[["RP"]]  + parms$iR/parms$cpIR*parms$relIR$P  + tvrP
  # here plant uptake as absolute parameter
  dITot <-  +parms$iIN + sum(fracPhiU*x$units$N) + minN + PhiNTvr -
    (parms$kINPlant*x$tot["IN"] + leachN + immoN)
  dIN <-  +parms$iIN*parms$relII + fracPhiU  + minN*relSN + PhiNTvr*x$rel[["BN"]] -
    (parms$kINPlant*x$tot["IN"] + leachN + immoN)*x$rel[["IN"]]
  dIP <-  +parms$iIP*parms$relIIP + fracPhiPU + minP*relSP + PhiPTvr*x$rel[["BP"]] -
    (parms$kIPPlant*x$tot["IP"] + leachP + immoP)*x$rel[["IP"]]
  dResp <- respB*relSC + respTvr*x$rel[["BC"]]
  dLeachN <- leachN*x$rel[["IN"]]
  dLeachP <- leachP*x$rel[["IP"]]
  .tmp.f.displaySumDeriv <- function(){
    c( BC = sum(dBC*x$units$C), RC = sum(dRC*x$units$C), LC = sum(dLC*x$units$C)
       , RN = sum(dRN*x$units$N), LN = sum(dLN*x$units$N), IN = sum(dIN*x$units$N)
       , RP = sum(dRP*x$units$P), LP = sum(dLP*x$units$P), IP = sum(dIP*x$units$P)
       , resp = sum(dResp*x$units$C
       , leachN = sum(dLeachN*x$units$N), leachP = sum(dLeachP*x$units$P))
       , alpha = as.numeric(dAlpha)
  )}
  #
  if (isTRUE(parms$isFixedS)) {
    # scenario of fixed substrate
    dRC[] <- dLC[] <- dRN[] <- dLN[] <- dRP[] <- dLP[] <- dIN[] <- dIP[] <- 0
  } else if (isTRUE(parms$isTvrNil)) {
    # scenario of enzymes and biomass not feeding back to R
    # subtract from R again an regard in mass balance check
    tvrExC <- tvrC; tvrExN <- tvrN; tvrExP <- tvrP
    dRC[] <- dRC - tvrC
    dRN[] <- dRN - tvrN
    dRP[] <- dRP - tvrP
  }
  #
  # make sure same order as stateNames in x (C, N P, scalars)
  resDeriv <- structure(as.numeric(
    c( dBC, dRC, dLC, dResp
       , dBN, dRN, dLN, dIN, dLeachN
       , dBP, dRP, dLP, dIP, dLeachP
       , dAlpha))
    , names = names(x$stateVec(x)) )
  if (any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
  #
  # checking the mass balance of fluxes
  # biomass mass balance
  if (diff( unlist(c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps )  stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(sN = uNOrg, usage = synE/parms$cnE + synB/parms$cnB + PhiNB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (diff( unlist(
    c(uP = uPOrg, usage = synE/parms$cpE + synB/parms$cpB + PhiPB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance P error")
  if ((sum(dBC*x$units$C) - parms$cnB*sum(dBN*x$units$N)) > sqrEps) stop(
    "biomass CN error")
  if ((sum(dBC*x$units$C) - parms$cpB*sum(dBP*x$units$P)) > sqrEps) stop(
    "biomass CP error")
  # biomass turnover mass balance
  if (diff( unlist(c(tvrB + tvrBPred, usage = respTvr + tvrBOrg )))^2 > sqrEps )  stop(
    "biomass turnover mass balance C error")
  if (diff( unlist(c(
    tvrB/cnB + tvrBPred/cnB
    , usage = PhiNTvr + tvrBOrg/cnB
    #, usage = PhiNTvr + cW*tvrBOrg/cnBW + (1 - cW)*tvrBOrg/cnBL
  )))^2 > sqrEps )  stop(
    "biomass turnover mass balance N error")
  if (diff( unlist(c(
    tvrB/cpB + tvrBPred/cpB
    , usage = PhiPTvr + tvrBOrg/cpB
    #, usage = PhiNTvr + cW*tvrBOrg/cnBW + (1 - cW)*tvrBOrg/cnBL
  )))^2 > sqrEps )  stop(
    "biomass turnover mass balance P error")
  if (!isTRUE(parms$isFixedS)) {
    .bookmarkDCBalance <- function(){}
    if (any(abs(
      (dBC + dRC + dLC + tvrExC + dResp) -
      (parms$iR*parms$relIR$C + parms$iL*parms$relIL$C)
      ) > sqrEps))  stop("mass balance dC error")
    if (any(abs(
      (dBN  + dRN + dLN + dIN + tvrExN) -
      (parms$iR/parms$cnIR*parms$relIR$N  + parms$iL/parms$cnIL*parms$relIL$N +
       parms$iIN*parms$relII - parms$kINPlant*x$tot["IN"]*x$rel[["IN"]] - dLeachN)
    ) > sqrEps))  stop("mass balance dN error")
    if (any(abs(
      (dBP  + dRP + dLP + dIP + tvrExP) -
      (parms$iR/parms$cpIR*parms$relIR$P  + parms$iL/parms$cpIL*parms$relIL$P +
       parms$iIP*parms$relIIP - parms$kIPPlant*x$tot["IP"]*x$rel[["IP"]] - dLeachP)
    ) > sqrEps))  stop("mass balance dP error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) {
    resDeriv[names(x$frac[["RC"]])] <- resDeriv[names(x$frac[["RN"]])] <-
      resDeriv[names(x$frac[["RP"]])] <- 0 }
  if (isTRUE(parms$isFixedL)) {
    resDeriv[names(x$frac[["LC"]])] <- resDeriv[names(x$frac[["LN"]])] <-
      resDeriv[names(x$frac[["LP"]])] <- 0 }
  if (isTRUE(parms$isFixedI)) { resDeriv[names(x$frac[["IN"]])] <- 0 }
  if (isTRUE(parms$isFixedIP)) { resDeriv[names(x$frac[["IP"]])] <- 0 }
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
  PhiNBU <- PhiNB + PhiNU
  # total mineralization flux including microbial turnover
  PhiNTotal <- PhiNBU + PhiNTvr
  #
  if (isTRUE(parms$isRecover) ) recover()
  list( resDeriv,  c(
     respTotal = as.numeric(resp)
     , respO = as.numeric(respO)
     , respB = as.numeric(respB)
     , respTvr = as.numeric(respTvr)
    , PhiNTotal = as.numeric(PhiNTotal)
    , PhiNB = as.numeric(PhiNB), PhiNU = as.numeric(PhiNU)
    , PhiNTvr = as.numeric(PhiNTvr)
    , PhiNBU = as.numeric(PhiNBU)
    , immoNPot = as.numeric(immoNPot)
    , PhiPTotal = as.numeric(PhiPB + PhiPU + PhiNTvr)
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
    , starvB = as.numeric(starvB)
    , tvrB = as.numeric(tvrB)
    , tvrBPred = as.numeric(tvrBPred)
    , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    , revRN = as.numeric(revRN)
    , revLN = as.numeric(revLN)
    , synB = as.numeric(synB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    , CsynBP = as.numeric(CsynBP)
    , uptakeC = as.numeric(uC)
    , synC = as.numeric(sC)
    , synN = as.numeric(sN)
    , decNL = as.numeric(decL/cnL), decNR = as.numeric(decR/cnR)
    , decNE = as.numeric(decNE), decNB = as.numeric(decNB)
    #, ER = as.numeric(ER), EL = as.numeric(EL)
  , structure(as.numeric(resp*relSC)
              , names = paste("respTotal",names(x$units$C), sep = "_"))
  , structure(as.numeric(immoN*x$rel[["IN"]])
              , names = paste("immoN",names(x$units$N), sep = "_"))
  , structure(as.numeric(minN*x$rel[["BN"]])
              , names = paste("minN",names(x$units$N), sep = "_"))
  , structure(as.numeric(immoP*x$rel[["IP"]])
              , names = paste("immoP",names(x$units$P), sep = "_"))
  , structure(as.numeric(minP*x$rel[["BP"]])
              , names = paste("minP",names(x$units$P), sep = "_"))
  ))
}


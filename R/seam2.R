#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSeam2 <- function(
  ### Soil Enzyme Allocation Model with explicit immobilization
  t, x, parms
){
  x <- pmax(unlist(x), 1e-16)      # no negative masses
  #
  ER <- x["ER"]
  EL <- x["EL"]
  cnR <- x["R"]/x["RN"]
  cnL <- x["L"]/x["LN"]
  cnE <- parms$cnE #alphaC*parms$cnER + (1-alphaC)*parms$cnEL
  cnB <- parms$cnB
  #
  rM <- parms$m * x["B"]          # maintenance respiration
  tvrB <- parms$tau*x["B"]        # microbial turnover
  tvrER <- parms$kNR*ER
  tvrEL <- parms$kNL*EL
  #
  decRp <- parms$kR * x["R"]
  decLp <- parms$kL * x["L"]
  limER <- ER / (parms$kmR + ER)
  limEL <- EL / (parms$kmL + EL)
  decR <- decRp * limER
  decL <- decLp * limEL
  #
  tvrERecycling <- parms$kNB*(tvrER + tvrEL)
  uC <- decC <- decR + decL + tvrERecycling
  #
  synE <- parms$aE * x["B"]     # total enzyme production per microbial biomass
  #synE <- parms$aEu * uC       # total enzyme production per uptake
  # C required for biomass and growth respiration under C limitation
  CsynBC <- uC - synE/parms$eps - rM
  #
  # Nitrogen balance
  decN <- decR/cnR + decL/cnL + tvrERecycling/parms$cnE
  # plants get at maximum half of the decomposed organic N
  plantNUp <- pmin(parms$plantNUp, decN/2)
  # mineralization due to soil heterogeneity (Manzoni 08)
  PhiU <- (1 - parms$nu) * (decN - plantNUp)
  # immobilization flux
  immoPot <- parms$iB * x["I"]
  uNSubstrate <- (decN - plantNUp - PhiU)    # plant uptake also of organic N
  uNPot <-  immoPot + uNSubstrate
  NsynBN <- uNPot - synE/cnE
  # C required for biomass growth and growth respiration under N limitation
  CsynBN <- (NsynBN*cnB)/parms$eps
  #
  # compute C available for biomass, this time without accounting immobilization flux
  NsynBNSubstrate <- uNSubstrate - synE/cnE
  CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
  #
  # N limited taking immobilization into account
  isLimN <- CsynBN <  CsynBC
  # N limited based on substrate uptake (without accounting for immobilization)
  isLimNSubstrate <-  CsynBNSubstrate < CsynBC
  CsynB = if (isLimN ) CsynBN else CsynBC
  #
  # average cN reqruired according to enzyme allocation accoording to C efficiencies
  #synE <- uC*parms$aE
  #uCGrowth <- uC - synE/parms$eps - rM    # C for both growth and growth respiration
  if (CsynB > 0) {
      synB <- parms$eps*CsynB
      rG <- (1 - parms$eps)*CsynB
  }else{
      synB <- CsynB    # negative
      rG <- 0
  }
  # imbalance fluxes
  respSynE <- (1 - parms$eps)/parms$eps * synE
  respO <- uC - (synE + respSynE + synB + rG + rM)
  # imbalance mineralization/immobilization flux substrates only
  PhiB <- uNSubstrate - (synE/parms$cnE + synB/parms$cnB)
  # imbalance taking into account potential immobilization
  MmImb <- uNPot - (synE/parms$cnE + synB/parms$cnB)
  #PhiB2 <- MmImb - immoPot # should equal PhiB
  PhiBU <- PhiB + PhiU
  respTvr <- (1 - parms$epsTvr) * tvrB
  PhiTvr <- respTvr/parms$cnB
  # net microbial mineralization/imm when taking into account uptake mineralization
  PhiBU <- PhiB + PhiU
  # total mineralization flux including microbial turnover
  PhiTotal <- PhiBU + PhiTvr
  #
  # Revenue strategy
  revLC <- decLp / (parms$kNL) / (parms$kmL + EL)
  revRC <- decRp / (parms$kNR) / (parms$kmR + ER)
  revLN <- revLC * cnE/cnL
  revRN <- revRC * cnE/cnR
  alphaC <- revRC / (revLC + revRC)
  alphaN <- revRN / (revLN + revRN)
  #
  if (isTRUE(parms$isAlphaMatch) || isTRUE(parms$isAlphaEnzMax)) {
    synB0 <- max(0, synB)
    #cnOpt <- (parms$cnE*synE + parms$cnB*synB0)/(synE+synB0)  # optimal biomass ratio
    cnOpt <- (synE + synB0)/(synE/cnE + synB0/parms$cnB)  # optimal biomass ratio
    # calcMatchAlphaEnz
    alpha <- alphaMatch <- calcMatchAlphaSeam2(
      E = ER + EL, decPotR = decRp, decPotL = decLp, rMaint = rM
      ,cnOpt = cnOpt
      , cnR = cnR, cnL = cnL
      , parms = parms
      , imm = max(0, immoPot - PhiU)  # match strategy for N-gain by immobilization
      #, imm = 0               # match strategy not relying on potential immobilization?
    )
    if (isTRUE(parms$isAlphaEnzMax)) {
      # 0 for PhiB < 0.05, i.e. immobilization, i.e. C limitation
      w <- getGradualWeight(PhiB,0.005)
      alpha <- (1 - w)*0.5 + w*alphaMatch
    }
  } else {
    alpha <- balanceAlphaBetweenCNLimitations(
      alphaC, alphaN, CsynBN, CsynBC
      , NsynBC = parms$eps*CsynBC/cnB, NsynBN)
  }
  #c(alphaC, alphaN, alpha)
  if (isTRUE(parms$isAlphaFix) )  alpha <- 0.5
  #
  # tvr feeding back to R pool, assume that N in SOM for resp (by epsTvr) is mineralized
  tvrC <- +parms$epsTvr*tvrB + (1 - parms$kNB)*(tvrER + tvrEL)
  tvrN <- +parms$epsTvr*tvrB/parms$cnB + (1 - parms$kNB)*(tvrER + tvrEL)/parms$cnE
  tvrExC <- 0     # fluxes leaving the system
  tvrExN <- 0
  #
  leach <- parms$l*x["I"]
  #
  dB <- synB - tvrB
  dER <- +alpha*synE  - tvrER
  dEL <- +(1 - alpha)*synE  - tvrEL
  dL <- -decL + parms$iL
  dLN <- -decL/cnL  + parms$iL/parms$cnIL
  dR <- -decR + parms$iR + tvrC
  dRN <- -decR/cnR + parms$iR/parms$cnIR + tvrN
  #dI <- +parms$iI +MmB +PhiTvr -(parms$kIPlant+parms$l)*x["I"]
  # plant uptake as absolute parameter
  dI <- +parms$iI - parms$kIPlant - leach + PhiB + PhiU + PhiTvr
  #if (dI > 0.01 ) recover()
  #
  if (isTRUE(parms$isFixedS)) {
    # scenario of fixed substrate
    dR <- dL <- dRN <- dLN <- dI <- 0
  } else if (isTRUE(parms$isTvrNil)) {
    # scenario of enzymes and biomass not feeding back to R
    dR <- +parms$iR - decR
    dRN <- +parms$iR/parms$cnIR - decR/cnR
    tvrExC <- tvrC
    tvrExN <- tvrN
  }
  respB <- respSynE + rG + rM + respO
  resp <- respB + respTvr
  #
  resDeriv <- structure(as.numeric(
    c( dB, dER, dEL, dR, dRN, dL, dLN, dI))
    , names = c("dB","dER","dEL","dR","dRN","dL","dLN","dI"))
  if (any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
  sqrEps <- sqrt(.Machine$double.eps)
  # parms$iL - (decL + dL)
  # parms$iR + +tvrC -(decR + dR)
  #
  if (diff( unlist(c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps ) stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(uN = uNSubstrate, usage = synE/parms$cnE + synB/parms$cnB + PhiB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dER + dEL + dR + dL + tvrExC + resp, parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + (dER + dEL)/parms$cnE  + dRN + dLN + dI + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL - plantNUp  + parms$iI -
         parms$kIPlant - parms$l*x["I"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
  }
  # keeping R,L, or I constant
  if (isTRUE(parms$isFixedR)) { resDeriv["dR"] <- resDeriv["dRN"] <-  0   }
  if (isTRUE(parms$isFixedL)) { resDeriv["dL"] <- resDeriv["dLN"] <-  0   }
  if (isTRUE(parms$isFixedI)) { resDeriv["dI"] <-  0   }
  #
  if (isTRUE(parms$isRecover) ) recover()
  list( resDeriv, c(
    respO = as.numeric(respO)
    #, MmB = as.numeric(MmB)
    , PhiB = as.numeric(PhiB), PhiU = as.numeric(PhiU)
    , PhiTvr = as.numeric(PhiTvr)
    , PhiBU = as.numeric(PhiBU), PhiTotal = as.numeric(PhiTotal)
    , immoPot = as.numeric(immoPot), MmImb = as.numeric(MmImb)
    , alpha = as.numeric(alpha)
    , alphaC = as.numeric(alphaC), alphaN = as.numeric(alphaN)
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , limER = as.numeric(limER), limEL = as.numeric(limEL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , resp = as.numeric(resp)
    , respB = as.numeric(respB), respTvr = as.numeric(respTvr)
    , tvrB = as.numeric(tvrB)
    , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    , revRN = as.numeric(revRN), revLN = as.numeric(revLN)
    , pCsyn = as.numeric(CsynBC / CsynBN)
    , CsynReq = as.numeric(CsynBN), Csyn = as.numeric(CsynBC)
    , pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    , NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    , dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    , dI = as.numeric(dI)
    #    wCLim = (CsynBN/CsynBC)^delta
    #    wNLim = (parms$eps*CsynBC/cnB / NsynBN)^delta
    , uC = as.numeric(uC), synB = as.numeric(synB)
    , decC = as.numeric(decC), decN = as.numeric(decN)
    , plantNUp = plantNUp
  ))
}

balanceAlphaBetweenCNLimitations <- function(
  ### compute balance between alphaC and alphaN based on C and N based biomass synthesis
  alphaC, alphaN, CsynBN, CsynBC, NsynBC, NsynBN, delta = 200
  , ...   ##<< allow passing further parameters on overloading
){
  # if there is only small potential of immobilizalization, do a smooth
  # transition between alphaC and alphaN
  # # min to avoid  + Inf
  wCLim = min( .Machine$double.xmax, (CsynBN/CsynBC)^delta )
  wNLim = min( .Machine$double.xmax, (NsynBC/NsynBN)^delta )
  alpha <- (wCLim*alphaC + wNLim*alphaN) / (wCLim + wNLim)
  alpha
}

.tmp.f <- function(
  ###  balance alphaC and alphaN based on current and potential immobilization fluxes
  alphaC, alphaN, isLimN, isLimNSubstrate, immoAct, immoPot
){
  # .depr.balanceAlphaLargeImmobilization
  ##details<<
  ## deprecated because its superceded by balanceAlphaBetweenCNLimitations
  ## where biosynthesis already takes into account the
  ## potential immobilization
  if (isLimN ) return(alphaN)
  if (isLimNSubstrate) {
    # overall C limited, but only with acquiring N by immobilization
    # increase proportion into N aquiring enzymes as the proportion of immobilization
    # to its potential increases
    pN <- immoAct / immoPot
    # increase proportion into N aquiring enzymes as the proportion of biomass synthesis
    # that is possible due to immobilization
    # pN <- (CSynB - CsynBNSubstrate)/CsynBNSubstrate
    alpha <- alphaC + pN*(alphaN - alphaC)
    return(alpha)
  }
  alphaC
}

plotResSeam2 <- function(
  ### plotting results of Seam2
  res, legendPos = "topleft"
  , cls = c("B","ER","EL","respO","PhiTotal")
  , xlab = "time (yr)"
  , ylab = "gC/m2 , gN/m2"
  , subsetCondition = TRUE
  , ...
){
  #res$B100 <- res$B/100
  res$ER_10 <- res$ER*10
  res$EL_10 <- res$EL*10
  res$ER_1000 <- res$ER*1000
  res$EL_1000 <- res$EL*1000
  res$alpha100 <- res$alpha*100
  res$B100 <- res$B * 100
  res$B10 <- res$B * 10
  res$Rr <- res$R / res$R[1]  * 100 #max(res$R, na.rm = TRUE)
  res$Lr <- res$L / res$L[1]  * 100 #max(res$L, na.rm = TRUE)
  res$I100 <- res$I  * 100
  res$mPhiB10 <- -res$PhiB  * 10
  res$mPhiBU10 <- -res$PhiBU  * 10
  res$mPhiTotal10 <- -res$PhiTotal  * 10
  #res$mMmB10 <- -res$MmB  * 10
  #res$MmB100 <- res$MmB  * 100
  #res$MmB1000 <- res$MmB  * 1000
  res$cnR10 <- res$cnR  * 10
  res$dR <- c(NA, diff(res$R)/diff(res$time))
  res$dS <- c(NA, diff(res$R + res$L)/diff(res$time))
  #cls <- c("B100","ER_10","EL_10","respO","PhiTotal")
  #cls <- c("B100","ER_10","EL_10","respO","PhiTotal","eff")
  #cls <- c("B1000","ER_10","EL_10","PhiTotal")
  #cls <- c("B","E","respO","PhiTotal","R","L")
  matplot( res[subsetCondition,1], res[subsetCondition,cls], type = "l", lty = 1:20
           , col = 1:20, xlab = xlab, ylab = "", ...)
  mtext( ylab, 2, line = 2.5, las = 0)
  legend(legendPos, inset = c(0.01,0.01), legend = cls, lty = 1:20, col = 1:20)
}


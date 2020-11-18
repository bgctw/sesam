#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSeam3b <- function(
  ### Soil Enzyme Allocation Model with state variable alpha and element-weighted cost
  t, x, parms
  , fBalanceAlpha = balanceAlphaBetweenCNLimitationsExp
){
  x <- pmax(unlist(x), 1e-16)      # no negative masses
  #
  ER <- x["ER"]
  EL <- x["EL"]
  cnR <- x["R"]/x["RN"]
  cnL <- x["L"]/x["LN"]
  cnE <- parms$cnE #alphaC*parms$cnER + (1-alphaC)*parms$cnEL
  cnB <- parms$cnB
  B <- x["B"]
  #
  rM <- parms$m * B          # maintenance respiration
  tvrB <- parms$tau*B        # microbial turnover
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
  synE <- parms$aE * B     # total enzyme production per microbial biomass
  #synE <- parms$aEu * uC       # total enzyme production per uptake
  # C required for biomass and growth respiration under C limitation
  CsynBCt <- uC - synE/parms$eps - rM
  CsynBC <- if(CsynBCt > 0) parms$eps*CsynBCt else CsynBCt
  #
  # Nitrogen balance
  decN <- decR/cnR + decL/cnL + tvrERecycling/parms$cnE
  # plants get at maximum half of the decomposed organic N
  plantNUpOrg <- if (!is.null(parms$plantNUpOrg))
    pmin(parms$plantNUpOrg, decN/2) else 0
  # mineralization due to soil heterogeneity (Manzoni 08)
  PhiNU <- (1 - parms$nuN) * (decN - plantNUpOrg)
  # immobilization flux
  immoNPot <- parms$iBN * x["IN"]
  uNSubstrate <- (decN - plantNUpOrg - PhiNU) # plant uptake also of organic N
  uNPot <-  immoNPot + uNSubstrate
  NsynBN <- uNPot - synE/cnE
  # C required for biomass growth under N limitation
  CsynBN <- NsynBN*cnB
  #
  # compute C available for biomass, now without accounting immobilization flux
  NsynBNSubstrate <- uNSubstrate - synE/cnE
  CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
  #
  # N limited taking immobilization into account
  isLimN <- CsynBN <  CsynBC
  # N limited based on substrate uptake (without accounting for immobilization)
  isLimNSubstrate <-  CsynBNSubstrate < CsynBC
  synB = if (isLimN ) CsynBN else CsynBC
  #
  rG <- if (synB > 0) (1 - parms$eps)/parms$eps*synB else 0
  # imbalance fluxes
  respSynE <- (1 - parms$eps)/parms$eps * synE
  respO <- uC - (synE + respSynE + synB + rG + rM)
  # imbalance mineralization/immobilization flux substrates only
  PhiNB <- uNSubstrate - (synE/parms$cnE + synB/parms$cnB)
  # imbalance taking into account potential immobilization
  MmImb <- uNPot - (synE/parms$cnE + synB/parms$cnB)
  #PhiB2 <- MmImb - immoNPot # should equal PhiNB
  PhiNBU <- PhiNB + PhiNU
  respTvr <- (1 - parms$epsTvr) * tvrB
  PhiNTvr <- respTvr/parms$cnB
  # net microbial min/imm when taking into account uptake mineralization
  PhiNBU <- PhiNB + PhiNU
  # total mineralization flux including microbial turnover
  PhiNTotal <- PhiNBU + PhiNTvr
  #
  # Revenue strategy
  revLC <- decLp / (parms$kNL) / (parms$kmL + EL)
  revRC <- decRp / (parms$kNR) / (parms$kmR + ER)
  revLN <- revLC * cnE/cnL
  revRN <- revRC * cnE/cnR
  alphaC <- revRC / (revLC + revRC)
  alphaN <- revRN / (revLN + revRN)
  #
  alphaTarget <- fBalanceAlpha(
    alphaC, alphaN, CsynBN, CsynBC
    , tauB = parms$tau*B
    , NsynBC = parms$eps*CsynBC/cnB, NsynBN)
  alpha <- x["alpha"]
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) *  (parms$tau + abs(synB)/B)
  #
  # tvr to R pool, assume that N in SOM for resp (by epsTvr) is mineralized
  tvrC <- +parms$epsTvr*tvrB + (1 - parms$kNB)*(tvrER + tvrEL)
  tvrN <- +parms$epsTvr*tvrB/parms$cnB +
    (1 - parms$kNB)*(tvrER + tvrEL)/parms$cnE
  tvrExC <- 0     # fluxes leaving the system
  tvrExN <- 0
  #
  leachN <- parms$lN*x["IN"]
  plantNUpPot <- parms$kINPlant*x["IN"]
  plantNUp <- if (!is.null(parms$plantNUpAbs)) {
    min(parms$plantNUpAbs, plantNUpPot)
  } else {
    plantNUpPot
  }
  #
  dB <- synB - tvrB
  dER <- +alpha*synE  - tvrER
  dEL <- +(1 - alpha)*synE  - tvrEL
  dL <- -decL + parms$iL
  dLN <- -decL/cnL  + parms$iL/parms$cnIL
  dR <- -decR + parms$iR + tvrC
  dRN <- -decR/cnR + parms$iR/parms$cnIR + tvrN
  #dIN <- +parms$iIN +MmB +PhiNTvr -(parms$kINPlant+parms$lN)*x["IN"]
  # plant uptake as absolute parameter
  dIN <- +parms$iIN - plantNUp - leachN + PhiNB + PhiNU + PhiNTvr
  #if (dIN > 0.01 ) recover()
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
  respB <- respSynE + rG + rM + respO
  resp <- respB + respTvr
  #
  resDeriv <- structure(as.numeric(
    c( dB, dER, dEL, dR, dRN, dL, dLN, dIN, dAlpha))
    , names = c("dB","dER","dEL","dR","dRN","dL","dLN","dIN","dAlpha"))
  if (any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
  sqrEps <- sqrt(.Machine$double.eps)
  # parms$iL - (decL + dL)
  # parms$iR + +tvrC -(decR + dR)
  #
  if (diff( unlist(c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps ) stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(uN = uNSubstrate, usage = synE/parms$cnE + synB/parms$cnB + PhiNB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dER + dEL + dR + dL + tvrExC + resp, parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + (dER + dEL)/parms$cnE  + dRN + dLN + dIN + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL + parms$iIN -
         plantNUp - plantNUpOrg - parms$lN*x["IN"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
  }
  # keeping R,L, or IN constant
  if (isTRUE(parms$isFixedR)) { resDeriv["dR"] <- resDeriv["dRN"] <-  0   }
  if (isTRUE(parms$isFixedL)) { resDeriv["dL"] <- resDeriv["dLN"] <-  0   }
  if (isTRUE(parms$isFixedI)) { resDeriv["dIN"] <-  0   }
  #
  if (isTRUE(parms$isRecover) ) recover()
  list( resDeriv, c(
    respO = as.numeric(respO)
    #, MmB = as.numeric(MmB)
    , PhiNB = as.numeric(PhiNB), PhiNU = as.numeric(PhiNU)
    , PhiNTvr = as.numeric(PhiNTvr)
    , PhiNBU = as.numeric(PhiNBU), PhiNTotal = as.numeric(PhiNTotal)
    , immoNPot = as.numeric(immoNPot), MmImb = as.numeric(MmImb)
    , alphaTarget = as.numeric(alphaTarget)
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
    , pNsyn = as.numeric(NsynBN / (CsynBC/cnB) )
    , NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    , dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    , dIN = as.numeric(dIN)
    #    wCLim = (CsynBN/CsynBC)^delta
    #    wNLim = (parms$eps*CsynBC/cnB / NsynBN)^delta
    , uC = as.numeric(uC), synB = as.numeric(synB)
    , decC = as.numeric(decC), decN = as.numeric(decN)
    , plantNUp = plantNUp
  ))
}

balanceAlphaBetweenCNLimitationsExp <- function(
  ### compute balance between alphaC and alphaN based on E_i biomass synthesis
  alphaC, alphaN, CsynBN, CsynBC, ..., delta = 20, tauB
){
  ##details<< if there is only small potential of immobilizalization,
  ## do a smooth
  ## transition between alphaC and alphaN.
  ## The formulation based on e^difference instead of the ratio
  ## also works for negative biomass synthesis
  ## but requires a scaling factor, tauB.
  ##seealso<< \code{\link{balanceAlphaBetweenCNLimitations}}
  # # min to avoid  + Inf
  wCLim = pmin( .Machine$double.xmax,
               exp( delta/tauB*(CsynBN - CsynBC)))
  wNLim = pmin( .Machine$double.xmax,
               exp( delta/tauB*(CsynBC - CsynBN)))
  alpha <- (wCLim*alphaC + wNLim*alphaN) / (wCLim + wNLim)
  alpha
}

getKmKNParsFromKmN <- function(
  ### get consistent parameters of the km and kN variant from SESAM kmN pars
  p            ##<< parameter list
  , kN = 60    ##<< turnover time of enzyme in 1/year
){
  pSeam <- within(p,{
    # in Sesam kmN = km * kN
    kN <- kNL <- kNR <- kN
    km <- kmL <- kmR <- kmN/kN
    kmN <- NULL
  })
  ##value<< parameter list with kmN replaced by components kNL, kNR, kmL, kmR
}

getSeamStateFromSesam <- function(
  ### get consistent state of the SEAM model from SESAM parameterization
  x0      ##<< state of SESAM
  ,parms  ##<< parameter list with entries aE, kNL and kNR
){
  c(x0
    , EL = unname((1 - x0["alpha"])*parms$aE*x0["B"]/parms$kNL)
    , ER = unname((x0["alpha"])*parms$aE*x0["B"]/parms$kNR)
    # make sure to have same order as derivative of Seam3
  )[c("B","ER","EL","R","RN","L","LN","IN","alpha")]
  ##value<< SEAM3 state vector with ER and EL set to quasi steady state
}

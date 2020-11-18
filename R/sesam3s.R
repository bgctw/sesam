#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam3s <- function(
  ### Soil Enzyme Steady Allocation model
  t,x,parms
){
  ##details<<
  ## deprecated, superseeded by sesam3a
  ## model corresponding to Seam3 with enzyme levels computed by quasi steady state
  ## Alpha as an explicit state variable that changes with turnover
  ## Simplified computation of target based on current revenue based on current alpha
  x <- pmax(unlist(x),1e-16)      # no negative masses
  # compute steady state enzyme levels for N and for C limitation
  decRp <- parms$kR * x["R"]
  decLp <- parms$kL * x["L"]
  cnR <- x["R"]/x["RN"]
  cnL <- x["L"]/x["LN"]
  cnE <- parms$cnE #alphaC*parms$cnER + (1-alphaC)*parms$cnEL
  cnB <- parms$cnB
  B <- x["B"]
  #ETot <- parms$aE * B / parms$kN
  rM <- parms$m * B          # maintenance respiration
  tvrB <- parms$tau*B        # microbial turnover
  if (isTRUE(parms$isEnzymeMassFlux)) {
    synE <- parms$aE * B       # total enzyme production per microbial biomass
    # growth respiration associated with enzyme production
    respSynE <- (1 - parms$eps)/parms$eps * synE
  } else {
    # neglect mass fluxes via enzymes
    synE <- 0
    respSynE <- 0
  }
  # for revenue, account for enzyme investments also if negleting mass fluxes
  synERev <- parms$aE * B
  # declare variables that will be computed/overidden in computeAlphaDependingFluxes
  # else operator '<<-' will override bindings in global environment
  tvrER <- tvrEL  <- decR <- decL <- tvrERecycling <- uC <-
    CsynBC <- decN <- plantNUp <- PhiNU <- immoNPot <- uNSubstrate <- uNPot <-
    NsynBN <- CsynBN <- isLimN <-
    CsynB <- PhiNB <- synB <- rG <- NA_real_
  computeAlphaDependingFluxes <- function(alpha){
    # compute fluxes that depend on alpha
    tvrER <<- alpha * synE
    tvrEL <<- (1 - alpha) * synE
    #
    decR <<- decRp * alpha*synERev/(parms$km*parms$kN + alpha*synERev)
    decL <<- decLp * (1 - alpha)*synERev/(parms$km*parms$kN + (1 - alpha)*synERev)
    #
    tvrERecycling <<- parms$kNB*(tvrER + tvrEL)
    uC <<- decR + decL + tvrERecycling
    # C required for biomass and growth respiration under C limitation
    CsynBC <<- uC - synE/parms$eps - rM
    #
    # Nitrogen balance
    decN <<- decR/cnR + decL/cnL + tvrERecycling/parms$cnE
    # plants get at maximum half of the decomposed organic N
    plantNUp <<- pmin(parms$plantNUp, decN/2)
    # mineralization due to soil heterogeneity (Manzoni 08)
    PhiNU <<- (1 - parms$nuN) * (decN - plantNUp)
    # immobilization flux
    immoNPot <<- parms$iBN * x["IN"]
    uNSubstrate <<- (decN - plantNUp - PhiNU)	# plant uptake also of organic N
    # potential N uptake from both inorganic and organic
    uNPot <<-  immoNPot + uNSubstrate
    # N required for biomass growth (after invest into enzymes)
    NsynBN <<- uNPot - synE/cnE
    # C required for biomass growth and associated growth resp under N limitation
    CsynBN <<- (NsynBN*cnB)/parms$eps
    #
    # whether is microbial N limited (with taking immobilization into account)
    isLimN <<- CsynBN <  CsynBC
    CsynB  <<- if (isLimN ) CsynBN else CsynBC
    if (CsynB > 0) {
      synB <<- parms$eps*CsynB
      rG <<- (1 - parms$eps)*CsynB
    } else {
      # with negative  biomass change, do not assign growth respiration
      synB <<- CsynB
      rG <<- 0
    }
    # imbalance mineralization/immobilization flux
    PhiNB <<- uNSubstrate - (synE/parms$cnE + synB/parms$cnB)
    alpha	##value<< return the given argument
  }
  #
  alpha <- computeAlphaDependingFluxes(x["alpha"])
  alphaC <- computeSesam3sAllocationPartitioning(
    dR = decRp, dL = decLp, B = B
    , kmkN = parms$km*parms$kN, aE = parms$aE
    , alpha = alpha
  )
  alphaN <- computeSesam3sAllocationPartitioning(
    dR = decRp/cnR, dL = decLp/cnL, B = B
    , kmkN = parms$km*parms$kN, aE = parms$aE
    , alpha = alpha
  )
  alphaTarget <- balanceAlphaBetweenCNLimitationsExp(
    alphaC, alphaN, CsynBN, CsynBC, tauB = parms$tau*B
    )
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) * (parms$tau + abs(synB)/B)
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC - (synE + respSynE + synB + rG + rM)
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
  plantNup <- parms$kINPlant*x["IN"]
  if (!is.null(parms$plantNUpAbs)) plantNup <- plantNup + parms$plantNUpAbs
  #
  dB <- synB - tvrB
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  # here plant uptake as absolute parameter
  dIN <-  +parms$iIN  - plantNup  - leachN  + PhiNU  + PhiNB  + PhiNTvr
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
  respB <- respSynE + rG + rM + respO
  resp <- respB + respTvr
  if (diff( unlist(c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps )  stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(uN = uNSubstrate, usage = synE/parms$cnE + synB/parms$cnB + PhiNB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dR + dL + tvrExC + resp,    parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + dRN + dLN + dIN + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL - plantNUp  + parms$iIN -
         plantNup - parms$lN*x["IN"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) { resDeriv["dR"] <- resDeriv["dRN"] <-  0   }
  if (isTRUE(parms$isFixedL)) { resDeriv["dL"] <- resDeriv["dLN"] <-  0   }
  if (isTRUE(parms$isFixedI)) { resDeriv["dIN"] <-  0   }
  #
  # further computations just for output for tacking the system
  ER <- alpha * parms$aE * B / parms$kN
  EL <- (1 - alpha) * parms$aE * B / parms$kN
  limER <- ER / (parms$kmR + ER)
  limEL <- EL / (parms$kmL + EL)
  revRC <- decRp / (parms$km*parms$kN + alphaC*synERev)
  revLC <- decLp / (parms$km*parms$kN + (1 - alphaC)*synERev)
  revRN <- decRp/cnR / (parms$km*parms$kN + alphaN*synERev)
  revLN <- decLp/cnL / (parms$km*parms$kN + alphaN*synERev)
  # net mic mineralization/immobilization when accounting uptake mineralization
  PhiNBU <- PhiNB + PhiNU
  # total mineralization flux including microbial turnover
  PhiNTotal <- PhiNBU + PhiNTvr
  # do not match in other limitation
  # c(alphaC, revRC/(revRC + revLC)); c(alphaN, revRN/(revRN + revLN))
  # compute C available for biomass, this time without accounting immobilization flux
  NsynBNSubstrate <- uNSubstrate - synE/cnE
  CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
  # N limited based on substrate uptake (without accounting for immobilization)
  isLimNSubstrate <-  CsynBNSubstrate < CsynBC
  #
  if (isTRUE(parms$isRecover) ) recover()
  list( resDeriv, c(
    respO = as.numeric(respO)
    , ER = as.numeric(ER), EL = as.numeric(EL)
    #, MmB = as.numeric(MmB)
    , PhiNB = as.numeric(PhiNB), PhiNU = as.numeric(PhiNU)
    , PhiNTvr = as.numeric(PhiNTvr)
    , PhiNBU = as.numeric(PhiNBU), PhiNTotal = as.numeric(PhiNTotal)
    , immoNPot = as.numeric(immoNPot)
    , alphaTarget = as.numeric(alphaTarget)
    , alphaC = as.numeric(alphaC), alphaN = as.numeric(alphaN)
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , limER = as.numeric(limER), limEL = as.numeric(limEL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , resp = as.numeric(resp), respB = as.numeric(respB)
    , respTvr = as.numeric(respTvr)
    , tvrB = as.numeric(tvrB)
    , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    , revRN = as.numeric(revRN)
    , revLN = as.numeric(revLN)
    , pCsyn = as.numeric(CsynBC / CsynBN), CsynReq = as.numeric(CsynBN)
    , Csyn = as.numeric(CsynBC)
    , pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    , NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    , dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    , dIN = as.numeric(dIN)
    , uC = as.numeric(uC), synB = as.numeric(synB)
    , decN = as.numeric(decN)
  ))
}


#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam4b <- function(
  ### Soil Enzyme Steady Allocation model with detailed turnover
  t,x,parms
){
  # derived from Sesam4a but with element-weighted investment
  xOrig <- x
  x <- pmax(unlist(x),1e-16)      # no negative masses
  # compute steady state enzyme levels for N and for C limitation
  dRPot <- parms$kR * x["RC"]
  dLPot <- parms$kL * x["LC"]
  dLPPot <- parms$kLP * x["LP"] # potential biomineralization
  dRPPot <- parms$kRP * x["RP"]
  immoNPot <- parms$iBN * x["IN"]
  immoPPot <- parms$iBP * x["IP"]
  cnR <- x["RC"]/x["RN"]
  cnL <- x["LC"]/x["LN"]
  cnE <- parms$cnE
  cnB <- parms$cnB
  cnBW <- parms$cnBW
  cnBL <- parms$cnBL
  cpR <- x["RC"]/x["RP"]
  cpL <- x["LC"]/x["LP"]
  cpE <- parms$cpE
  cpB <- parms$cpB
  cpBW <- parms$cpBW
  cW <- parms$cW
  if (((cW == 1) || (cW == 0)) && ((cnB != cnBW) || (cpB != cpBW))) stop(
    "If ceB != ceBW, organic turnover must be partitioned between R and DOM.",
    "cW=1 or cW=0 are not allowed.")
  cnBL <- if (cW == 1 ) cnB else (1 - cW)/(1/cnB - cW/cnBW)
  cpBL <- if (cW == 1 ) cpB else (1 - cW)/(1/cpB - cW/cpBW)
  alpha <- cbind(L = x["alphaL"], R = x["alphaR"], LP = NA, RP = x["alphaRP"])[1,]
  alpha["LP"] <- 1 - sum(alpha, na.rm = TRUE)
  B <- x["BC"]
  aeB <- parms$aE*B        # aeB without associanted growth respiration
  kmN <- parms$kmN #parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * alpha["L"]*aeB/(kmN + alpha["L"]*aeB)
  decR <- dRPot * alpha["R"]*aeB/(kmN + alpha["R"]*aeB)
  decLP <- dLPPot * alpha["LP"]*aeB/(kmN + alpha["LP"]*aeB)
  decRP <- dRPPot * alpha["RP"]*aeB/(kmN + alpha["RP"]*aeB)
  #
  # microbial turnover , partly feeds back to uptake
  tvrB <- parms$tau*B      # without predation/mineralization
  tvrBPred <- if (B < parms$B0) 0 else parms$tauP*(B - parms$B0)*B
  tvrBOrg <- tvrB + parms$epsP*tvrBPred
  tvrBMin <- (1 - parms$epsP)*tvrBPred
  #
  tvrERecycling <- parms$kNB*synE
  decNLR <- decL/cnL + decR/cnR
  decNE <- tvrERecycling/cnE
  decNB <- tvrBOrg*(1 - cW)/cnBL
  uNOrg <- parms$nuN*(decNLR + decNE + decNB)
  #uNOrg <- parms$nuN*(decL/cnL + decR/cnR + tvrERecycling/cnE +
  #tvrBOrg*(1 - cW)/cnBL)
  uPOrg <- parms$nuP*(decL/cpL + decR/cpR + tvrERecycling/cpE +
                        tvrBOrg*(1 - cW)/cpBL)
  uC <- decL + decR + tvrERecycling + tvrBOrg*(1 - cW)
  CsynBCt <- uC - rM - synE/parms$eps
  CsynBC <- if(CsynBCt > 0) parms$eps*CsynBCt else CsynBCt
  CsynBN <- cnB*(uNOrg + immoNPot - synE/cnE)
  CsynBP <- cpB*(uPOrg + immoPPot - synE/cpE)
  synB <- synB <- min(CsynBC, CsynBN, CsynBP)
  rG <- if (synB > 0) (1 - parms$eps)/parms$eps*synB else 0
  PhiNB <- uNOrg - synB/cnB - synE/cnE
  PhiPB <- uPOrg - synB/cpB - synE/cpE
  # community composition and enzyme allocation
  limE <- computeElementLimitations(
    cbind(C = CsynBC, N = CsynBN, P = CsynBP)[1,]
    , tauB = max(0.01,parms$tau)*B + tvrBPred) # avoid division by zero
  alphaTarget <- computeSesam4bAllocationPartitioning(
    dS = cbind(L = dLPot, R = dRPot)[1,]
    ,dSP = cbind(LP = dLPPot, RP = dRPPot)[1,]
    , B = B
    ,kmkN = kmN, aE =  parms$aE
    ,alpha = alpha
    ,limE = limE
    ,betaN = cbind(L = cnL, R = cnR, E = parms$cnE)[1,]
    ,betaP = cbind(L = cpL, R = cpR, E = parms$cpE)[1,]
  )
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) * (synB + tvrB + tvrBPred)/B
  #

  # imbalance fluxes of microbes and predators
  # (consuming part of microbial turnover)
  respO <- uC - (synE/parms$eps + synB + rG + rM)
  respTvr <- (1 - parms$epsP) * tvrBPred
  # assuming same cnRatio of predators to equal cn ratio of microbes
  PhiNTvr <- respTvr/parms$cnB
  PhiPTvr <- respTvr/parms$cpB
  #
  # tvr that feeds R pool, assume higher C:N ratios of cell walls
  tvrC <-  +cW*tvrBOrg   + (1 - parms$kNB)*synE
  tvrN <-  +cW*tvrBOrg/parms$cnBW   + (1 - parms$kNB)*synE/parms$cnE
  tvrP <-  +cW*tvrBOrg/parms$cpBW   + (1 - parms$kNB)*synE/parms$cpE
  # fluxes leaving the system (will be set in scen where trv does not feed back)
  tvrExC <- tvrExN <- tvrExP <- 0
  #
  leachN <- parms$lN*x["IN"]
  plantNUpPot <- parms$kINPlant*x["IN"]
  plantNUp <- if (!is.null(parms$plantNUpAbs)) {
    min(parms$plantNUpAbs, plantNUpPot)
  } else {
    plantNUpPot
  }
  PhiNU <- (1 - parms$nuN)*(decL/cnL + decR/cnR + tvrERecycling/cnE +
                            tvrBOrg*(1 - cW)/cnBL)
  leachP <- parms$lP*x["IP"]
  PhiPU <- (1 - parms$nuP)*(decL/cpL + decR/cpR + tvrERecycling/cpE +
                              tvrBOrg*(1 - cW)/cpBL)
  immoN <- max(0,-PhiNB); minN <- max(0,PhiNB)
  immoP <- max(0,-PhiPB); minP <- max(0,PhiPB)
  respB <- (synE)/parms$eps*(1 - parms$eps)  + rG + rM + respO
  resp <- respB + respTvr
  #
  dB <- synB - tvrB - tvrBPred
  if ((xOrig["BC"] <= 1e-16) && (dB < 0)) dB <- 0
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dLP <- -decL/cpL   + parms$iL/parms$cpIL
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  dRP <- -decR/cpR  + parms$iR/parms$cpIR  + tvrP
  # here plant uptake as absolute parameter
  dIN <-  +parms$iIN  - plantNUp  - leachN  + PhiNU  + PhiNB  + PhiNTvr
  dIP <-  +parms$iIP  - parms$kIPPlant*x["IP"]  - leachP  + PhiPU  + PhiPB  + PhiPTvr
  dResp <- resp
  dLeachN <- leachN
  dLeachP <- leachP
  #
  if (isTRUE(parms$isFixedS)) {
    # scenario of fixed substrate
    dR <- dL <- dRN <- dLN <- dRP <- dLP <- dIN <- dIP <- 0
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
    c(dB, dR, dRN, dRP, dL, dLN, dLP, dIN, dIP,
      dAlpha["L"], dAlpha["R"], dAlpha["RP"],
      dResp, dLeachN, dLeachP))
    ,names = c("BC","RC","RN","RP","LC","LN","LP","IN","IP"
               ,"alphaL","alphaR","alphaRP"
               , "resp", "leachN", "leachP"))[names(x)]
  if (any(!is.finite(resDeriv))) stop("encountered nonFinite derivatives")
  sqrEps <- sqrt(.Machine$double.eps)
  # parms$iL - (decL + dL)
  # parms$iR + tvrC -(decR + dR)
  #
  # checking the mass balance of fluxes
  # biomass mass balance
  if (diff( unlist(
    c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps )  stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(uN = uNOrg, usage = synE/parms$cnE
      + synB/parms$cnB + PhiNB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (diff( unlist(
    c(uP = uPOrg, usage = synE/parms$cpE
      + synB/parms$cpB + PhiPB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  # biomass turnover mass balance
  if (diff( unlist(
    c(tvrB + tvrBPred, usage = respTvr + tvrBOrg )))^2 > sqrEps )  stop(
    "biomass turnover mass balance C error")
  if (diff( unlist(c(
    tvrB/cnB + tvrBPred/cnB
    , usage = PhiNTvr + tvrBOrg/cnB
    #, usage = PhiNTvr + cW*tvrBOrg/cnBW + (1 - cW)*tvrBOrg/cnBL
  )))^2 > sqrEps )  stop(
    "biomass turnover mass balance N error")
  if (diff( unlist(c(
    (tvrB + tvrBPred)/cpB
    , usage = PhiPTvr + cW*tvrBOrg/cpBW + (1 - cW)*tvrBOrg/cpBL
     )))^2 > sqrEps )  stop("biomass turnover mass balance P error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dR + dL + tvrExC + dResp,    parms$iR + parms$iL )
      ))^2 >
      sqrEps )  stop("mass balance dC error")
    if (diff(unlist(
      c( dB/parms$cnB  + dRN + dLN + dIN + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL + parms$iIN -
         plantNUp - dLeachN)
      ))^2 > .Machine$double.eps )  stop("mass balance dN error")
    if (diff(unlist(
      c( dB/parms$cpB  + dRP + dLP + dIP + tvrExP
         , parms$iR/parms$cpIR  + parms$iL/parms$cpIL + parms$iIP -
         parms$kIPPlant*x["IP"] - dLeachP)))^2 >
      .Machine$double.eps )  stop("mass balance dP error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) {
    resDeriv["RC"] <- resDeriv["RN"] <- resDeriv["RP"] <- 0 }
  if (isTRUE(parms$isFixedL)) {
    resDeriv["LC"] <- resDeriv["LN"] <- resDeriv["LP"] <- 0 }
  if (isTRUE(parms$isFixedI)) {resDeriv["IN"] <- 0   }
  if (isTRUE(parms$isFixedIP)) {resDeriv["IP"]  <- 0   }
  #
  # further computations just for output for tacking the system
  limZ <- alpha * aeB / (parms$kmN + alpha * aeB)
  # net mic mineralization/immobilization when accounting uptake mineralization
  PhiNBU <- PhiNB + PhiNU
  # total mineralization flux including microbial turnover
  PhiNTotal <- PhiNBU + PhiNTvr
  # do not match in other limitation
  # compute C available for biomass, this time without accounting
  # immobilization flux
  NsynBNSubstrate <- uNOrg - synE/cnE
  CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
  # N limited based on substrate uptake (without accounting for immobilization)
  isLimNSubstrate <-  CsynBNSubstrate < CsynBC
  #
  if (isTRUE(parms$isRecover) ) recover()
  aux <- c(
    respTotal = as.numeric(resp)
    , respO = as.numeric(respO)
    , respB = as.numeric(respB)
    , respTvr = as.numeric(respTvr)
    #, ER = as.numeric(ER), EL = as.numeric(EL)
    #, MmB = as.numeric(MmB)
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
    , structure(limE, names = paste0("lim",names(limE)))
    , structure(limZ, names = paste0("limZ",names(limZ)))
    , structure(alphaTarget, names = paste0("alphaTarget",names(alphaTarget)))
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , cpR = as.numeric(cpR), cpL = as.numeric(cpL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , tvrB = as.numeric(tvrB)
    , tvrBPred = as.numeric(tvrBPred)
    , synB = as.numeric(synB) # may be negative
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    , CsynBP = as.numeric(CsynBP)
    , uptakeC = as.numeric(uC)
    , uptakeNOrg = as.numeric(uNOrg)
    , decNL = as.numeric(decL/cnL), decNR = as.numeric(decR/cnR)
    , decNE = as.numeric(decNE), decNB = as.numeric(decNB)
    #, pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    #, NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    #, dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    #, dIN = as.numeric(dIN)
    #, uC = as.numeric(uC), synB = as.numeric(synB)
    #, decN = as.numeric(decN)
  )
  list( resDeriv,  aux )
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
computeSesam4bAllocationPartitioning <- function(
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
  ,betaN    ##<< numeric vector of C/N ratios of substrates and enzymes
  ,betaP   ##<< numeric vector of C/P ratios of substrates (L,R) and enzymes (E)
){
  if (any(sort(c(names(dS),names(dSP))) != sort(names(alpha)))) stop(
    "expected names in c(dS, dSP) to match names(alpha) but got ",
    c(names(dS),names(dSP)), " versus ", names(alpha))
  p0 <- 0
  synEnz <- aE*B
  invest <- alpha*synEnz *
    (limE["C"] + limE["N"]/betaN["E"] + limE["P"]/betaP["E"])
  #cost <- (kmkN + alpha*synEnz) *
  #  (limE["C"] + limE["N"]/betaN["E"] + limE["P"]/betaP["E"])
  returnS <- structure(sapply(names(dS),function(S){
    dS[[S]] * alpha[[S]]*synEnz / (kmkN + alpha[[S]]*synEnz)*
      (limE["C"] + limE["N"]/betaN[[S]] + limE["P"]/betaP[[S]])
  }, USE.NAMES = FALSE), names = names(dS))
  returnSP <- structure(sapply(names(dSP),function(S){
    dSP[[S]]*((alpha[[S]]*synEnz + p0)/(kmkN + alpha[[S]]*synEnz + p0) -
                (p0)/(kmkN + p0))*limE["P"]
  }, USE.NAMES = FALSE), names = names(dSP))
  # make sure that ordering of return components is the same as in alpha
  revenue <- c(returnS, returnSP)[names(alpha)]/pmax(1e-12,invest)
  alphaTarget <- revenue/sum(revenue)
  alphaTarget
  ##value<< numeric vector of revenue-optimal alpha for enzymes L,R,LP,RP
}

.tmp.f <- function(){
  revL = dS["L"]/(kmkN + (1 - alpha["R"])*aE*B)
  revR = dS["R"]/(kmkN + alpha["R"]*aE*B)
  c(revL, revR)
  alphaTarget = revR/(revL + revR)

}


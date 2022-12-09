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
  dLPPot <- parms$kLP * x["LP"] # potential biomineralization attacking P part
  dRPPot <- parms$kRP * x["RP"]
  immoNPot <- parms$iBN * x["IN"]
  immoPPot <- parms$iBP * x["IP"]
  cnR <- x["RC"]/x["RN"]
  cnL <- x["LC"]/x["LN"]
  cnE <- parms$cnE
  cnB <- parms$cnB
  cnBW <- parms$cnBW
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
  alpha["LP"] <- max(0,1 - sum(alpha, na.rm = TRUE))
  if (any(alpha < -1e-8)) stop(
    "expected all alpha in [0..1] but was ", paste0(alpha,sep = ","))
  #alpha <- pmax(alpha,0)      # no negative community compositions
  B <- x["BC"]
  aeB <- parms$aE*B        # aeB without associanted growth respiration
  kmN <- parms$kmN #parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * alpha["L"]*aeB/(kmN + alpha["L"]*aeB)
  decR <- dRPot * alpha["R"]*aeB/(kmN + alpha["R"]*aeB)
  decLP <- dLPPot * (alpha["LP"]*aeB + parms$pELP)/(kmN + alpha["LP"]*aeB + parms$pELP)
  decRP <- dRPPot * (alpha["RP"]*aeB + parms$pERP)/(kmN + alpha["RP"]*aeB + parms$pERP)
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
  decP <- decL/cpL + decR/cpR + tvrERecycling/cpE
  decPP <- decLP + decRP
  uPOrg <- parms$nuP*(decP + tvrBOrg*(1 - cW)/cpBL)
  uPP <- parms$nuPP*(decPP)
  uC <- decL + decR + tvrERecycling + tvrBOrg*(1 - cW)
  CsynBCt <- uC - rM - synE/parms$eps
  CsynBC <- if(CsynBCt > 0) parms$eps*CsynBCt else CsynBCt
  CsynBN <- cnB*(uNOrg + immoNPot - synE/cnE)
  CsynBP <- cpB*(uPOrg + uPP + immoPPot - synE/cpE)
  synB <- synB <- min(CsynBC, CsynBN, CsynBP)
  rG <- if (synB > 0) (1 - parms$eps)/parms$eps*synB else 0
  PhiNB <- uNOrg - synB/cnB - synE/cnE
  PhiPB <- uPOrg + uPP - synB/cpB - synE/cpE
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
  #
  leachP <- parms$lP*x["IP"]
  PhiPU <-  (1 - parms$nuP)*(decP + tvrBOrg*(1 - cW)/cpBL) +
    (1 - parms$nuPP)*(decPP)
  immoN <- max(0,-PhiNB); minN <- max(0,PhiNB)
  immoP <- max(0,-PhiPB); minP <- max(0,PhiPB)
  respB <- (synE)/parms$eps*(1 - parms$eps)  + rG + rM + respO
  resp <- respB + respTvr
  #
  dB <- synB - tvrB - tvrBPred
  if ((xOrig["BC"] <= 1e-16) && (dB < 0)) dB <- 0
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dLP <- -decL/cpL   + parms$iL/parms$cpIL - decLP
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  dRP <- -decR/cpR  + parms$iR/parms$cpIR  + tvrP - decRP
  # here plant uptake as absolute parameter
  dIN <-  +parms$iIN  - plantNUp  - leachN  + PhiNU  + PhiNB  + PhiNTvr
  dIP <-  +parms$iIP  - parms$kIPPlant*x["IP"]  - leachP  +
    PhiPU  + PhiPB  + PhiPTvr
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
  if (any(!is.finite(resDeriv))) {
    stop("encountered nonFinite derivatives")
  }
  sqrEps <- sqrt(.Machine$double.eps)
  # parms$iL - (decL + dL)
  # parms$iR + tvrC -(decR + dR)
  #
  # checking the mass balance of fluxes
  # biomass mass balance
  if (diff( unlist(
    c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps )
    stop("biomass mass balance C error")
  if (diff( unlist(
    c(uN = uNOrg, usage = synE/parms$cnE
      + synB/parms$cnB + PhiNB )))^2 >
    .Machine$double.eps^{1/4}) {
    stop("biomass mass balance N error")
  }
  if (diff( unlist(
    c(uP = uPOrg + uPP, usage = synE/parms$cpE
      + synB/parms$cpB + PhiPB )))^2 >
    .Machine$double.eps)
    stop("biomass mass balance P error")
  # biomass turnover mass balance
  if (diff( unlist(
    c(tvrB + tvrBPred, usage = respTvr + tvrBOrg )))^2 > sqrEps )
    stop("biomass turnover mass balance C error")
  if (diff( unlist(c(
    tvrB/cnB + tvrBPred/cnB
    , usage = PhiNTvr + tvrBOrg/cnB
    #, usage = PhiNTvr + cW*tvrBOrg/cnBW + (1 - cW)*tvrBOrg/cnBL
  )))^2 > sqrEps )
    stop("biomass turnover mass balance N error")
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
      ))^2 > .Machine$double.eps^{1/4} ){
      stop("mass balance dN error")
    }
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


computeElementLimitations_sesam2 <- function(
  ### compute element limitations, superseded by computeElementLimitations
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
.tmp.f <- function(){
  w_sesam2 = computeElementLimitations_sesam2(CsynBE, tauB)
  wELimNorm - w_sesam2
}

computeElementLimitations <- function(
  ### compute element limitations
  CsynBE   ##<< numeric vector of carbon available for biomass synthesis
  , tauB   ##<< scalar typical microbial turnover flux for scaling
  , delta = 40  ##<< scalar smoothing factor, the higher, the steeper the transition
  , max_w = 12  ##<< maximum log(weight) to prevent numerical problems with large epxonents
){
  ##details<< Select the alpha corresponding to the smallest CsynBE.
  ## However, if elemental limitations are close,
  ## do a smooth transition between corresponding smallest alpha values
  synB = min(CsynBE)
  wELim <- exp(pmin(max_w, -delta/tauB*(CsynBE - synB)) )
  wELimNorm <- structure(wELim/sum(wELim), names=names(CsynBE))
  ##value<< numeric vector: fractional element limitations with names
  ## corresponding to names of CsynBE
  wELimNorm
}

#' @export
computeSesam4bAllocationPartitioning <- function(
  ###  allocation partitioning alpha of four enzymes and biomineralization
  dS    ##<< numeric vector (L,R) of potential depolymerization C-fluxes
  ,dSP  ##<< numeric vector of potential biomineralization P-fluxes
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
  ,e_P = 0
){
  #limE <- limEP <- c(C = 0, N = 0, P = 1)
  #limE <- limE0 <- c(C = 0, N = 0.408838431792989, P = 0.591161568206928  )
  #limE <- limE0 <- c(C = 0, N = 0.226012593299406, P = 0.773987406700594  )
  #limE <- limE0 <- c(C = 0.000435098291453755, N = 0.563867465009725, P = 0.435697436698822  )
  #limE["N"] <- limE["N"] + limE["P"]; limE["P"] <- 0
  #limE <- c(C = 0, N = 0.563867465009725, P = 0.435697436698822)
  #limE <- limEN <- c(C = 0, N = 1, P = 0)
  if (any(sort(c(names(dS),names(dSP))) != sort(names(alpha)))) stop(
    "expected names in c(dS, dSP) to match names(alpha) but got ",
    c(names(dS),names(dSP)), " versus ", names(alpha))
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
    dSP[[S]]*((alpha[[S]]*synEnz + e_P)/(kmkN + alpha[[S]]*synEnz + e_P) -
                (e_P)/(kmkN + e_P))*limE["P"]
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
  (alphaTarget_old = revR/(revL + revR))
}


compute_eweighted_potential <- function(
  dS,    ##<< numeric vector (L,R) of potential depolymerization C-fluxes
  dSP,   ##<< numeric vector of potential biomineralization P-fluxes
  limE,      ##<< elemental limiations
  betaN,    ##<< numeric vector of C/N ratios of substrates and enzymes
  betaP   ##<< numeric vector of C/P ratios of substrates (L,R) and enzymes (E)
) {
  c(
    dS["L"]*(limE["C"]+limE["N"]/betaN["L"]+limE["P"]/betaP["L"]),
    dS["R"]*(limE["C"]+limE["N"]/betaN["R"]+limE["P"]/betaP["R"]),
    #P = unname((dSP["L"]/betaP["L"]+dSP["R"]/betaP["R"])*limE["P"])
    # dSP already in P units
    P = unname(limE["P"] * (dSP["L"]+dSP["R"]))
  )
}

#' @export
computeSesam4bOptimalAllocationPartitioning <- function(
  ###  allocation partitioning alpha of four enzymes and biomineralization
  dL, dR, dP, ##<< elemental-limitation weighted potential decomposition
  params,    ##<< list of model parameters with entries kmN, aE, e_P
  B, synB
){
  p = within(params, B<-B, synB <- synB)
  dumax = sort(du_dalpha(c(L=0,R=0,P=0), dL, dR, dP, p), decreasing = TRUE)
  if(names(dumax[1]) == "L") {
    du1 = du_dalphaS(1, dL, p, e_P=0)
    if (du1 > dumax[2]) return(c(L=1,R=0,P=0))
    if (names(dumax[2]) == "R"){
      alphaL = calc_alpha3_optLR(dL, dR, p)
      du2 = du_dalphaS(alphaL, dL, p, e_P = 0)
      if (du2 > dumax[3]) return(c(L=alphaL,R=(1-alphaL),P=0))
    } else { # umax2 must be P
      alphaL = calc_alphaS_optSP(dL,dP,p)
      du2 = du_dalphaS(alphaL, dL, p, e_P = 0)
      if (du2 > dumax[3]) return(c(L=alphaL,R=0,P=(1-alphaL)))
    }
    alpha3 = calc_alpha3_optLRP(dL,dR,dP,p)
    return(alpha3)
  }
  if(names(dumax[1]) == "R") {
    du1 = du_dalphaS(1, dR, p, e_P = 0)
    if (du1 > dumax[2]) return(c(L=0,R=1,P=0))
    if (names(dumax[2]) == "L"){
      alphaL = calc_alpha3_optLR(dL, dR, p)
      du2 = du_dalphaS(alphaL, dL, p, e_P = 0)
      if (du2 > dumax[3]) return(c(L=alphaL,R=(1-alphaL),P=0))
    } else { # umax2 must be P
      alphaR = calc_alphaS_optSP(dR,dP,p)
      du2 = du_dalphaS(alphaR, dR, p, e_P = 0)
      if (du2 > dumax[3]) return(c(L=0,R=alphaR,P=(1-alphaR)))
    }
    alpha3 = calc_alpha3_LRP(dL,dR,dP,p)
    return(alpha3)
  }
  if(names(dumax[1]) == "P") {
    du1 = du_dalphaS(1, dP, p)
    if (du1 > dumax[2]) return(c(L=0,R=0,P=1))
    if (names(dumax[2]) == "L"){
      alphaL = calc_alphaS_optSP(dL,dP,p)
      du2 = du_dalphaS(alphaL, dL, p, e_P = 0)
      if (du2 > dumax[3]) return(c(L=alphaL,R=0,P=(1-alphaL)))
    } else if (names(dumax[2]) == "R") {
      alphaR = calc_alphaS_optSP(dR,dP,p)
      du2 = du_dalphaS(alphaR, dR, p, e_P = 0)
      if (du2 > dumax[3]) return(c(L=0,R=alphaR,P=(1-alphaR)))
    }
    alpha3 = calc_alpha3_optLRP(dL,dR,dP,p)
    return(alpha3)
  }
  stop("unhandled umax: ", dumax)
  ##value<< numeric vector of revenue-optimal alpha for enzymes L,R,LP,RP
}

du_dalpha <- function(
  ### compute derivatives of total return given an allocation
  alpha3,
  dL, dR, dP, ##<< elemental-limitation weighted potential decomposition
  p  ## list with entries kmN, aE, e_P, and B
){
  if (is.matrix(alpha3)) {
    du <- cbind(
      du_dalphaS(alpha3[,"L", drop=FALSE], dL, p, e_P = 0),
      du_dalphaS(alpha3[,"R", drop=FALSE], dR, p, e_P = 0),
      du_dalphaS(alpha3[,"P", drop=FALSE], dP, p)
    )
    colnames(du) <- c("L","R","P")
    du
  } else {
    c(
      L = unname(du_dalphaS(alpha3["L"], dL, p, e_P = 0)),
      R = unname(du_dalphaS(alpha3["R"], dR, p, e_P = 0)),
      L = unname(du_dalphaS(alpha3["P"], dP, p))
    )
  }
}
du_dalphaS <- function(
  ### compute derivatives of total return given an allocation
  alpha, ##<< allocation to specific enzyme
  dS,    ##<< elemental-limitation weighted potential decomposition
  p,     ##<< list with entries kmN, aE, e_P, and B
  e_P = p$e_P
){
  ##details<< Set e_P to zero for depolymerizing enzymes
  aeB = p$aE * p$B
  aeB * p$kmN * dS / (e_P + p$kmN + alpha*aeB)^2
}

calc_alpha3_optLR <- function(dL, dR, p, alphaP=0){
  dLmdR = dL - dR
  alphaL <- unname(ifelse(dLmdR == 0, 0.5, {
    aeB = p$aE*p$B
    A = aeB*(1-alphaP) * dL + p$kmN*(dL+dR)
    D = (aeB*(1-alphaP) + 2*p$kmN)*sqrt(dL*dR)
    d = aeB*dLmdR
    alpha0 <- (A-D)/d  # only the second root is reasonable
    alpha <- pmin(1,pmax(0,alpha0))
  }))
  alphaL
}
calc_alpha3_optLRP <- function(dL,dR, dP, p){
  aeB = p$aE*p$B
  A = 2*aeB*dL^(3/2)*dP*sqrt(dR) - aeB*dL^2*dP - aeB*dL*dP*dR + 4*dL^(3/2)*dP*sqrt(dR)*p$kmN - dL^3*p$e_P - dL^3*p$kmN -
    2*dL^2*dP*p$kmN + 2*dL^2*dR*p$e_P + 2*dL^2*dR*p$kmN - 2*dL*dP*dR*p$kmN - dL*dR^2*p$e_P - dL*dR^2*p$kmN
  D = sqrt(dP)*(aeB + p$e_P + 3*p$kmN)*sqrt(-2*dL^(9/2)*sqrt(dR) + 4*dL^(7/2)*dR^(3/2) - 2*dL^(5/2)*dR^(5/2) + dL^5 - dL^4*dR - dL^3*dR^2 + dL^2*dR^3)
  B = aeB*(2*dL^(3/2)*dP*sqrt(dR) + dL^3 - dL^2*dP - 2*dL^2*dR - dL*dP*dR + dL*dR^2)
  alphaP <- (A+D)/B
  alphaL <- calc_alpha3_optLR(dL,dR,p,alphaP)
  alphaR <- 1-alphaL-alphaP
  alpha=cbind(L=alphaL,R=alphaR,P=alphaP)
  alpha
}
calc_alphaS_optSP <- function(dS, dP, p){
  dLmdR = dS - dP
  alphaL <- unname(ifelse(dLmdR == 0, 0.5, {
    aeB = p$aE*p$B
    A = (aeB+p$e_P+p$kmN)*dL + p$kmN*dP
    D = sqrt(dL*dP)*(aeB + p$e_P + 2*p$kmN)
    d = aeB*dLmdR
    alpha0 <- (A-D)/d  # only the second root is reasonable
    alphaL <- pmin(1,pmax(0,alpha0))
  }))
  alphaL
}




u_decomp <- function(
  ### compute the derivative of revenue of a depolymerizing enzyme
  alpha, ##<< allocation to this enzyme
  dLorR, ##<< elemental-limitation weighted potential decomposition
  p,      ##<< list with entries kmN, aE and B
  dE=1  ##<< elemental-limitation weighted enzyme investment
){
  dLorR/dE/(p$kmN + alpha * p$aE * p$B)
}
u_biomin <- function(alpha, dP, p, dE=1){
  dP*p$kmN/dE/((p$kmN + p$e_P)^2 + alpha * p$aE * p$B * ((p$kmN + p$e_P)))
}





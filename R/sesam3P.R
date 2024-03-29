#library(deSolve)

# gC/m2 and gN/m2, /yr

#' @export
derivSesam3P <- function(
  ### Soil Enzyme Steady Allocation model including phosporous dynamics
  t,x,parms
){
  ##details<<
  ## derivSesam3P (Sesam3a extended by Phosphorous) extended by biomineralization
  ## Compared to Sesam3P, here
  ## - the biomineralizing enzmye attacks both L and R pool
  ## - the biomineralizing enzyme is also produced by plants at rate e_P
  x_neg <- x # remember original potential negative values
  x <- pmax(unlist(x),1e-16)      # no negative masses
  # compute steady state enzyme levels for N and for C limitation
  dRPot <- parms$kR * x["R"]
  dLPot <- parms$kL * x["L"]
  immoNPot <- parms$iBN * x["IN"]
  immoPPot <- parms$iBP * x["IP"]
  cnR <- x["R"]/x["RN"]
  cnL <- x["L"]/x["LN"]
  cnE <- parms$cnE
  cnB <- parms$cnB
  cpR <- x["R"]/x["RP"]
  cpL <- x["L"]/x["LP"]
  cpE <- parms$cpE
  cpB <- parms$cpB
  if (is.null(parms$cpm) ) stop("need to provide parameter cpm (C:P for which dec_LP decreases to 1/2.")
  lim_LP <- 1/(1 + cpL/parms$cpm)
  lim_RP <- 1/(1 + cpR/parms$cpm)
  dLPPot <- parms$kLP * x["LP"] * lim_LP # potential biomineralization in P units
  dRPPot <- parms$kRP * x["RP"] * lim_RP
  # one state variable less, here P, because alphas sum to one
  alpha <- cbind(L = x["alphaL"], R = x["alphaR"], P = NA)[1,]
  # take care for small negative alphas and renormalize to sum to 1
  alpha["P"] <- max(0, 1 - sum(alpha, na.rm = TRUE))
  alpha <- alpha/sum(alpha)
  B <- x["B"]
  aeB <- parms$aE*B        # aeB without associated growth respiration
  kmNL <- if (is.null(parms$kmNL)) parms$kmN else parms$kmNL
  kmNR <- if (is.null(parms$kmNR)) parms$kmN else parms$kmNR
  kmNP <- if (is.null(parms$kmNP)) parms$kmN else parms$kmNP
  #kmN <- parms$kmN #parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  tvrB <- parms$tau*B      # microbial turnover
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * alpha["L"]*aeB/(kmNL + alpha["L"]*aeB)
  decR <- dRPot * alpha["R"]*aeB/(kmNR + alpha["R"]*aeB)
  limEnzP <- (parms$e_P + alpha["P"]*aeB)/(kmNP + parms$e_P + alpha["P"]*aeB)
  decLP_P <- dLPPot * limEnzP # P units
  decRP_P <- dRPPot * limEnzP
  #
  tvrERecycling <- parms$kNB*synE
  uNOrg <- parms$nuN*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  uPOrg <- parms$nuP*(decL/cpL + decR/cpR + tvrERecycling/cpE)
  uC <- decL + decR + tvrERecycling
  CsynBCt <- uC - rM - synE/parms$eps
  CsynBC <- if (CsynBCt > 0) parms$eps*CsynBCt else CsynBCt
  CsynBN <- cnB*(uNOrg + immoNPot - synE/cnE)
  CsynBP <- cpB*(uPOrg + immoPPot - synE/cpE)
  synB <- min(CsynBC, CsynBN, CsynBP)
  if (synB > 0) {
    rG <- (1 - parms$eps)/parms$eps * synB
  } else {
    rG <- 0 # with negative biomass change, do growth respiration
  }
  PhiNB <- uNOrg - synB/cnB - synE/cnE
  PhiPB <- uPOrg - synB/cpB - synE/cpE
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC - (synE/parms$eps + synB + rG + rM)
  respTvr <- (1 - parms$epsTvr) * tvrB
  # assuming same cnRatio of predators to equal cn ratio of microbes
  PhiNTvr <- respTvr/parms$cnB
  PhiPTvr <- respTvr/parms$cpB
  #
  # tvr that feeds R pool, assume that N in SOM for resp (by epsTvr) is mineralized
  tvrC <-  +parms$epsTvr*tvrB   + (1 - parms$kNB)*synE
  tvrN <-  +parms$epsTvr*tvrB/parms$cnB   + (1 - parms$kNB)*synE/parms$cnE
  tvrP <-  +parms$epsTvr*tvrB/parms$cpB   + (1 - parms$kNB)*synE/parms$cpE
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
  plantPUpPot <- parms$kIPPlant*x["IP"]
  plantPUp <- if (!is.null(parms$plantPUpAbs)) {
    min(parms$plantPUpAbs, plantPUpPot)
  } else {
    plantPUpPot
  }
  PhiNU <- (1 - parms$nuN)*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  leachP <- parms$lP*x["IP"]
  PhiPU <- (1 - parms$nuP)*(decL/cpL + decR/cpR + tvrERecycling/cpE)
  #
  # community composition and enzyme allocation
  limE <- computeElementLimitations(
    cbind(C = CsynBC, N = CsynBN, P = CsynBP)[1,]
    , tauB = parms$tau*B
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )
  p_uNmic <- immoNPot/(immoNPot + plantNUp)
  p_uPmic <- immoPPot/(immoPPot + plantPUp)
  nuTZ <- if (isTRUE(parms$useOneNuT)) c(C=1, N=1, P=1) else c(
    C = parms$eps,
    N = unname(parms$nuN+(1-parms$nuN)*p_uNmic),
    P = unname(parms$nuP+(1-parms$nuP)*p_uPmic)
  )
  dAlpha <- if (isTRUE(parms$isFixedAlpha)) {
    c(L=0, R=0, P=0)
  } else if (isTRUE(parms$isRelativeAlpha)) {
    dAlpha_rel <- calc_dAlphaP_relative_plant(
      alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
      cnL, cnR, parms$cnB, cpL, cpR, parms$cpB,
      kmNZ = c(L=kmNL, R = kmNR, P = kmNP),
      nuTZ = nuTZ
    )
  } else if (isTRUE(parms$isOptimalAlpha)) {
    res_dAlpha_opt <- calc_dAlphaP_optimal(
      alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
      cnL, cnR, parms$cnB, cpL, cpR, parms$cpB,
      nuTZ = nuTZ
    )
    res_dAlpha_opt[[1]]
  } else {
    res_dAlpha <- calc_dAlphaP_propto_du(
      alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
      cnL, cnR, parms$cnB, cpL, cpR, parms$cpB,
      kmNZ = c(L=kmNL, R = kmNR, P = kmNP),
      nuTZ = nuTZ
    )
    res_dAlpha[[1]]
  }
  if (x["alphaR"] <= 1e-16){
    tmp = 1
  }
  # make sure to not decrease already negative alpha - workaround solver
  dAlpha[c("L","R")] = ifelse(
    x_neg[c("alphaL","alphaR")] < 0,
    -x_neg[c("alphaL","alphaR")], dAlpha[c("L","R")])
  dB <- synB - tvrB
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dLP <- -decL/cpL -decLP_P + parms$iL/parms$cpIL
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  dRP <- -decR/cpR  -decRP_P + parms$iR/parms$cpIR  + tvrP
  # here plant uptake as absolute parameter
  dIN <-  +parms$iIN  - plantNUp  - leachN  + PhiNU  + PhiNB  + PhiNTvr
  dIP <-  +parms$iIP  - plantPUp  - leachP  + PhiPU  + PhiPB  + PhiPTvr +
    decLP_P +decRP_P
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
    c( dB, dR, dRN, dRP, dL, dLN, dLP, dIN, dIP
       , dAlpha["L"], dAlpha["R"]))
    ,names = c("dB","dR","dRN","dRP","dL","dLN","dLP","dIN","dIP"
               ,"dAlphaL","dAlphaR"))
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
  if (diff( unlist(
    c(uP = uPOrg, usage = synE/parms$cpE + synB/parms$cpB + PhiPB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance P error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dR + dL + tvrExC + resp,    parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + dRN + dLN + dIN + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL + parms$iIN -
         plantNUp - parms$lN*x["IN"])))^2 >
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
  if (isTRUE(parms$isFixedI)) { resDeriv["dIN"] <-  resDeriv["dIP"] <- 0   }
  #
  # further computations just for output for tacking the system
  limZ <- alpha * aeB / (parms$kmN + alpha * aeB)
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
  list( resDeriv,  c(
     resp = as.numeric(resp)
     , respO = as.numeric(respO)
     , respB = as.numeric(respB)
     , respTvr = as.numeric(respTvr)
     #, ER = as.numeric(ER), EL = as.numeric(EL)
    #, MmB = as.numeric(MmB)
    , PhiNTotal = as.numeric(PhiNTotal)
    , PhiNB = as.numeric(PhiNB), PhiNU = as.numeric(PhiNU)
    , PhiNTvr = as.numeric(PhiNTvr)
    , PhiNBU = as.numeric(PhiNBU)
    , plantNUp = as.numeric(plantNUp)
    , immoNPot = as.numeric(immoNPot)
    , PhiPTotal = as.numeric(PhiPB + PhiPU + PhiPTvr)
    , PhiPB = as.numeric(PhiPB), PhiPU = as.numeric(PhiPU)
    , PhiPTvr = as.numeric(PhiPTvr)
    , PhiPBU = as.numeric(PhiPB + PhiPU)
    , immoPPot = as.numeric(immoPPot)
    , structure(limE, names = paste0("lim",names(limE)))
    , structure(limZ, names = paste0("limZ",names(limZ)))
    #, structure(alphaTarget, names = paste0("alphat",names(alphaTarget)))
    , alphaP = alpha[["P"]]
    , dAlphaR = dAlpha[["R"]]
    , dAlphaR2 = resDeriv[["dAlphaR"]]
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , cpR = as.numeric(cpR), cpL = as.numeric(cpL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , tvrB = as.numeric(tvrB)
    , synB = as.numeric(synB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    , CsynBP = as.numeric(CsynBP)
    #, pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    #, NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    #, dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    #, dIN = as.numeric(dIN)
    , uptakeC = as.numeric(uC)
    #, decN = as.numeric(decN)
    , decLP_P = as.numeric(decLP_P), decRP_P = as.numeric(decRP_P)
    , lim_LP = unname(lim_LP)
    , lim_RP = unname(lim_RP)
  ))
}

# moved computation of elemental limitations to sesam4b

calc_dAlphaP_relative_plant <- function(
    ### compute time derivative of alpha with Relative approach
  alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
  cnL, cnR, cnB, cpL, cpR, cpB,
  kmNZ = c(L=parms$kmN, R = parms$kmN, P = parms$kmN),
  nuTZ = setNames(rep(1.0, length(limE)), names(limE))  ##<< proportion of biomass synthesis
){
  alphaTarget <- computeSesam4bAllocationPartitioning(
    dS = cbind(R = dRPot, L = dLPot)[1,]
    ,dSP = c(P = as.numeric(dRPPot+dLPPot))
    , B = B
    ,kmkN = parms$kmN, aE =  parms$aE
    ,alpha = alpha
    ,limE = limE
    ,betaN = c(L = unname(cnL), R = unname(cnR), B=unname(cnB), E = parms$cnE)
    ,betaP = c(L = unname(cpL), R = unname(cpR), B=unname(cpB), E = parms$cpE)
    ,e_P = parms$e_P
    ,kmNZ = kmNZ
    ,nuTZ = nuTZ
  )
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) *  (parms$tau + abs(synB)/B)
}

calc_dAlphaP_optimal <- function(
    ### compute time derivative of alpha with Optimal approach
    alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
    cnL, cnR, cnB, cpL, cpR, cpB,
    nuTZ = setNames(rep(1.0, length(limE)), names(limE))  ##<< proportion of biomass synthesis
) {
  betaB <- c(C=1, N=cnB, P=cpB)
  omega_L <- compute_elemental_weightfactor(limE, c(C=1, N=cnL, P=cpL), betaB, nuTZ=nuTZ)
  omega_R <- compute_elemental_weightfactor(limE, c(C=1, N=cnR, P=cpR), betaB, nuTZ=nuTZ)
  omega_P <- compute_elemental_weightfactor(limE["P"], c(P=1), betaB["P"], nuTZ=nuTZ[["P"]])
  # dSw <- compute_eweighted_potential(
  #   dS = cbind(R = dRPot, L = dLPot)[1,]
  #   ,dSP = c(L = unname(dLPPot), R=unname(dRPPot))
  #   ,limE = limE
  #   ,betaN = cbind(L = cnL, R = cnR, cnB= cnB, E = parms$cnE)[1,]
  #   ,betaP = cbind(L = cpL, R = cpR, cpB= cpB, E = parms$cpE)[1,]
  # )
  alphaTarget <- computeSesam4bOptimalAllocationPartitioning(
    dLPot*omega_L, dRPot*omega_R, (dLPPot+dRPPot)*omega_P,
    params=parms, B, synB)
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) *  (parms$tau + abs(synB)/B)
  list(dAlpha=dAlpha, alphaTarget = alphaTarget)
}

# calc_dAlphaP_propto_du_alpha <- function(
#     ### changes in community additionally being proportional to alpha
#     alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
#     cnL, cnR, cpL, cpR,
#     kmNZ = c(L=parms$kmN, R = parms$kmN, P = parms$kmN)
# ){
#   dL <- dLPot * (limE["C"] + limE["N"]/cnL + limE["P"]/cpL)
#   dR <- dRPot * (limE["C"] + limE["N"]/cnR + limE["P"]/cpR)
#   #dP <- limE["P"]*(dLPPot/cpL+dRPPot/cpR)
#   dP <- limE["P"]*(dLPPot+dRPPot) # potential fluxes here already in P units
#   aeB <- parms$aE*B
#   du <- c(
#     L = unname(aeB*kmNZ[["L"]]*dL/(kmNZ[["L"]] + alpha["L"]*aeB)^2),
#     R = unname(aeB*kmNZ[["R"]]*dR/(kmNZ[["R"]] + alpha["R"]*aeB)^2),
#     P = unname(aeB*kmNZ[["P"]]*dP/(parms$e_P + kmNZ[["P"]] + alpha["P"]*aeB)^2)
#   )
#   mean_du <- sum(alpha*du) #weighted.mean(du,alpha)
#   dud = alpha*(du - mean_du)/mean_du
#   dalpha = (parms$tau + abs(synB)/B) * dud
#   if (abs(sum(dalpha)) > 1e-12) {
#     tmp = 1
#   }
#   dalpha
#   list(dalpha=dalpha, dud=dud, du=du, dS=c(L=unname(dL),R=unname(dR),P=unname(dP)),is=is)
# }

calc_dAlphaP_propto_du <- function(
    ### compute time derivative of alpha with Derivative approach
  alpha, dRPot, dLPot, dRPPot, dLPPot, synB, B, parms, limE,
    cnL, cnR, cnB, cpL, cpR, cpB,
    kmNZ = c(L=parms$kmN, R = parms$kmN, P = parms$kmN),
    nuTZ = setNames(rep(1.0, length(limE)), names(limE))  ##<< proportion of biomass synthesis
) {
  # # try expressing elemental-weighted potential in biomass C units (*ceB)
  # dL <- dLPot * (limE["C"] + limE["N"]/cnL*cnB + limE["P"]/cpL*cpB)
  # dR <- dRPot * (limE["C"] + limE["N"]/cnR*cnB + limE["P"]/cpR*cpB)
  # #dP <- limE["P"]*(dLPPot/cnL+dRPPot/cpR)
  # #dLPPot already in P units
  # dP <- limE["P"]*(dLPPot+dRPPot)*cpB
  betaB <- c(C=1, N=cnB, P=cpB)
  omega_L <- compute_elemental_weightfactor(limE, c(C=1, N=cnL, P=cpL), betaB, nuTZ=nuTZ)
  omega_R <- compute_elemental_weightfactor(limE, c(C=1, N=cnR, P=cpR), betaB, nuTZ=nuTZ)
  omega_P <- compute_elemental_weightfactor(limE["P"], c(P=1), betaB["P"], nuTZ=nuTZ[["P"]])
  dL <- dLPot * omega_L
  dR <- dRPot * omega_R
  #dLPPot already in P units
  dP <- (dLPPot+dRPPot) * omega_P
  aeB <- parms$aE*B
  du <- c(
    L = unname(aeB*kmNZ[["L"]]*dL/(kmNZ[["L"]] + alpha["L"]*aeB)^2),
    R = unname(aeB*kmNZ[["R"]]*dR/(kmNZ[["R"]] + alpha["R"]*aeB)^2),
    P = unname(aeB*kmNZ[["P"]]*dP/(parms$e_P + kmNZ[["P"]] + alpha["P"]*aeB)^2)
  )
  # only change alpha if its larger than zero or if the change is positive
  # update the mean_du to only include the participating enzymes
  mean_du_prev <- -Inf; mean_du <- mean(du)
  Z0 <- setNames(rep(FALSE, length(alpha)), names(alpha)) # enzymes not in Z0
  while (mean_du_prev != mean_du) {
    Z0 <- Z0 | (du - mean_du < -alpha*mean_du)
    mean_du_prev <- mean_du
    mean_du <- sum(du[!Z0])/(sum(!Z0) + sum(alpha[Z0]))
  }
  dud = (du - mean_du)/mean_du
  #the formulation for optimal allocation is on an absolute scale rather than
  #a relative scale. Putting optimal allocation to a relative scale would
  #result in deviding (alpha_Target - alpha) by mean(alpha), which is 1/n_enz
  #because sum(alpha)=1. To compensate, here we multiply by (1/n_enz = 1/sum(is))
  #dud = (du - mean_du)/mean_du/sum(is)
  dud[Z0] <- -alpha[Z0]
  dalpha = (parms$tau + abs(synB)/B) * dud
  dalpha
  list(dalpha=dalpha, dud=dud, du=du, dS=c(L=unname(dL),R=unname(dR),P=unname(dP)),Z0=Z0)
}



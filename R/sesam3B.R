#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam3B <- function(
  ### Soil Enzyme Steady Allocation model with biomass in QSS
  t,x,parms
){
  ##details<<
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
  #ETot <- parms$aE * B / parms$kN
  # potential immobilization flux
  immoPot <- parms$iB * x["I"]
  # steady state biomass
  B <- sesam3BSteady(dLC = decLp, dRC = decRp
                     , dLN = decLp/cnL, dRN = decRp/cnR
                     , alpha = x["alpha"], parms = parms
                     , immNPot = immoPot
                     )
  rM <- parms$m * B          # maintenance respiration
  tvrB <- parms$tau*B        # microbial turnover
  # neglect mass fluxes via enzymes
  synE <- 0
  respSynE <- 0
  # for revenue, account for enzyme investments also if negleting mass fluxes
  synERev <- parms$aE * B
  # compute fluxes that depend on alpha
  alpha <- x["alpha"]
  tvrER <- alpha * synE
  tvrEL <- (1 - alpha) * synE
  #
  decR <- decRp * alpha*synERev/(parms$km*parms$kN + alpha*synERev)
  decL <- decLp * (1 - alpha)*synERev/(parms$km*parms$kN + (1 - alpha)*synERev)
  #
  tvrERecycling <- parms$kNB*(tvrER + tvrEL)
  uC <- decR + decL + tvrERecycling
  # C required for biomass and growth respiration under C limitation
  CsynBC <- uC - synE/parms$eps - rM
  #
  # Nitrogen balance
  decN <- decR/cnR + decL/cnL + tvrERecycling/parms$cnE
  # plants get at maximum half of the decomposed organic N
  plantNUp <- pmin(parms$plantNUp, decN/2)
  # mineralization due to soil heterogeneity (Manzoni 08)
  PhiU <- (1 - parms$nu) * (decN - plantNUp)
  uNSubstrate <- (decN - plantNUp - PhiU)	# plant uptake also of organic N
  # potential N uptake from both inorganic and organic
  uNPot <-  immoPot + uNSubstrate
  # N required for biomass growth (after invest into enzymes)
  NsynBN <- uNPot - synE/cnE
  # c(NsynBN, parms$tau/parms$cnB*B)
  # C required for biomass growth and associated growth resp under N limitation
  CsynBN <- (NsynBN*cnB)/parms$eps
  #
  # whether is microbial N limited (with taking immobilization into account)
  isLimN <- CsynBN <  CsynBC
  CsynB  <- if (isLimN ) CsynBN else CsynBC
  if (CsynB > 0) {
    synB <- parms$eps*CsynB
    rG <- (1 - parms$eps)*CsynB
  } else {
    # with negative  biomass change, do not assign growth respiration
    synB <- CsynB
    rG <- 0
  }
  # imbalance mineralization/immobilization flux
  PhiB <- uNSubstrate - (synE/parms$cnE + synB/parms$cnB)
  #
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
  alphaTarget <- balanceAlphaBetweenCNLimitations(
    alphaC, alphaN, CsynBN, CsynBC
    , NsynBC = parms$eps*CsynBC/cnB, NsynBN)
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) * parms$tau
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC - (synE + respSynE + synB + rG + rM)
  respTvr <- (1 - parms$epsTvr) * tvrB
  # assuming same cnRatio of predators to equal cn ratio of microbes
  PhiTvr <- respTvr/parms$cnB
  #
  # tvr that feeds R pool, assume that N in SOM for resp (by epsTvr) is mineralized
  tvrC <-  +parms$epsTvr*tvrB   + (1 - parms$kNB)*(tvrER  + tvrEL)
  tvrN <-  +parms$epsTvr*tvrB/parms$cnB   + (1 - parms$kNB)*(tvrER  + tvrEL)/parms$cnE
  # fluxes leaving the system (will be set in scen where trv does not feed back)
  tvrExC <- tvrExN <- 0
  #
  leach <- parms$l*x["I"]
  #
  dB <- synB - tvrB
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  # here plant uptake as absolute parameter
  dI <-  +parms$iI  - parms$kIP  - leach  + PhiU  + PhiB  + PhiTvr
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
  #
  resDeriv <- structure(as.numeric(
    c( dR, dRN, dL, dLN, dI, dAlpha))
    ,names = c("dR","dRN","dL","dLN","dI","dAlpha"))
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
    c(uN = uNSubstrate, usage = synE/parms$cnE + synB/parms$cnB + PhiB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dR + dL + tvrExC + resp,    parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + dRN + dLN + dI + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL - plantNUp  + parms$iI -
         parms$kIP - parms$l*x["I"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) { resDeriv["dR"] <- resDeriv["dRN"] <-  0   }
  if (isTRUE(parms$isFixedL)) { resDeriv["dL"] <- resDeriv["dLN"] <-  0   }
  if (isTRUE(parms$isFixedI)) { resDeriv["dI"] <-  0   }
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
  PhiBU <- PhiB + PhiU
  # total mineralization flux including microbial turnover
  PhiTotal <- PhiBU + PhiTvr
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
    , B = as.numeric(B), dB = as.numeric(dB)
    , ER = as.numeric(ER), EL = as.numeric(EL)
    #, MmB = as.numeric(MmB)
    , PhiB = as.numeric(PhiB), PhiU = as.numeric(PhiU)
    , PhiTvr = as.numeric(PhiTvr)
    , PhiBU = as.numeric(PhiBU), PhiTotal = as.numeric(PhiTotal)
    , immoPot = as.numeric(immoPot)
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
    , dI = as.numeric(dI)
    , uC = as.numeric(uC), synB = as.numeric(synB)
    , decN = as.numeric(decN)
  ))
}

sesam3BSteady <- function(
  ### compute quasi steady sate of microbial biomass given C anc N fluxes
  dLC        ##<< potential decomposition flux from L decomposition in C units
  , dRC      ##<< potential decomposition flux from L decomposition in C units
  , dLN      ##<< potential decomposition flux from L decomposition in N units
  , dRN      ##<< potential decomposition flux from L decomposition in N units
  , alpha    ##<< microbial community partitioning
  , parms    ##<< list of parmaters (aE, eps, km, kN, m, tau, cnB)
  , immNPot  ##<< maximum immobilization rate in N units
){
  BC <- sesam3BSteadyClim(dLC, dRC, alpha, parms)
  BN <- sesam3BSteadyNlim(dLN, dRN, alpha, parms, immNPot = immNPot)
  min(BC, BN)
}

.depr.sesam3BSteadyClim <- function(
  ### compute quasi steady sate of microbial biomass given
  dL
  , dR
  , alpha
  , parms
){
  aE = parms[["aE"]]
  eps = parms[["eps"]]
  kmkN = parms[["km"]] * parms[["kN"]]
  m = parms[["m"]]
  tau = parms[["tau"]]
  kmkN2 = kmkN*kmkN
  Bsteady <-  aE*alpha^2*dL*eps + aE*alpha^2*dR*eps - aE*alpha*dL*eps -
    (aE*alpha*dR*eps + eps*kmkN*m + kmkN*tau -
    sqrt(aE^2*alpha^4*dL^2*eps^2 + 2*aE^2*alpha^4*dL*dR*eps^2 +
           aE^2*alpha^4*dR^2*eps^2 - 2*aE^2*alpha^3*dL^2*eps^2 -
           4*aE^2*alpha^3*dL*dR*eps^2 - 2*aE^2*alpha^3*dR^2*eps^2 +
           aE^2*alpha^2*dL^2*eps^2 + 2*aE^2*alpha^2*dL*dR*eps^2 +
           aE^2*alpha^2*dR^2*eps^2 + 4*aE*alpha^3*dL*eps^2*kmkN*m +
           4*aE*alpha^3*dL*eps*kmkN*tau - 4*aE*alpha^3*dR*eps^2*kmkN*m -
           4*aE*alpha^3*dR*eps*kmkN*tau - 6*aE*alpha^2*dL*eps^2*kmkN*m -
           6*aE*alpha^2*dL*eps*kmkN*tau + 6*aE*alpha^2*dR*eps^2*kmkN*m +
           6*aE*alpha^2*dR*eps*kmkN*tau + 2*aE*alpha*dL*eps^2*kmkN*m +
           2*aE*alpha*dL*eps*kmkN*tau - 2*aE*alpha*dR*eps^2*kmkN*m -
           2*aE*alpha*dR*eps*kmkN*tau + 4*alpha^2*eps^2*kmkN2*m^2 +
           8*alpha^2*eps*kmkN2*m*tau + 4*alpha^2*kmkN2*tau^2 -
           4*alpha*eps^2*kmkN2*m^2 - 8*alpha*eps*kmkN2*m*tau -
           4*alpha*kmkN2*tau^2 + eps^2*kmkN2*m^2 +
           2*eps*kmkN2*m*tau + kmkN2*tau^2)) /
  (2*aE*alpha*(alpha*eps*m + alpha*tau - eps*m - tau))
}
sesam3BSteadyClim <- function(
  ### compute quasi steady sate of microbial biomass given
  dL
  , dR
  , alpha
  , parms
){
  aE = parms[["aE"]]
  eps = parms[["eps"]]
  kmkN = parms[["km"]] * parms[["kN"]]
  m = parms[["m"]]
  tau = parms[["tau"]]
  tem <- (tau/eps+m)
  a <- -tem*alpha*(1 - alpha)*aE*aE
  b <- aE*aE*alpha*(1 - alpha)*(dL + dR) - tem*kmkN*aE
  c <- kmkN*aE*((1-alpha)*dL + alpha*dR) - tem*kmkN*kmkN
  B <- max(c(0,solveSquare(a,b,c)))
  #solveSquare(a,b,c)
}
.tmp.f <- function(){
  #breakpoint after computing B root
  decLN <- dLN*(1 - alpha)*aE*B/(kmkN + (1 - alpha)*aE*B)
  decRN <- dRN*(alpha)*aE*B/(kmkN + (alpha)*aE*B)
  tvrBN <- tauN*B
  decLN + decRN + immNPot - tvrBN
  c1 <- (kmkN+alpha*aE*B)*(kmkN+(1-alpha)*aE*B)
  immNPot*c1
  decLN/c1 + decRN/c1
  .tmp.f <- function(){
    c1b <- B^2*alpha*(1-alpha)*aE2 + B*kmkN*aE + kmkN2
    c(c1,c1b) # correct
  }
  c(decLN, decRN, immNPot, decLN + decRN + immNPot, tvrBN)
  c(B^3*a + B^2*b + B*c + d, (decLN + decRN + immNPot - tvrBN)*c1)

}

sesam3BSteadyNlim <- function(
  ### compute quasi steady sate of microbial biomass given
  dLN        ##<< potential decomposition flux from L decomposition in N units
  , dRN      ##<< potential decomposition flux from L decomposition in N units
  , alpha
  , parms
  , immNPot
){
  betaB <- parms[["cnB"]]
  aE = parms[["aE"]]
  aE2 <- aE*aE
  kmkN = parms[["km"]] * parms[["kN"]]
  tau = parms[["tau"]]
  nu = parms[["nu"]]  # N efficiency during DON uptake
  tauN = tau/betaB/nu
  immoNu = immNPot/nu
  kmkN2 = kmkN*kmkN
  a <- -tauN * alpha*(1 - alpha)*aE2
  b <- aE2*alpha*(1 - alpha)*(dLN + dRN + immoNu) - tauN*kmkN*aE
  c <- aE*kmkN*((1 - alpha)*dLN + alpha*dRN + immoNu) - tauN*kmkN2
  d <- kmkN2*immoNu
  # see test_modSesam3B where the three roots are plotted.
  # Only the second one is positive and increases with increasing immNPot
  B <- solveCubic(a,b,c,d)
  ### numeric scalar: biomass to be in steady state with N balance
  Re(B[2])
}
.tmp.f <- function(){
  #B <- Re(B[2])
  #breakpoint after computing B root
  decLN <- dLN*(1 - alpha)*aE*B/(kmkN + (1 - alpha)*aE*B)
  decRN <- dRN*(alpha)*aE*B/(kmkN + (alpha)*aE*B)
  tvrBN <- tauN*B
  decLN + decRN + immNPot - tvrBN
  c1 <- (kmkN+alpha*aE*B)*(kmkN+(1-alpha)*aE*B)
  immNPot*c1
  decLN/c1 + decRN/c1
  .tmp.f <- function(){
    c1b <- B^2*alpha*(1-alpha)*aE2 + B*kmkN*aE + kmkN2
    c(c1,c1b) # correct
  }
  c(decLN, decRN, immNPot, decLN + decRN + immNPot, tvrBN)
  c(B^3*a + B^2*b + B*c + d, (decLN + decRN + immNPot - tvrBN)*c1)

}

solveSquare <- function(
  ### provides the solutions of eq. 0 = a x^2 + b x + c
  a,b,c   ##<< coefficients of the square equation
){
  p2 <- b/a/2
  q <- c/a
  D <- sqrt(p2*p2 -q)
  x0 <- -p2 +c(+1,-1)*D
  ### complex vector of length 2: roots
  x0
}
attr(solveSquare,"ex") <- function(){
  x1 <- -1
  x2 <- 2
  a <- 2
  b <- -(x1+x2)*a
  c <- x1*x2*a
  #c(a,b,c)
  x <- seq(x1 - 0.2, x2 + 0.2, length.out = 31)
  y2 <- (x - x1)*(x - x2)
  plot( y2 ~ x, type = "l")
  abline(h = 0, col = "grey")
  abline(v = c(x1,x2), col = "grey", lty = "dotted")
  (x0 <- solveSquare(a,b,c))
  rbind(c(x1,x2), sort(x0))
}

solveCubic <- function(
  ### provides the solutions of eq. 0 = a x^3 + b x^2 + c x + d
  a,b,c,d   ##<< coefficients of the cubic equation
){
  # http://www.mathe.tu-freiberg.de/~hebisch/seminar1/kubik.html
  p <- 3*a*c - b*b
  q <- 2*b*b*b - 9*a*b*c + 27*a*a*d
  sqD <- 4*as.complex(q*q + 4*p*p*p)^(1/2)
  u <- as.complex(-4*q + sqD)^(1/3)/2
  v <- as.complex(-4*q - sqD)^(1/3)/2
  y1 <- u + v
  D2 <- sqrt(3)*1i/2
  y2 <- -(u + v)/2 + (u - v)*D2
  y3 <- -(u + v)/2 - (u - v)*D2
  x0 <- (c(y1,y2,y3) - b)/3/a
  ### complex vector of length 3: roots
  x0
}
attr(solveCubic,"ex") <- function(){
  x1 <- -1
  x2 <- 1.8
  x3 <-  2
  a <- 1
  b <- -x1 - (x2 + x3)
  c <- x2*x3 + x1*(x2 + x3)
  d <- -x1*x2*x3
  x <- seq(x1 - 0.2, x3 + 0.2, length.out = 31)
  y <- a*x^3 + b*x^2 + c*x + d
  y2 <- (x - x1)*(x - x2)*(x - x3)
  plot( y2 ~ x, type = "l")
  lines( y ~ x, col = "blue")
  abline(h = 0, col = "grey")
  abline(v = c(x1,x2,x3), col = "grey", lty = "dotted")
  (x0 <- solveCubic(a,b,c,d))
  rbind(c(x1,x2,x3), sort(Re(x0)))
}





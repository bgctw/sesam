#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam3B <- function(
  ### Soil Enzyme Steady Allocation model with biomass in quasi steady state
  t,x,parms
){
  ##details<<
  ## Simplified version of sesam3s
  ## model corresponding to Seam3 with enzyme levels computed by quasi steady state
  ## Alpha as an explicit state variable that changes with turnover
  ## Simplified computation of target based on current revenue based on current alpha
  x <- pmax(unlist(x),1e-16)      # no negative masses
  # compute steady state enzyme levels for N and for C limitation
  dRPot <- parms$kR * x["R"]
  dLPot <- parms$kL * x["L"]
  immoPot <- parms$iB * x["I"]
  cnR <- x["R"]/x["RN"]
  cnL <- x["L"]/x["LN"]
  cnE <- parms$cnE
  cnB <- parms$cnB
  alpha <- x["alpha"]
  # steady state biomass
  B <- sesam3BSteady(dLC = dLPot, dRC = dRPot
                     , dLN = dLPot/cnL, dRN = dRPot/cnR
                     , alpha = alpha, parms = parms
                     , immNPot = immoPot
  )
  aeB <- synERev <- parms$aE*B        # aeB without associanted growth respiration
  kmN <- parms$kmN  #parms$km*parms$kN
  rM <- parms$m*B          # maintenance respiration
  tvrB <- parms$tau*B      # microbial turnover
  synE <- if (isTRUE(parms$isEnzymeMassFlux)) aeB else 0
  #
  # enzyme limitations of decomposition
  decL <- dLPot * (1 - alpha)*aeB/(kmN + (1 - alpha)*aeB)
  decR <- dRPot * alpha*aeB/(kmN + alpha*aeB)
  #
  tvrERecycling <- parms$kNB*synE
  uNOrg <- parms$nu*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  uC <- decL + decR + tvrERecycling
  CsynBC <- uC - rM - synE/parms$eps
  CsynBN <- cnB/parms$eps*(uNOrg + immoPot - synE/cnE)
  CsynB <- min(CsynBC, CsynBN)
  if (CsynB > 0) {
    synB <- parms$eps*CsynB
    rG <- (1 - parms$eps)*CsynB
  } else {
    synB <- CsynB # with negative biomass change, do growth respiration
    rG <- 0
  }
  PhiB <- uNOrg - synB/cnB - synE/cnE
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
  alphaTarget <- balanceAlphaBetweenCNLimitationsExp(
    alphaC, alphaN, CsynBN, CsynBC, tauB = parms$tau*B  )
  # microbial community change as fast as microbial turnover
  dAlpha <- (alphaTarget - alpha) *  (parms$tau + abs(synB)/B)
  #
  # imbalance fluxes of microbes and predators (consuming part of microbial turnover)
  respO <- uC - (synE/parms$eps + synB + rG + rM)
  respTvr <- (1 - parms$epsTvr) * tvrB
  # assuming same cnRatio of predators to equal cn ratio of microbes
  PhiTvr <- respTvr/parms$cnB
  #
  # tvr that feeds R pool, assume that N in SOM for resp (by epsTvr) is mineralized
  tvrC <-  +parms$epsTvr*tvrB   + (1 - parms$kNB)*synE
  tvrN <-  +parms$epsTvr*tvrB/parms$cnB   + (1 - parms$kNB)*synE/parms$cnE
  # fluxes leaving the system (will be set in scen where trv does not feed back)
  tvrExC <- tvrExN <- 0
  #
  leach <- parms$l*x["I"]
  plantNUpPot <- parms$kIPlant*x["I"]
  plantNUp <- if (!is.null(parms$plantNUpAbs)) {
    min(parms$plantNUpAbs, plantNUpPot)
  } else {
    plantNUpPot
  }
  PhiU <- (1 - parms$nu)*(decL/cnL + decR/cnR + tvrERecycling/cnE)
  #
  dB <- synB - tvrB
  dL <- -decL  + parms$iL
  dLN <- -decL/cnL   + parms$iL/parms$cnIL
  dR <- -decR  + parms$iR  + tvrC
  dRN <- -decR/cnR  + parms$iR/parms$cnIR  + tvrN
  # here plant uptake as absolute parameter
  dI <-  +parms$iI  - plantNUp  - leach  + PhiU  + PhiB  + PhiTvr
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
  respB <- (synE)/parms$eps*(1 - parms$eps)  + rG + rM + respO
  resp <- respB + respTvr
  if (diff( unlist(c(uC = uC, usage = respB + synB + synE )))^2 > sqrEps )  stop(
    "biomass mass balance C error")
  if (diff( unlist(
    c(uN = uNOrg, usage = synE/parms$cnE + synB/parms$cnB + PhiB )))^2 >
    .Machine$double.eps)  stop("biomass mass balance N error")
  if (!isTRUE(parms$isFixedS)) {
    if (diff(unlist(
      c(dB + dR + dL + tvrExC + resp,    parms$iR + parms$iL )))^2 >
      sqrEps )  stop("mass balance C error")
    if (diff(unlist(
      c( dB/parms$cnB  + dRN + dLN + dI + tvrExN
         , parms$iR/parms$cnIR  + parms$iL/parms$cnIL + parms$iI -
         plantNUp - parms$l*x["I"])))^2 >
      .Machine$double.eps )  stop("mass balance dN error")
  }
  #
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) { resDeriv["dR"] <- resDeriv["dRN"] <-  0   }
  if (isTRUE(parms$isFixedL)) { resDeriv["dL"] <- resDeriv["dLN"] <-  0   }
  if (isTRUE(parms$isFixedI)) { resDeriv["dI"] <-  0   }
  #
  # further computations just for output for tacking the system
  if (!is.null(parms$kN)) {
    ER <- alpha * parms$aE * B / parms$kN
    EL <- (1 - alpha) * parms$aE * B / parms$kN
  } else {
    kN = 60       ##<< /yr enzyme turnover 60 times a year, each 6 days
    # for output of enzyme levels only
    ER <- alpha * parms$aE * B / kN
    EL <- (1 - alpha) * parms$aE * B / kN
  }
  limER <- alpha*synERev/(parms$kmN + alpha*synERev)
  limEL <- (1 - alpha)*synERev/(parms$kmN + (1 - alpha)*synERev)
  revRC <- dRPot / (parms$kmN + alphaC*aeB)
  revLC <- dLPot / (parms$kmN + (1 - alphaC)*aeB)
  revRN <- dRPot/cnR / (parms$kmN + alphaN*aeB)
  revLN <- dLPot/cnL / (parms$kmN + alphaN*aeB)
  # net mic mineralization/immobilization when accounting uptake mineralization
  PhiBU <- PhiB + PhiU
  # total mineralization flux including microbial turnover
  PhiTotal <- PhiBU + PhiTvr
  # do not match in other limitation
  # c(alphaC, revRC/(revRC + revLC)); c(alphaN, revRN/(revRN + revLN))
  # compute C available for biomass, this time without accounting immobilization flux
  NsynBNSubstrate <- uNOrg - synE/cnE
  CsynBNSubstrate <- (NsynBNSubstrate*cnB)/parms$eps
  # N limited based on substrate uptake (without accounting for immobilization)
  isLimNSubstrate <-  CsynBNSubstrate < CsynBC
  #
  if (isTRUE(parms$isRecover) ) recover()
  list( resDeriv, c(
     B = as.numeric(B), dB = as.numeric(dB)
     , resp = as.numeric(resp)
     , respO = as.numeric(respO)
     , respB = as.numeric(respB)
     , respTvr = as.numeric(respTvr)
     #, ER = as.numeric(ER), EL = as.numeric(EL)
    #, MmB = as.numeric(MmB)
    , PhiTotal = as.numeric(PhiTotal)
    , PhiB = as.numeric(PhiB), PhiU = as.numeric(PhiU)
    , PhiTvr = as.numeric(PhiTvr)
    , PhiBU = as.numeric(PhiBU)
    , immoPot = as.numeric(immoPot)
    , alphaTarget = as.numeric(alphaTarget)
    , alphaC = as.numeric(alphaC), alphaN = as.numeric(alphaN)
    , cnR = as.numeric(cnR), cnL = as.numeric(cnL)
    , limER = as.numeric(limER), limEL = as.numeric(limEL)
    , decR = as.numeric(decR), decL = as.numeric(decL)
    , tvrB = as.numeric(tvrB)
    , revRC = as.numeric(revRC), revLC = as.numeric(revLC)
    , revRN = as.numeric(revRN)
    , revLN = as.numeric(revLN)
    , CsynB = as.numeric(CsynB)
    , CsynBC = as.numeric(CsynBC)
    , CsynBN = as.numeric(CsynBN)
    #, pNsyn = as.numeric(NsynBN / (parms$eps*CsynBC/cnB) )
    #, NsynReq = as.numeric(CsynBC/cnB), Nsyn = as.numeric(NsynBN)
    #, dR = as.numeric(dR), dL = as.numeric(dL), dB = as.numeric(dB)
    #, dI = as.numeric(dI)
    #, uC = as.numeric(uC), synB = as.numeric(synB)
    #, decN = as.numeric(decN)
    , plantNUp = plantNUp
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

sesam3BSteadyClim <- function(
  ### compute quasi steady sate of microbial biomass given
  dL
  , dR
  , alpha
  , parms
){
  aE = parms[["aE"]]
  eps = parms[["eps"]]
  kmkN = parms$kmN #parms[["km"]] * parms[["kN"]]
  m = parms[["m"]]
  tau = parms[["tau"]]
  kappaE = parms[["kNB"]]
  #tem <- tau/eps + m  # for disregarding enzyme mass fluxes
  tem <- tau/eps + m + (1/eps - kappaE)*aE
  a <- -tem*alpha*(1 - alpha)*aE*aE
  b <- aE*aE*alpha*(1 - alpha)*(dL + dR) - tem*kmkN*aE
  c <- kmkN*aE*((1 - alpha)*dL + alpha*dR) - tem*kmkN*kmkN
  #B <- max(c(0,solveSquare(a,b,c)))
  B <- max(c(1e-16,solveSquare(a,b,c))) # sustain minimal biomass
  #solveSquare(a,b,c)
}
.tmp.f <- function(){
  #B <- Re(B[2])
  #breakpoint after computing B root
  decL <- dL*(1 - alpha)*aE*B/(kmkN + (1 - alpha)*aE*B)
  decR <- dR*(alpha)*aE*B/(kmkN + (alpha)*aE*B)
  decE <- kappaE*aE*B
  tvrB <- tau*B
  maint <- m*B
  synE <- aE*B/eps
  c1 <- (kmkN + alpha*aE*B)*(kmkN + (1 - alpha)*aE*B)
  .tmp.f <- function(){
    c1b <- B^2*alpha*(1 - alpha)*aE*aE + B*kmkN*aE + kmkN*kmkN
    c(c1,c1b) # correct
  }
  c(decL, decR, decE = decE, maint = maint
    , synB = eps*(decL + decR + decE - synE - maint), tvrB = tvrB)
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
  kmkN = parms$kmN #parms[["km"]] * parms[["kN"]]
  tau = parms[["tau"]]
  nu = parms[["nu"]]  # N efficiency during DON uptake
  kappaE = parms[["kNB"]]
  betaE = parms[["cnE"]]
  #tauN = tau/betaB/nu # for neglecting enzyme mass fluxes
  tauN = tau/betaB/nu + (1/nu -  kappaE)*aE/betaE
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
  synEN <- aE*B/betaB
  decEN <- kappaE*synE
  tvrBN <- tau*B/betaB
  c1 <- (kmkN + alpha*aE*B)*(kmkN + (1 - alpha)*aE*B)
  .tmp.f <- function(){
    c1b <- B^2*alpha*(1 - alpha)*aE2 + B*kmkN*aE + kmkN2
    c(c1,c1b) # correct
  }
  c(decLN, decRN, decEN, immNPot,
    nu*(decLN + decRN + decEN) + immNPot - synEN, tvrBN)
  c(B^3*a + B^2*b + B*c + d, (decLN + decRN + immNPot - tvrBN)*c1)
}

solveSquare <- function(
  ### provides the solutions of eq. 0 = a x^2 + b x + c
  a,b,c   ##<< coefficients of the square equation
){
  p2 <- b/a/2
  q <- c/a
  D <- sqrt(p2*p2 - q)
  x0 <- -p2 + c(+1,-1)*D
  ### complex vector of length 2: roots
  x0
}
attr(solveSquare,"ex") <- function(){
  x1 <- -1
  x2 <- 2
  a <- 2
  b <- -(x1 + x2)*a
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





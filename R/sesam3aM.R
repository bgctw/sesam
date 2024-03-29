#library(deSolve)

# gC/m2 and gN/m2, /yr

derivSesam3aM <- function(
  ### Soil Enzyme Steady Allocation model
  t,x,parms
){
  ##details<<
  ## Matrix formulation variant of sesam3a
  ## see inst/docu_sesam2/Sesam3MatrixForm.Rmd
  x0 <- pmax(unlist(x),1e-16)      # no negative masses
  alpha <- x["alpha"]
  B <- x["B"]
  cnR <- x["R"]/x["RN"]
  cnL <- x["L"]/x["LN"]
  cnE <- parms$cnE
  cnB <- parms$cnB
  x <- c(x0, BN = setNames(B/cnB,NULL))  # need BN explicitly
  # transfer matrix
  T <- diag(-1, nrow = length(x)); colnames(T) <- rownames(T) <- names(x)
  kappaE <- parms$kNB
  #
  #-- carbon fluxes
  lL <- (1 - alpha)*parms$aE / (parms$kmN + (1 - alpha)*parms$aE*B)
  lR <- alpha*parms$aE / (parms$kmN + alpha*parms$aE*B)
  decL <- parms$kL*lL # not multiplied by B and L yet
  decR <- parms$kR*lR
  outL <- decL*B
  outR <- decR*B
  #
  # difference of enzmye production cost and enzTvr uptake per B
  mE <- parms$aE*(1/parms$eps - kappaE)
  CsynBCt <- decL*x["L"] + decR*x["R"] - mE - parms$m
  CsynBC <- ifelse(CsynBCt > 0, parms$eps*CsynBCt, CsynBCt)
  #
  #-- N fluxes
  outLN <- outL
  outRN <- outR
  # plant N uptake can have an upper absolute limit, compute rate per IN
  outI <- min(parms$kINPlant, parms$plantNUpAbs/x["IN"]) + parms$lN + parms$iBN
  #
  uN <- parms$nuN*(decL*x["LN"] + decR*x["RN"] + kappaE*parms$aE/cnE)*cnB +
    parms$iBN*x["IN"]/(B/cnB)
  NsynBN <- uN - parms$aE*cnB/cnE
  CsynBN <- NsynBN  # *B_N*cnB/BC: same proportion, see Sesam3MatrixForm.Rmd
  #
  #-- depending on synB both C anc N
  synB <- min(CsynBC, CsynBN)
  #
  rG <- if (synB > 0) (1 - parms$eps)/parms$eps*synB else 0
  rO <- CsynBCt - (synB + rG)
  outB <- parms$tau + parms$m + mE + rG + rO
  T["L","B"] <- T["R","B"] <- 1
  T["B","R"] <- (parms$epsTvr*parms$tau + (1 - kappaE)*parms$aE) / outB
  #
  MImb <- NsynBN - synB # *B/beta_B/B_N: see Sesam3MatrixForm.Rmd
  tvrEN <- parms$aE*cnB/cnE
  outBNEnzR <- tvrEN*(1 - kappaE)
  outBNEnzI <- (1 - parms$nuN)*tvrEN*kappaE
  outBN <-  parms$tau + outBNEnzR + outBNEnzI + MImb
  T["LN","BN"] <- T["RN","BN"] <- parms$nuN
  T["LN","IN"] <- T["RN","IN"] <- 1 - parms$nuN
  T["BN","RN"] <- (parms$epsTvr*parms$tau + outBNEnzR)/outBN
  T["BN","IN"] <- ((1 - parms$epsTvr)*parms$tau + MImb + outBNEnzI)/outBN
  T["IN","BN"] <- parms$iBN / outI
  #
  if (isTRUE(parms$isTvrNil)) {
    # scenario of enzymes and biomass not feeding back to R
    T["B","R"] <- T["BN","RN"] <- 0
  }
  #
  # I,N here refer to inputs and transfer matrices of Sierra
  I <- c(B = 0, R = parms$iR, RN = parms$iR/parms$cnIR
         , L = parms$iL, LN = parms$iL/parms$cnIL
         , IN = parms$iIN, alpha = 0, BN = 0)
  N <- setNames(c(outB, outR, outRN, outL, outLN, outI, 0, outBN), names(x))
  loss <- N*x
  dx <- setNames(as.vector(t(T)%*%(N*as.matrix(x, ncol = 1))), names(x)) + I
  #
  #--- community composition
  dRPot <- parms$kR * x["R"]
  dLPot <- parms$kL * x["L"]
  alphaC <- computeSesam3sAllocationPartitioning(
    dR = dRPot, dL = dLPot, B = B
    , kmkN = parms$kmN, aE = parms$aE
    , alpha = alpha
  )
  alphaN <- computeSesam3sAllocationPartitioning(
    dR = dRPot/cnR, dL = dLPot/cnL, B = B
    , kmkN = parms$kmN, aE = parms$aE
    , alpha = alpha
  )
  # not that CsynBE must be multiplied with biomass B here
  alphaTarget <- balanceAlphaBetweenCNLimitationsExp(
    alphaC, alphaN, CsynBN*B, CsynBC*B, tauB = parms$tau*B  )
  # microbial community change as fast as microbial turnover
  dAlpha <- dx["alpha"] <- (alphaTarget - alpha) *  (parms$tau + abs(synB))
  #
  if (isTRUE(parms$isFixedS)) {
    # scenario of fixed substrate
    dx["R"] <- dx["L"] <- dx["RN"] <- dx["LN"] <- dx["IN"] <- 0
  }
  #
  if (any(!is.finite(dx))) stop("encountered nonFinite derivatives")
  sqrEps <- sqrt(.Machine$double.eps)
  #
  if (dx["B"] - cnB*dx["BN"] > sqrEps) stop(
    "error in computing microbial biomass change stoichiometry")
  sT <- rowSums(T)
  if (any(abs(sT[c("R","RN","L","LN","BN")]) > sqrEps)) stop(
    "mass balance error of conserved pools")
  # allowing scenarios with holding some pools fixed
  if (isTRUE(parms$isFixedR)) { dx["R"] <- dx["RN"] <-  0   }
  if (isTRUE(parms$isFixedL)) { dx["L"] <- dx["LN"] <-  0   }
  if (isTRUE(parms$isFixedI)) { dx["IN"] <-  0   }
  #
  if (isTRUE(parms$isRecover) ) recover()
  # for matrix formulation, explicitly tracked BN, but do not return here.
  list( dx[1:length(x0)], c(
  ))
}

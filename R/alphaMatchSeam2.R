calcMatchAlphaSeam2 <- function(
  ### optimal enzyme allocation for given optimal biomass ratio and enzyme levels
  E           ##<< numeric scalar: sum of current enzyme amounts
  , decPotR   ##<< numeric scalar: current potential decomposition of residue pool
    ## (with out accounting for enzyme limitations)
  , decPotL   ##<< numeric scalar: current potential decomposition of litter pool
    ## (with out accounting for enzyme limitations)
  , rMaint    ##<< numeric scalar: maintenance respiration
  , parms     ##<< list of parameter values including substrate affinities k
    ## mR, kmL, and anabolic efficiency eps
  , cnOpt = parms$cnB ##<< numeric scalar: target C/N ratio of microbial biomass
  , cnR = parms$cnIR    ##<< numeric scalar: current C/N ratio of residue pool
  , cnL = parms$cnIL    ##<< numeric scalar: current C/N ratio of litter pool
  , imm = 0   ##<< immobilization flux
){
  ##details<< Differs from calcMatchAlphaSeam by accounting for immobilization flux \code{imm}
  kmR <- parms$kmR
  kmL <- parms$kmL
  eps <- parms$eps
  # simplified and faster calculation for kmL==kmR, else uncomment the second line below
  if (parms$kmL != parms$kmR ) stop("calcMatchAlphaSeam2: only valid for parameters kmL==kmR ")
  km <- parms$kmL
  alpha <-
    # with kmL == kmR
    -(-E^2*cnL*cnOpt*cnR*imm - E^2*cnL*cnOpt*decPotR + E^2*cnL*cnR*decPotL*eps +
        E^2*cnL*cnR*decPotR*eps - E^2*cnL*cnR*eps*rMaint - E^2*cnOpt*cnR*decPotL -
        E*cnL*cnOpt*decPotR*kmR - E*cnL*cnR*decPotL*eps*kmR + E*cnL*cnR*decPotR*eps*kmR +
        E*cnOpt*cnR*decPotL*kmR +
        sqrt(
          E^2*cnL^2*cnOpt^2*cnR^2*imm^2 + 2*E^2*cnL^2*cnOpt^2*cnR*decPotR*imm +
            E^2*cnL^2*cnOpt^2*decPotR^2 - 2*E^2*cnL^2*cnOpt*cnR^2*decPotL*eps*imm -
            2*E^2*cnL^2*cnOpt*cnR^2*decPotR*eps*imm + 2*E^2*cnL^2*cnOpt*cnR^2*eps*imm*rMaint -
            2*E^2*cnL^2*cnOpt*cnR*decPotL*decPotR*eps - 2*E^2*cnL^2*cnOpt*cnR*decPotR^2*eps +
            2*E^2*cnL^2*cnOpt*cnR*decPotR*eps*rMaint + E^2*cnL^2*cnR^2*decPotL^2*eps^2 +
            2*E^2*cnL^2*cnR^2*decPotL*decPotR*eps^2 - 2*E^2*cnL^2*cnR^2*decPotL*eps^2*rMaint +
            E^2*cnL^2*cnR^2*decPotR^2*eps^2 - 2*E^2*cnL^2*cnR^2*decPotR*eps^2*rMaint +
            E^2*cnL^2*cnR^2*eps^2*rMaint^2 + 2*E^2*cnL*cnOpt^2*cnR^2*decPotL*imm +
            2*E^2*cnL*cnOpt^2*cnR*decPotL*decPotR - 2*E^2*cnL*cnOpt*cnR^2*decPotL^2*eps -
            2*E^2*cnL*cnOpt*cnR^2*decPotL*decPotR*eps +
            2*E^2*cnL*cnOpt*cnR^2*decPotL*eps*rMaint +
            E^2*cnOpt^2*cnR^2*decPotL^2 + 4*E*cnL^2*cnOpt^2*cnR^2*imm^2*kmR +
            6*E*cnL^2*cnOpt^2*cnR*decPotR*imm*kmR + 2*E*cnL^2*cnOpt^2*decPotR^2*kmR -
            6*E*cnL^2*cnOpt*cnR^2*decPotL*eps*imm*kmR -
            6*E*cnL^2*cnOpt*cnR^2*decPotR*eps*imm*kmR +
            8*E*cnL^2*cnOpt*cnR^2*eps*imm*kmR*rMaint -
            4*E*cnL^2*cnOpt*cnR*decPotL*decPotR*eps*kmR -
            4*E*cnL^2*cnOpt*cnR*decPotR^2*eps*kmR +
            6*E*cnL^2*cnOpt*cnR*decPotR*eps*kmR*rMaint +
            2*E*cnL^2*cnR^2*decPotL^2*eps^2*kmR +
            4*E*cnL^2*cnR^2*decPotL*decPotR*eps^2*kmR -
            6*E*cnL^2*cnR^2*decPotL*eps^2*kmR*rMaint +
            2*E*cnL^2*cnR^2*decPotR^2*eps^2*kmR -
            6*E*cnL^2*cnR^2*decPotR*eps^2*kmR*rMaint +
            4*E*cnL^2*cnR^2*eps^2*kmR*rMaint^2 +
            6*E*cnL*cnOpt^2*cnR^2*decPotL*imm*kmR +
            4*E*cnL*cnOpt^2*cnR*decPotL*decPotR*kmR -
            4*E*cnL*cnOpt*cnR^2*decPotL^2*eps*kmR -
            4*E*cnL*cnOpt*cnR^2*decPotL*decPotR*eps*kmR +
            6*E*cnL*cnOpt*cnR^2*decPotL*eps*kmR*rMaint +
            2*E*cnOpt^2*cnR^2*decPotL^2*kmR + 4*cnL^2*cnOpt^2*cnR^2*imm^2*kmR^2 +
            4*cnL^2*cnOpt^2*cnR*decPotR*imm*kmR^2 + cnL^2*cnOpt^2*decPotR^2*kmR^2 -
            4*cnL^2*cnOpt*cnR^2*decPotL*eps*imm*kmR^2 -
            4*cnL^2*cnOpt*cnR^2*decPotR*eps*imm*kmR^2 +
            8*cnL^2*cnOpt*cnR^2*eps*imm*kmR^2*rMaint +
            2*cnL^2*cnOpt*cnR*decPotL*decPotR*eps*kmR^2 -
            2*cnL^2*cnOpt*cnR*decPotR^2*eps*kmR^2 +
            4*cnL^2*cnOpt*cnR*decPotR*eps*kmR^2*rMaint + cnL^2*cnR^2*decPotL^2*eps^2*kmR^2 -
            2*cnL^2*cnR^2*decPotL*decPotR*eps^2*kmR^2 - 4*cnL^2*cnR^2*decPotL*eps^2*kmR^2*rMaint +
            cnL^2*cnR^2*decPotR^2*eps^2*kmR^2 - 4*cnL^2*cnR^2*decPotR*eps^2*kmR^2*rMaint +
            4*cnL^2*cnR^2*eps^2*kmR^2*rMaint^2 + 4*cnL*cnOpt^2*cnR^2*decPotL*imm*kmR^2 -
            2*cnL*cnOpt^2*cnR*decPotL*decPotR*kmR^2 - 2*cnL*cnOpt*cnR^2*decPotL^2*eps*kmR^2 +
            2*cnL*cnOpt*cnR^2*decPotL*decPotR*eps*kmR^2 +
            4*cnL*cnOpt*cnR^2*decPotL*eps*kmR^2*rMaint +
            cnOpt^2*cnR^2*decPotL^2*kmR^2) * abs(E)) /
    (2*E^2*cnL*cnOpt*cnR*imm + 2*E^2*cnL*cnOpt*decPotR - 2*E^2*cnL*cnR*decPotL*eps -
       2*E^2*cnL*cnR*decPotR*eps + 2*E^2*cnL*cnR*eps*rMaint + 2*E^2*cnOpt*cnR*decPotL)
  # allowed kmL != kmR
  #-(-E^2*cnL*cnOpt*cnR*imm - E^2*cnL*cnOpt*decPotR + E^2*cnL*cnR*decPotL*eps +E^2*cnL*cnR*decPotR*eps - E^2*cnL*cnR*eps*rMaint - E^2*cnOpt*cnR*decPotL - E*cnL*cnOpt*cnR*imm*kmL + E*cnL*cnOpt*cnR*imm*kmR - E*cnL*cnOpt*decPotR*kmL - E*cnL*cnR*decPotL*eps*kmR + E*cnL*cnR*decPotR*eps*kmL - E*cnL*cnR*eps*kmL*rMaint + E*cnL*cnR*eps*kmR*rMaint + E*cnOpt*cnR*decPotL*kmR + sqrt(E^2*cnL^2*cnOpt^2*cnR^2*imm^2 + 2*E^2*cnL^2*cnOpt^2*cnR*decPotR*imm + E^2*cnL^2*cnOpt^2*decPotR^2 - 2*E^2*cnL^2*cnOpt*cnR^2*decPotL*eps*imm - 2*E^2*cnL^2*cnOpt*cnR^2*decPotR*eps*imm + 2*E^2*cnL^2*cnOpt*cnR^2*eps*imm*rMaint - 2*E^2*cnL^2*cnOpt*cnR*decPotL*decPotR*eps - 2*E^2*cnL^2*cnOpt*cnR*decPotR^2*eps + 2*E^2*cnL^2*cnOpt*cnR*decPotR*eps*rMaint + E^2*cnL^2*cnR^2*decPotL^2*eps^2 + 2*E^2*cnL^2*cnR^2*decPotL*decPotR*eps^2 - 2*E^2*cnL^2*cnR^2*decPotL*eps^2*rMaint + E^2*cnL^2*cnR^2*decPotR^2*eps^2 - 2*E^2*cnL^2*cnR^2*decPotR*eps^2*rMaint + E^2*cnL^2*cnR^2*eps^2*rMaint^2 + 2*E^2*cnL*cnOpt^2*cnR^2*decPotL*imm + 2*E^2*cnL*cnOpt^2*cnR*decPotL*decPotR - 2*E^2*cnL*cnOpt*cnR^2*decPotL^2*eps - 2*E^2*cnL*cnOpt*cnR^2*decPotL*decPotR*eps + 2*E^2*cnL*cnOpt*cnR^2*decPotL*eps*rMaint + E^2*cnOpt^2*cnR^2*decPotL^2 + 2*E*cnL^2*cnOpt^2*cnR^2*imm^2*kmL + 2*E*cnL^2*cnOpt^2*cnR^2*imm^2*kmR + 4*E*cnL^2*cnOpt^2*cnR*decPotR*imm*kmL + 2*E*cnL^2*cnOpt^2*cnR*decPotR*imm*kmR + 2*E*cnL^2*cnOpt^2*decPotR^2*kmL - 2*E*cnL^2*cnOpt*cnR^2*decPotL*eps*imm*kmL - 4*E*cnL^2*cnOpt*cnR^2*decPotL*eps*imm*kmR - 4*E*cnL^2*cnOpt*cnR^2*decPotR*eps*imm*kmL - 2*E*cnL^2*cnOpt*cnR^2*decPotR*eps*imm*kmR + 4*E*cnL^2*cnOpt*cnR^2*eps*imm*kmL*rMaint + 4*E*cnL^2*cnOpt*cnR^2*eps*imm*kmR*rMaint - 2*E*cnL^2*cnOpt*cnR*decPotL*decPotR*eps*kmL - 2*E*cnL^2*cnOpt*cnR*decPotL*decPotR*eps*kmR - 4*E*cnL^2*cnOpt*cnR*decPotR^2*eps*kmL + 4*E*cnL^2*cnOpt*cnR*decPotR*eps*kmL*rMaint + 2*E*cnL^2*cnOpt*cnR*decPotR*eps*kmR*rMaint + 2*E*cnL^2*cnR^2*decPotL^2*eps^2*kmR + 2*E*cnL^2*cnR^2*decPotL*decPotR*eps^2*kmL + 2*E*cnL^2*cnR^2*decPotL*decPotR*eps^2*kmR - 2*E*cnL^2*cnR^2*decPotL*eps^2*kmL*rMaint - 4*E*cnL^2*cnR^2*decPotL*eps^2*kmR*rMaint + 2*E*cnL^2*cnR^2*decPotR^2*eps^2*kmL - 4*E*cnL^2*cnR^2*decPotR*eps^2*kmL*rMaint - 2*E*cnL^2*cnR^2*decPotR*eps^2*kmR*rMaint + 2*E*cnL^2*cnR^2*eps^2*kmL*rMaint^2 + 2*E*cnL^2*cnR^2*eps^2*kmR*rMaint^2 + 2*E*cnL*cnOpt^2*cnR^2*decPotL*imm*kmL + 4*E*cnL*cnOpt^2*cnR^2*decPotL*imm*kmR + 2*E*cnL*cnOpt^2*cnR*decPotL*decPotR*kmL + 2*E*cnL*cnOpt^2*cnR*decPotL*decPotR*kmR - 4*E*cnL*cnOpt*cnR^2*decPotL^2*eps*kmR - 2*E*cnL*cnOpt*cnR^2*decPotL*decPotR*eps*kmL - 2*E*cnL*cnOpt*cnR^2*decPotL*decPotR*eps*kmR + 2*E*cnL*cnOpt*cnR^2*decPotL*eps*kmL*rMaint + 4*E*cnL*cnOpt*cnR^2*decPotL*eps*kmR*rMaint + 2*E*cnOpt^2*cnR^2*decPotL^2*kmR + cnL^2*cnOpt^2*cnR^2*imm^2*kmL^2 + 2*cnL^2*cnOpt^2*cnR^2*imm^2*kmL*kmR + cnL^2*cnOpt^2*cnR^2*imm^2*kmR^2 + 2*cnL^2*cnOpt^2*cnR*decPotR*imm*kmL^2 + 2*cnL^2*cnOpt^2*cnR*decPotR*imm*kmL*kmR + cnL^2*cnOpt^2*decPotR^2*kmL^2 - 2*cnL^2*cnOpt*cnR^2*decPotL*eps*imm*kmL*kmR - 2*cnL^2*cnOpt*cnR^2*decPotL*eps*imm*kmR^2 - 2*cnL^2*cnOpt*cnR^2*decPotR*eps*imm*kmL^2 - 2*cnL^2*cnOpt*cnR^2*decPotR*eps*imm*kmL*kmR + 2*cnL^2*cnOpt*cnR^2*eps*imm*kmL^2*rMaint + 4*cnL^2*cnOpt*cnR^2*eps*imm*kmL*kmR*rMaint + 2*cnL^2*cnOpt*cnR^2*eps*imm*kmR^2*rMaint + 2*cnL^2*cnOpt*cnR*decPotL*decPotR*eps*kmL*kmR - 2*cnL^2*cnOpt*cnR*decPotR^2*eps*kmL^2 + 2*cnL^2*cnOpt*cnR*decPotR*eps*kmL^2*rMaint + 2*cnL^2*cnOpt*cnR*decPotR*eps*kmL*kmR*rMaint + cnL^2*cnR^2*decPotL^2*eps^2*kmR^2 - 2*cnL^2*cnR^2*decPotL*decPotR*eps^2*kmL*kmR - 2*cnL^2*cnR^2*decPotL*eps^2*kmL*kmR*rMaint - 2*cnL^2*cnR^2*decPotL*eps^2*kmR^2*rMaint + cnL^2*cnR^2*decPotR^2*eps^2*kmL^2 - 2*cnL^2*cnR^2*decPotR*eps^2*kmL^2*rMaint - 2*cnL^2*cnR^2*decPotR*eps^2*kmL*kmR*rMaint + cnL^2*cnR^2*eps^2*kmL^2*rMaint^2 + 2*cnL^2*cnR^2*eps^2*kmL*kmR*rMaint^2 + cnL^2*cnR^2*eps^2*kmR^2*rMaint^2 + 2*cnL*cnOpt^2*cnR^2*decPotL*imm*kmL*kmR + 2*cnL*cnOpt^2*cnR^2*decPotL*imm*kmR^2 - 2*cnL*cnOpt^2*cnR*decPotL*decPotR*kmL*kmR - 2*cnL*cnOpt*cnR^2*decPotL^2*eps*kmR^2 + 2*cnL*cnOpt*cnR^2*decPotL*decPotR*eps*kmL*kmR + 2*cnL*cnOpt*cnR^2*decPotL*eps*kmL*kmR*rMaint + 2*cnL*cnOpt*cnR^2*decPotL*eps*kmR^2*rMaint + cnOpt^2*cnR^2*decPotL^2*kmR^2)*abs(E))/(2*E^2*cnL*cnOpt*cnR*imm + 2*E^2*cnL*cnOpt*decPotR - 2*E^2*cnL*cnR*decPotL*eps - 2*E^2*cnL*cnR*decPotR*eps + 2*E^2*cnL*cnR*eps*rMaint + 2*E^2*cnOpt*cnR*decPotL)
  if (alpha < 0 ) alpha <- 0
  if (alpha > 1 ) alpha <- 1
  alpha
}
attr(calcMatchAlphaSeam2, "ex") <- function(){
  if (FALSE) {
    # parms0 from simSteady2.R
    # gC/m2 and gN/m2, /yr
    parms0 <- list(
      cnB = 11
      ,cnE = 3.1     ##<< Sterner02: Protein (Fig. 2.2.), high N investment (low P)
      ,cnIR = 7      ##<< between micr and enzyme signal
      ,cnIL = 30     ##<< N poor substrate, here no inorganic N supply, need to
        # have low C/N for CO2 increase scenario
      ,kN = 60       ##<< /yr enzyme turnover 60 times a year, each 6 days ->
        # fast priming pulses
      ,km = 0.05     ##<< enzyme half-saturation constant, in magnitude of
        # enzymes, determined by kN
      ,kNB = 0.8     ##<< amount of recycling enzyme turnover by biomass
        # (added to uptake instead of R)
      ,kR = 1/(10)   ##<< 1/(x years) # to demonstrate changes on short time scale
      ,kL = 5        ##<< 1/(x years) # formerly 1 year
      ,aE = 0.001*365 ##<< C biomass allocated to enzymes 1/day /microbial biomass
      ,m = 0.005*365 ##<< maintenance respiration rate   1/day /microbial biomass,
        # Bogedom Fig. 1
      ,tau = 1/60*365 ##<< biomass turnover rate (12 days)
      ,eps = 0.5     ##<< carbon use efficiency for growth respiration
      ,epsTvr = 0.3  ##<< carbon use efficiency of microbial tvr (part by
        # predators which respire and corresponding amount of N must be mineralized)
      ,iR = 0        ##<< input modelled explicitely
      ,iL = 400      ##<< g/m2 input per year (half NPP)
      ,plantNUp = 0
      ,useFixedAlloc = FALSE    ##<< set to true to use Fixed allocation (alpha = 0.5)
      ,iP = 10.57    # 0.0289652*365  ##<< plant uptake iP I
      ,iB = 25       #0.0110068*365   ##<< immobilization flux iB I
      ,iI = 0        ##<< input of mineral N
      ,l = 0.96      #0.00262647*365       ##<< leaching rate of mineralN l I
      ,kIP = 10.57  #0.0289652*365          ##<< plant uptake iP I
      ,iB = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iB I
      ,iI = 0        ##<< input of mineral N
      ,l = 0   #0.00262647*365       ##<< leaching rate of mineralN l I
      ,nu = 0.9     # microbial N use efficiency
      ,useAlphaCUntilNLimited = TRUE      ##<< do not decrease investment into
        #C enzmyes when NSubstrateLimtited, but only when N-Limited
    )
    parms0 <- within(parms0,{
      kmR <- kmL <- km
      epsR <- epsL <- eps
      cnER <- cnEL <- cnE
      kNR <- kNL <- kN
    })
    parms <- parms0
    E = 16; decPotR = 1000; decPotL = 1000; respMaint = 0.4
    p <- parms <- parms0; cnOpt  <-  p$cnB; cnR = 7; cnL = 70
    imm = 20
    alpha <- alpha2 <- calcMatchAlphaSeam2(
      E = E, decPotR = decPotR, decPotL = decPotL, rMaint = respMaint, parms = p
      , cnOpt = cnOpt, cnR = cnR, cnL = cnL)
    alpha <- alpha2i <- calcMatchAlphaSeam2(
      E = E, decPotR = decPotR, decPotL = decPotL, rMaint = respMaint, parms = p
      , cnOpt = cnOpt, cnR = cnR, cnL = cnL, imm = imm)
    decR <- decPotR*alpha*E/(p$km + alpha*E)
    decL <- decPotL*(1 - alpha)*E/(p$km + (1 - alpha)*E)
    synC <- p$eps*(decR + decL - respMaint)
    synN <- decR/cnR + decL/cnL + imm
    # check that C/N ratio of fluxes for biomass synthesis equal the target
    c(cnOpt, synC/synN)
  }
}




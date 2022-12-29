#require(testthat)
context("modSesam4b")

parms0 <- list(
  cnB = 7.16
  , cnBW = 10    ##<< C/N ratio of cell walls (that go to R, )
  ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
  #,cnE = 7.16
  ,cnIR = 4.5     ##<< between micr and enzyme signal
  ,cnIL = 30      ##<< N poor substrate
  #,kN = 0.05     ##<< (per day) enzyme turnover
  ##,kN = 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
  , kmN = 0.3*0.01*365  ##<< /yr enzyme turnover 1% turning over each day
  ,kNB = 0.8      ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of R)
  #,kR = 0.2      ##<< substrate decomposition rate N-rich (here simulate large N stock)
  #,kL = 1        ##<< substrate decomposition rate N-poor
  #,kR = 5e-3     ##<< substrate decomposition rate N-rich (here simulate large N stock)
  #,kL = 10e-3    ##<< substrate decomposition rate N-poor
  #,kL = 5e-2     ##<< substrate decomposition rate N-poor
  #,aE = 0.05     ##<< C-uptake allocated to enzymes
  #,kR = 1/(5*365)      ##<< 5 years
  #,kL = 1/(0.5*365)    ##<< 1/2 years
  ,kR = 1/(50)    ##<< 1/(x years)
  ,kL = 1/(1)     ##<< 1/(x years)
  #,aE = 0.003*365 ##<< C biomass allocated to enzymes gC/day /microbial biomass
  ,aE = 0.001*365 ##<< C biomass allocated to enzymes gC/day /microbial biomass
  ##,km = 0.3       ##<< enzyme half-saturation constant
  #,km = 0.03     ##<< enzyme half-saturation constant
  #,km = 14       ##<< enzyme half-saturation constant
  ,m = 0.02*365   ##<< maintenance respiration rate   gC/day /microbial biomass
  ,tau = 1/60*365 ##<< biomass turnover rate (12 days)
  ,eps = 0.5      ##<< carbon use efficiency
  ,epsPred = 0.3  ##<< carbon use efficiency of microbial tvr (predators respire)
  ,iR = 0         ##<< input modelled explicitely
  ,iL = 300       ##<< g/m2 input per year (half NPP)
  #,plantNUp = 300/70*1/4  ##<< plant N uptake balancing N inputs
  ,plantNUp = 0   ##<< plant N uptake balancing N inputs
  ,useFixedAlloc = FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
  ,kINPlant = 10.57 #0.0289652*365         ##<< plant uptake iP IN
  ,iBN = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iBN IN
  ,iIN = 0         ##<< input of mineral N
  ,lN = 0.96 #0.00262647*365       ##<< leaching rate of mineralN lN IN
  , nuN = 0.9       ##<< microbial N use efficiency
  #, isEnzymeMassFlux = FALSE  ##<< steady state B solution neglects enyzme mass fluxes
  , isEnzymeMassFlux = TRUE  ##<< steady state B solution accounts for enyzme mass fluxes
  , nuP = 0.3      ##<< microbial uptake of depolymerized P, (1-nuP) is mineralized
  , cpE = 1e5
  #, cpB = 40
  #, cpBW = 50
  , cpB = 1e8 # for not running into P limitation, here assume no P required in biomass
  , cpBW = 1e8
  , cpIR = 40
  , cpIL = 40*3
  , iBP = 0.38 * 10.57 # start with same as N
  , cW = 0.5 # proportion of cell wall within microbial biomass
  , B0 = 0  # minimal biomass, below which predation rate is zero
  , tauP = 0.1 # slope of predation rate with biomass
  , kSP = 1/30 # 1/x years
  , pESP = 0.01  # /year production of phospatase by plants, small to check balance
  , nuPP = 0.05 # most of the phosphats to to mineral pool
)
parms0 <- within(parms0,{
  #kmR <- kmL <- km
  #eps1 <- eps2 <- eps
  #cnER <- cnEL <- cnE
  #kNR <- kNL <- kN
  plantNUpAbs <- iL / cnIL	# same litter input as plant uptake
  kINPlant <- plantNUpAbs <- 0			# no plant uptake
  lP <- lN       # leaching rate of inorganic P equals that of N
  nuP <- nuN     # mineralization of P during decomposiition equals that of N
  kIPPlant <- kINPlant  # plant uptake rate of P equals that of N
  iIP <- lP      # assume no P inputs compensate for leaching
  kLP <- kRP <- kSP
  pELP <- pERP <- pESP
})
# for compatibility set
# C:N:P ratio of cell walls to that of biomass
# all cell walls (all goes to R)
# all predation turnover starting from zero biomass
parms0C <- within(parms0,{
  cnBW <- cnB
  cpBW <- cpB
  cW <- 1.0
  B0 <- 0
  tauP <- tau
})

parms <- parms0C

x0 <- x0Orig <- c( #aE = 0.001*365
  BC = 20                     ##<< microbial biomass
  , RC = 7000                 ##<< N rich substrate
  , RN = 7000/parms0$cnIR    ##<< N rich substrate N pool
  , RP = 7000/parms0$cpIR
  , LC = 200                  ##<< N poor substrate
  , LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  , LP = 200/parms0$cpIL     ##<< N poor substrate P pool
  , IN =  1                   ##<< inorganic N pool
  , IP =  1                  ##<< inorganic P pool
  , alphaL = 0.5              ##<< initial community composition
  , alphaR = 0.5              ##<< initial community composition
  , alphaRP = 0              ##<< initial community composition
  , resp = 0                 ##<< cumulated respiration
  , leachN = 0               ##<< cumulated N losses by leaching
  , leachP = 0               ##<< cumulated P losses by leaching
)
x <- x0

mapValues <- function(
  ### replace some of the entries in a vector
  x       ##<< vector where replacement takes place
  , from  ##<< original values
  , to    ##<< replacement values of the same order as in from
) {
  ##details<< only values in from are replaced that are found in x
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}
getX0NoC <- function(x0){
  names(x0) <- mapValues(names(x0), c("BC","RC","LC"), c("B","R","L"))
  x0[setdiff(names(x0),c("resp","leachN","leachP"))]
}
getX0NoP <- function(x0){
  x0[setdiff(names(x0),c("RP","LP","IP","dRP","dLP","dIP","resp","leachN","leachP"))]
}
getX0SingleAlpha <- function(x0){
  names(x0) <- mapValues(names(x0), c("alphaR"), c("alpha"))
  x0[setdiff(names(x0),c("alphaL","alphaRP"))]
}
getX0NoP(x0)

x0Nlim <- c( #aE = 0.001*365
  BC = 20                     ##<< microbial biomass
  , RC = 700                 ##<< N rich substrate
  , RN = 700/parms0$cnIR    ##<< N rich substrate N pool
  , RP = 700/parms0$cpIR
  , LC = 2000                  ##<< N poor substrate
  , LN = 2000/parms0$cnIL     ##<< N poor substrate N pool
  , LP = 2000/parms0$cpIL     ##<< N poor substrate P pool
  , IN =  0                   ##<< inorganic N pool
  , IP =  10                  ##<< inorganic P pool
  , alphaL = 0.5              ##<< initial community composition
  , alphaR = 0.5              ##<< initial community composition
  , alphaRP = 0              ##<< initial community composition
  , resp = 0                 ##<< cumulated respiration
  , leachN = 0               ##<< cumulated N losses by leaching
  , leachP = 0               ##<< cumulated P losses by leaching
)

# Allocation Partitioning --------------------------------------------------


test_that("computeSesam4bAllocationPartitioning carbon limited",{
  # allocate only few to R (and most to L)
  cnL <- x["LC"]/x["LN"]
  cnR <- x["RC"]/x["RN"]
  cpL <- x["LC"]/x["LP"]
  cpR <- x["RC"]/x["RP"]
  # with nearly pure C limitation, should match with old partitioning
  limE <- c(C = 1, N = 1e-5, P = 1e-5)
  alpha0 <- c(L = 0.6, R = 0.4, LP = 1e-5, RP = 1e-5)
  limE <- limE/sum(limE); alpha0 <- alpha0/sum(alpha0)
  alpha_noCB <- computeSesam4bAllocationPartitioning_noCNB(
    dS = cbind(R = parms$kR*x["RC"], L = parms$kL*x["LC"]) [1,]
    ,dSP = cbind(RP = parms$kRP*x["RC"], LP = parms$kLP*x["LC"])[1,]
    , B = x["BC"]
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    ,alpha = alpha0
    ,limE = limE
    ,betaN = cbind(L = cnL, R = cpR, B = parms$cnB, E = parms$cnE)[1,]
    ,betaP = cbind(L = cpL, R = cpR, B = parms$cpB, E = parms$cpE)[1,]
  )
  alpha <- computeSesam4bAllocationPartitioning(
    dS = cbind(R = parms$kR*x["RC"], L = parms$kL*x["LC"]) [1,]
    ,dSP = cbind(RP = parms$kRP*x["RC"], LP = parms$kLP*x["LC"])[1,]
    , B = x["BC"]
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    ,alpha = alpha0
    ,limE = limE
    ,betaN = cbind(L = cnL, R = cpR, B = parms$cnB, E = parms$cnE)[1,]
    ,betaP = cbind(L = cpL, R = cpR, B = parms$cpB, E = parms$cpE)[1,]
  )
  expect_equal(names(alpha), c("L", "R", "LP", "RP"))
  expect_equal(sum(alpha),1)
  alphaC <- structure(computeSesam3sAllocationPartitioning(
    parms$kR*x["RC"], parms$kL*x["LC"],B = x["BC"]
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    ,alpha = alpha0["R"]
  ), names = "R")
  expect_equal(alpha_noCB["R"], alphaC, tolerance = 0.05 )
})

test_that("computeSesam4bAllocationPartitioning nitrogen limited",{
  # allocate only few to R (and most to L)
  cnL <- x["LC"]/x["LN"]
  cnR <- x["RC"]/x["RN"]
  cpL <- x["LC"]/x["LP"]
  cpR <- x["RC"]/x["RP"]
  # with nearly pure N limitation, should match with old partitioning
  limE <- c(C = 1e-5, N = 1, P = 1e-5)
  alpha0 <- c(L = 0.6, R = 0.4, LP = 1e-5, RP = 1e-5)
  limE <- limE/sum(limE); alpha0 <- alpha0/sum(alpha0)
  alpha <- computeSesam4bAllocationPartitioning(
    dS = cbind(R = parms$kR*x["RC"], L = parms$kL*x["LC"]) [1,]
    ,dSP = cbind(RP = parms$kRP*x["RC"], LP = parms$kLP*x["LC"])[1,]
    , B = x["BC"]
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    ,alpha = alpha0
    ,limE = limE
    ,betaN = cbind(L = cnL, R = cnR, B = parms$cnB, E = parms$cnE)[1,]
    ,betaP = cbind(L = cpL, R = cpR, B = parms$cpB, E = parms$cpE)[1,]
  )
  alpha_noCNB <- computeSesam4bAllocationPartitioning_noCNB(
    dS = cbind(R = parms$kR*x["RC"], L = parms$kL*x["LC"]) [1,]
    ,dSP = cbind(RP = parms$kRP*x["RC"], LP = parms$kLP*x["LC"])[1,]
    , B = x["BC"]
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    ,alpha = alpha0
    ,limE = limE
    ,betaN = cbind(L = cnL, R = cnR, B = parms$cnB, E = parms$cnE)[1,]
    ,betaP = cbind(L = cpL, R = cpR, B = parms$cpB, E = parms$cpE)[1,]
  )
  expect_equal(names(alpha), c("L", "R", "LP", "RP"))
  expect_equal(sum(alpha),1)
  alphaN <- structure(computeSesam3sAllocationPartitioning(
    parms$kR*x["RC"]/cnR, parms$kL*x["LC"]/cnL,B = x["BC"]
    ,kmkN = parms$km*parms$kN, aE =  parms$aE
    ,alpha = alpha0["R"]
  ), names = "R")
  expect_equal(alpha_noCNB["R"], alphaN, tolerance = 0.05 )
})

test_that("computeElementLimitations",{
  #computeElementLimitations()
  scen <- "C - limited"
  CsynBE <- c(C = 40, N = 6000, P = 6000)
  limE <- computeElementLimitations(
    CsynBE, tauB = 1
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )#, ce, eps  )
  expect_equal(names(limE), names(CsynBE))
  expect_equal(limE["C"], c(C = 1), tolerance = 1e-2)
  #
  scen <- "N - limited"
  CsynBE <- c(C = 6000, N = 40, P = 6000)
  limE <- computeElementLimitations(
    CsynBE, tauB = 1
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )#, ce, eps  )
  expect_equal(limE["N"], c(N = 1), tolerance = 1e-2)
  #
  scen <- "P - limited"
  CsynBE <- c(C = 6000, N = 6000, P = 40)
  limE <- computeElementLimitations(
    CsynBE, tauB = 1
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )#, ce, eps  )
  expect_equal(limE["P"], c(P = 1), tolerance = 1e-2)
  #
  scen <- "CN - limited"
  CsynBE <- c(C = 40, N = 40, P = 6000)
  limE <- computeElementLimitations(
    CsynBE, tauB = 1
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )#, ce, eps  )
  expect_equivalent(sum(limE[c("C","N")]), 1, tolerance = 1e-2)
  expect_equivalent(limE["C"], limE["N"], tolerance = 1e-2)
  # really want exact co-limitaiton at same CsynBE
  #expect_equivalent(limE["C"], limE["N"]/parms$cnB, tolerance = 1e-2)
  #
  scen <- "CN - limited, slightly more N limted"
  CsynBE <- c(C = 40, N = 40 - 1e-1, P = 6000)
  limE <- computeElementLimitations(
    CsynBE, tauB = 1
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )#, ce, eps  )
  limE
  expect_equivalent(sum(limE[c("C","N")]), 1, tolerance = 1e-2)
  expect_true(limE["N"] > limE["C"])
  #
  scen <- "NP - limited"
  CsynBE <- c(C = 6000, N = 40, P = 40)
  limE <- computeElementLimitations(
    CsynBE, tauB = 1
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )#, ce, eps  )
  expect_equivalent(sum(limE[c("N","P")]), 1, tolerance = 1e-2)
  expect_equivalent(limE["P"], limE["N"], tolerance = 1e-2)
  #
  scen <- "equal co-limitation"
  CsynBE <- c(C = 40, N = 40, P = 40)
  limE <- computeElementLimitations(
    CsynBE, tauB = 1
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )#, ce, eps  )
  expect_equivalent(limE["C"], limE["N"], tolerance = 1e-2)
  expect_equivalent(limE["C"], limE["P"], tolerance = 1e-2)
  expect_equivalent(sum(limE), 1, tolerance = 1e-2)
  #
  .tmp.f <- function(){
    epsNs <- seq(-0.4,+0.4, length.out = 51)
    limN <- sapply(epsNs, function(epsN){
      CsynBE <- c(C = 40, N = 40 - epsN, P = 6000)
      alphaBalanced <- computeElementLimitations(
        CsynBE, tauB = 1, delta = 5
        , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
      )["N"]
      #CsynBE, tauB = 1)#, ce, eps  )
    })
    plot(limN ~ epsNs, type = "lN")
    limN10 <- sapply(epsNs, function(epsN){
      CsynBE <- c(C = 40, N = 40 - epsN, P = 6000)
      alphaBalanced <- computeElementLimitations(
        CsynBE, tauB = 1, delta = 10
        , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
      )["N"]
    })
    lines(limN10 ~ epsNs, lty = "dashed")
    limN20 <- sapply(epsNs, function(epsN){
      CsynBE <- c(C = 40, N = 40 - epsN, P = 6000)
      alphaBalanced <- computeElementLimitations(
        CsynBE, tauB = 1, delta = 20
        , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
      )["N"]
    })
    lines(limN20 ~ epsNs, lty = "dotted")
  }
  #
  scen <- "strongly C - limited (negative synthesis)"
  # negative carbon balance
  CsynBE <- c(C = -40, N = 4, P = 4)
  limE <- computeElementLimitations(
    CsynBE, tauB = 1
    , betaB = c(C=1, N = parms$cnB, P = parms$cpB)
  )#, ce, eps  )
  limE
  # think of accounting for negative biomass balance
  expect_equivalent(limE["C"], 1, tolerance = 1e-2)
})


# CN limitations----------------------------------------------------
testParmsScen <- function(parmsInit){
  parmsInit$use_noCNB_return <- TRUE
  ans0 <- derivSesam4b(0, x0, parms = within(
    parmsInit, {tauP <- tau/x0["BC"]; tau <- 0}))
  ans0E <- derivSesam3a(0, getX0SingleAlpha(getX0NoC(x0)), parms = within(
    parmsInit,{epsTvr <- epsPred}))
  names4A <- names(getX0SingleAlpha(getX0NoP(ans0[[1]])))
  expect_equal(getX0SingleAlpha(getX0NoP(ans0[[1]])),
               structure(ans0E[[1]], names = names4A), tolerance = 1e-6)
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0SingleAlpha(getX0NoP(getX0NoC(x0))), times, derivSesam3a
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0, times, derivSesam4b, parms = within(
      parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  rbind(getX0SingleAlpha(getX0NoP(getX0NoC(xETest[2:16]))), xEExp[2:8])
  getX0SingleAlpha(getX0NoP(getX0NoC(xETest[2:16]))) - xEExp[2:8]
  expect_equal(  getX0SingleAlpha(getX0NoP(getX0NoC(xETest[2:16]))), xEExp[2:8],
                 tolerance = 1e-1)
  #
  # N limitation
  ans0 <- derivSesam4b(0, x0Nlim, parms = within(
    parmsInit, {tauP <- tau/x0["BC"]; tau <- 0}))
  ans0E <- derivSesam3a(0, getX0SingleAlpha(getX0NoC(x0Nlim)), parms = within(
    parmsInit,{epsTvr <- epsPred}))
  names4A <- names(getX0SingleAlpha(getX0NoP(ans0[[1]])))
  rbind(getX0SingleAlpha(getX0NoP(ans0[[1]])),
        structure(ans0E[[1]], names = names4A))
  expect_equal(getX0SingleAlpha(getX0NoP(ans0[[1]])),
               structure(ans0E[[1]], names = names4A)
               , tolerance = 1e-4)
  #times <- seq(0,2100, length.out = 101)
  #times <- seq(0,0.2, length.out = 1001)
  #times <- seq(0,2, length.out = 1001)
  resExp <- as.data.frame(lsoda(
    getX0SingleAlpha(getX0NoP(getX0NoC(x0Nlim))), times, derivSesam3a
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Nlim, times, derivSesam4b
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  #xETest <- unlist(slice(resTest,842))
  #xETest <- unlist(slice(resTest,841))
  expect_equivalent( xETest["BC"], xEExp["B"], tolerance = 1e-4)
  expect_equivalent( xETest["alphaR"], xEExp["alpha"], tolerance = 1e-1)
  expect_equivalent( xETest["alphaTargetR"], xEExp["alphaTarget"],
                     tolerance = 1e-1)
  expect_equal( getX0SingleAlpha(getX0NoP(getX0NoC(xETest[2:15]))), xEExp[2:8],
                tolerance = 1e-1)
  #
  # NP co-limitation
  # make R more attractive for phosphorous limitation, omit biomineralization
  # near colimitation almost no change in alpha
  parmsInitPlim <- within(parmsInit,{ cpIL <- 220; cpIR <- 20; cpB <- cpBW <- 50; kLP <- kRP <- 0})
  #parmsInitPlim <- within(parmsInit,{ cpIL <- 330; cpIR <- 20; cpB <- cpBW <- 50})
  x0Plim <- x0Nlim; x0Plim["IP"] <- 0;
  x0Plim["LP"] <- x0Plim["LC"]/parmsInitPlim$cpIL
  x0Plim["RP"] <- x0Plim["RC"]/parmsInitPlim$cpIR
  ans0 <- derivSesam4b(0, x0Plim, parms = within(
    parmsInitPlim, {tauP <- tau/x0["BC"]; tau <- 0}))
  ans0E <- derivSesam3a(0, getX0SingleAlpha(getX0NoC(x0Plim)), parms = within(
    parmsInitPlim,{epsTvr <- epsPred}))
  c(ans0[[2]]["alphaTargetR"], ans0E[[2]]["alphaTarget"])
  expect_true(ans0[[2]]["alphaTargetR"] / ans0E[[2]]["alphaTarget"] > 1.1)
  #times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0SingleAlpha(getX0NoP(getX0NoC(x0Plim))), times, derivSesam3a
    , parms = within(parmsInitPlim, {epsTvr <- epsPred})))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Plim, times, derivSesam4b
    , parms = within(parmsInitPlim,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  c(xETest["alphaTargetR"], xEExp["alphaTarget"])
  c(xETest["BC"], xEExp["B"])
  #expect_true( xETest["alphaR"] > xEExp["alpha"]) # not with substrate feedback
  expect_true( xETest["BC"] < xEExp["B"]) # more limited
  #
  # Microbial starvation (negative biomass synthesis)
  x0Starv <- x0; x0Starv["BC"] <- 80 # initially very high biomass
  times <- seq(0, 2100, length.out = 2)
  # due to predation increases with biomass, biomass decline is faster with v4
  #times <- seq(0,2, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0SingleAlpha(getX0NoP(getX0NoC(x0Starv))), times, derivSesam3a
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Starv, times, derivSesam4b
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  expect_equal( getX0SingleAlpha(getX0NoP(getX0NoC(xETest[2:15]))),
                xEExp[2:8], tolerance = 1e-1)
}

### plotting results --------------------------------
.tmp.f <- function(){
  # plots of results
  library(dplyr)
  library(ggplot2)
  # res <- suppressWarnings(bind_rows(
  #   cbind(resTest, scen = "Test"), cbind(resTest2, scen = "Test2"), cbind(resExp, scen = "Exp")))
  res <- bind_rows(
    cbind(resTest, scen = "Test")
    ,cbind(rename(resExp, BC = B, LC = L, RC = R, alphaR = alpha), scen = "Exp")
    )
  res <- bind_rows(
    cbind(resTest, scen = "Test")
    #,cbind(rename(resExp, BC = B, LC = L, RC = R, alphaR = alpha), scen = "Exp")
  )
  ggplot(filter(res, time >= 0), aes(time, alphaR, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time >= 0), aes(time, alphaTargetR, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time >= 0), aes(time, limN, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time >= 0), aes(time, limC, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time >= 0), aes(time, limP, color = scen)) + geom_line(alpha = 0.5)

  ggplot(filter(res, time >= 0), aes(time, BC, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time >= 0), aes(time, synB, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time > 0), aes(time, IP, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time < 500 & time > 0), aes(time, alpha, color = scen)) + geom_point()
  ggplot(filter(res, time < 500), aes(time, BC, color = scen)) + geom_line()
  ggplot(filter(res, time >= 0), aes(time, LC, color = scen)) + geom_line()
  ggplot(filter(res, time >= 0), aes(time, RC, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, respO, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, ER, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, EL, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, alphaC, color = scen)) + geom_line()
  ggplot(filter(res, time < 5000), aes(time, alphaN, color = scen)) + geom_line()
  ggplot(filter(res, time > 01), aes(time, PhiNB, color = scen, linetype = scen)) + geom_line()
}

.tmp.f <- function(){
  parmsTmp <- parmsInit
  parmsTmp <- parmsInitPlim
  # inspect derivative in steady state of C limitation
  ans0 <- derivSesam4b(0, xETest[2:15], parms = within(
    parmsTmp, {tauP <- tau/xEExp["B"]; tau <- 0}))
  ans0E <- derivSesam3a(0, getX0SingleAlpha(getX0NoC(xETest[2:15])),
                        parms = within(parmsTmp,{epsTvr <- epsPred}))
  c(ans0[[2]]["alphaTargetR"], ans0E[[2]]["alphaTarget"])
  ans0[[2]][c("limC","limN","limP")]
  ans0E[[2]][c("limC","limN")]

  names4A <- names(getX0SingleAlpha(getX0NoP(ans0[[1]])))
  rbind(getX0SingleAlpha(getX0NoP(ans0[[1]])),
        structure(ans0E[[1]], names = names4A))
}


test_that("same as sesam3s for fixed substrates", {
  parmsInit <- parmsFixedS <- within(parms0C,{
    isFixedS <- TRUE
  })
  testParmsScen(parmsInit)
})


test_that("mass balance also with feedback to DOM", {
  parmsInit <- within(parms0C, {isFixedI <- isFixedIP <- TRUE})
  expect_error(
    ans0 <- derivSesam4b(0, x0, parms = within(parmsInit, {cW <- 1; cnBW <- 10}))
  )
  ans0 <- derivSesam4b(0, x0, parms = parmsInit)
})


test_that("same as sesam2 with substrate feedbacks", {
  parmsInit <- within(parms0C, {isFixedI <- isFixedIP <- TRUE})
  testParmsScen(parmsInit)
})


.tmp.f <- function(){
  test_that("compare balanceAlphaBetweenElementLimitationsMin", {
#  Warning messages: 1: In lsoda(x0CNLim, times, derivSesam4b, parms =
#  within(parmsInit,  : an excessive amount of work (> maxsteps ) was done, but
#  integration was not successful - increase maxsteps 2: In lsoda(x0CNLim,
#  times, derivSesam4b, parms = within(parmsInit,  :
  parmsInit <- within(parms0C, {isFixedI <- TRUE})
  # from C to N limitation
  x0CNLim <- x0; x0CNLim["IN"] <- 0
  #times <- seq(0,2100, length.out = 2)
  times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0SingleAlpha(getX0NoP(x0CNLim)), times, derivSesam3a
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0CNLim, times, derivSesam4b
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  # yields timeout when going near col-limitation
  xETest <- unlist(tail(resTest,1))
  resTest2 <- as.data.frame(lsoda(
    x0CNLim, times, derivSesam4b
    , parms = within(parmsInit,{
      tauP <- tau/xEExp["B"]; tau <- 0; isBalanceAlphaMin <- TRUE})))
  xETest2 <- unlist(tail(resTest2,1))
  expect_equal( xETest["BC"], xEExp["B"], tolerance = 1e-6)
  expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( getX0SingleAlpha(getX0NoP(getX0NoC(xETest[2:15]))), xEExp[2:8], tolerance = 1e-6)
})
}


# P Biomineralization -----------------------------------------------------
parms0PBiomin <- within(parms0,{
  cpIL <- 220; cpIR <- 20; cpB <- cpBW <- 50; kLP <- kRP <- 1/30
})
x0P <- x0Nlim
x0P["alphaRP"] <- 0.1
x0P["alphaL"] <- x0P["alphaL"]*(1 - x0P["alphaRP"])
x0P["alphaR"] <- x0P["alphaR"]*(1 - x0P["alphaRP"])

testParmsScenBiomin <- function(parmsInit){
  parmsInit$use_noCNB_return <- TRUE
  ans0 <- derivSesam4b(0, x0P, parms = within(
    parmsInit, {tauP <- tau/x0P["BC"]; tau <- 0}))
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  #times <- seq(0,2, length.out = 101)
  resTest <- as.data.frame(lsoda(
    x0, times, derivSesam4b, parms = within(
      parmsInit,{tauP <- tau/x0P["BC"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  #
  # N limitation
  ans0 <- derivSesam4b(0, x0Nlim, parms = within(
    parmsInit, {tauP <- tau/x0P["BC"]; tau <- 0}))
  resTest <- as.data.frame(lsoda(
    x0Nlim, times, derivSesam4b
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  #
  # NP co-limitation
  # make R more attractive for phosphorous limitation, omit biomineralization
  # near colimitation almost no change in alpha
  parmsInitPlim <- within(parmsInit,{ cpIL <- 220; cpIR <- 20; cpB <- cpBW <- 50; kLP <- kRP <- 0})
  #parmsInitPlim <- within(parmsInit,{ cpIL <- 330; cpIR <- 20; cpB <- cpBW <- 50})
  x0Plim <- x0Nlim; x0Plim["IP"] <- 0;
  x0Plim["LP"] <- x0Plim["LC"]/parmsInitPlim$cpIL
  x0Plim["RP"] <- x0Plim["RC"]/parmsInitPlim$cpIR
  ans0 <- derivSesam4b(0, x0Plim, parms = within(
    parmsInitPlim, {tauP <- tau/x0P["BC"]; tau <- 0}))
  #times <- seq(0,2100, length.out = 101)
  resTest <- as.data.frame(lsoda(
    x0Plim, times, derivSesam4b
    , parms = within(parmsInitPlim,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
}

test_that("same as sesam2 with biomineralization", {
  ans0 <- derivSesam4b(0, x0P, parms = within(
    parms0PBiomin, {tauP <- tau/x0P["BC"]; tau <- 0}))
  parmsInit <- within(parms0PBiomin,{isFixedS <- TRUE })
  #testParmsScenBiomin(parmsInit)  # TODO
})




#require(testthat)
#test_file("tests/testthat/test_modSesam4a.R")
context("modSesam4a")

parms0 <- list(
  cnB = 7.16
  , cnBW = 10    ##<< C/N ratio of cell walls (that go to R, )
  ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
  #,cnE = 7.16
  ,cnIR = 4.5     ##<< between micr and enzyme signal
  ,cnIL = 30      ##<< N poor substrate
  #,kN = 0.05     ##<< (per day) enzyme turnover
  ,kN = 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
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
  ,km = 0.3       ##<< enzyme half-saturation constant
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
  ,kIPlant = 10.57 #0.0289652*365         ##<< plant uptake iP I
  ,iB = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iB I
  ,iI = 0         ##<< input of mineral N
  ,l = 0.96 #0.00262647*365       ##<< leaching rate of mineralN l I
  , nu = 0.9       ##<< microbial N use efficiency
  #, isEnzymeMassFlux = FALSE  ##<< steady state B solution neglects enyzme mass fluxes
  , isEnzymeMassFlux = TRUE  ##<< steady state B solution accounts for enyzme mass fluxes
  , nuP = 0.3      ##<< microbial uptake of depolymerized P, (1-nuP) is mineralized
  , cpE = 50
  , cpB = 40
  , cpBW = 50
  , cpIR = 40
  , cpIL = 40*3
  , iBP = 0.38 * 10.57 # start with same as N
  , cW = 0.5 # proportion of cell wall within microbial biomass
  , B0 = 0  # minimal biomass, below which predation rate is zero
  , tauP = 0.1 # slope of predation rate with biomass
)
parms0 <- within(parms0,{
  kmR <- kmL <- km
  eps1 <- eps2 <- eps
  cnER <- cnEL <- cnE
  kNR <- kNL <- kN
  kIPlant <- iL / cnIL	# same litter input as plant uptake
  kIPlant <- 0			# no plant uptake
  lP <- l       # leaching rate of inorganic P equals that of N
  nuP <- nu     # mineralization of P during decomposiition equals that of N
  kIPPlant <- kIPlant  # plant uptake rate of P equals that of N
  iIP <- l      # assume no P inputs compensate for leaching
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
  B = 20                     ##<< microbial biomass
  , R = 7000                 ##<< N rich substrate
  , RN = 7000/parms0$cnIR    ##<< N rich substrate N pool
  , RP = 7000/parms0$cpIR
  , L = 200                  ##<< N poor substrate
  , LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  , LP = 200/parms0$cpIL     ##<< N poor substrate P pool
  , I =  1                   ##<< inorganic N pool
  , IP =  1                  ##<< inorganic P pool
  , alpha = 0.5              ##<< initial community composition
)
x <- x0
getX0NoP <- function(x0){x0[setdiff(names(x0),c("RP","LP","IP","dRP","dLP","dIP"))]}
getX0NoP(x0)

x0Nlim <- c( #aE = 0.001*365
  B = 20                     ##<< microbial biomass
  , R = 1000                 ##<< N rich substrate
  , RN = 1000/parms0$cnIR    ##<< N rich substrate N pool
  , RP = 1000/parms0$cpIR
  , L = 200                  ##<< N poor substrate
  , LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  , LP = 200/parms0$cpIL     ##<< N poor substrate P pool
  , I =  0                   ##<< inorganic N pool
  , IP =  10                  ##<< inorganic P pool
  , alpha = 0.5              ##<< initial community composition
)

testParmsScen <- function(parmsInit){
  ans0 <- derivSesam4a(0, x0, parms = within(
    parmsInit, {tauP <- tau/x0["B"]; tau <- 0}))
  ans0E <- derivSesam3a(0, x0, parms = within(
    parmsInit,{epsTvr <- epsPred}))
  expect_equal(getX0NoP(ans0[[1]]), ans0E[[1]], tolerance = 1e-6)
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0), times, derivSesam3s
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0, times, derivSesam4a
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-6)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  xEExp2 <- xETest[2:11]; xEExp2[names(xEExp[c(2:8)])] <- xEExp[c(2:8)]
  .tmp.f <- function(){
    derivSesam4a(0, xETest[c(2:11)], within(parmsInit, {isRecover <- TRUE; epsTvr <- epsPred}))
    derivSesam4a(0, xEExp2, within(parmsInit, {isRecover <- TRUE; tauP <- tau/xEExp["B"]; tau <- 0}))
    derivSesam3a(0, xEExp[c(2:8)], within(parmsInit, {isRecover <- TRUE; epsTvr <- epsPred}))
  }
  #
  # N limitation
  #times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0Nlim), times, derivSesam3s
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Nlim, times, derivSesam4a
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-6)
  expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  #
  # NP col-limitation
  x0Plim <- x0Nlim; x0Plim["IP"] <- 0
  #times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0Plim), times, derivSesam3s
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Plim, times, derivSesam4a
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  expect_true( xETest["alpha"] < xEExp["alpha"])
  # interestingly this leads to slightly higher biomass under colimitation
  # expect_true( xETest["B"] < xEExp["B"])
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-2)
  #expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
}

.tmp.f <- function(){
  # plots of results
  library(dplyr)
  library(ggplot2)
  res <- suppressWarnings(bind_rows(
    cbind(resTest, scen = "Test"), cbind(resExp, scen = "Exp")))
  ggplot(filter(res, time > 1), aes(time, B, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time > 0), aes(time, alpha, color = scen)) + geom_point(alpha = 0.5)
  ggplot(filter(res, time < 500 & time > 0), aes(time, alpha, color = scen)) + geom_point()
  ggplot(filter(res, time < 500), aes(time, B, color = scen)) + geom_line()
  ggplot(filter(res, time >= 0), aes(time, L, color = scen)) + geom_line()
  ggplot(filter(res, time >= 0), aes(time, R, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, respO, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, ER, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, EL, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, alphaC, color = scen)) + geom_line()
  ggplot(filter(res, time < 5000), aes(time, alphaN, color = scen)) + geom_line()
  ggplot(filter(res, time > 01), aes(time, PhiB, color = scen, linetype = scen)) + geom_line()
}

test_that("same as sesam3s for fixed substrates", {
  parmsFixedS <- within(parms0C,{
    isFixedS <- TRUE
  })
  testParmsScen(parmsFixedS)
})


test_that("mass balance also with feedback to DOM", {
  parmsInit <- within(parms0C, {isFixedI <- TRUE})
  ans0 <- derivSesam4a(0, x0, parms = within(parmsInit, {cW <- 1; B0 <- 0}))
  ans0 <- derivSesam4a(0, x0, parms = within(parmsInit, {cW <- 0; B0 <- 1e20}))
  ans0 <- derivSesam4a(0, x0, parms = parmsInit)
})


test_that("same as sesam2 with substrate feedbacks", {
  parmsInit <- within(parms0C, {isFixedI <- TRUE})
  testParmsScen(parmsInit)
})





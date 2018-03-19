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

# setup splitting all pools into SOM and amend
elements <- c("C","N","P")
units <- structure(lapply(elements, function(el) c(SOM = 1, amend = 1))
                   , names = elements)
parms0 <- within(parms0,{
  multiPoolFractions <- createMultiPoolFractions(
    units, setX = createSesam4setX(units))
  relIL <- c(SOM = 0, amend = 1)
  relIR <- c(SOM = 0, amend = 1)
})
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

x0 <- x0Orig <- c( #aE = 0.001*365
  B_SOM = 20                     ##<< microbial biomass
  , R_SOM = 7000                 ##<< N rich substrate
  , L_SOM = 200                  ##<< N poor substrate
  , RN_SOM = 7000/parms0$cnIR    ##<< N rich substrate N pool
  , LN_SOM = 200/parms0$cnIL     ##<< N poor substrate N pool
  , I_SOM =  1                   ##<< inorganic N pool
  , RP_SOM = 7000/parms0$cpIR
  , LP_SOM = 200/parms0$cpIL     ##<< N poor substrate P pool
  , IP_SOM =  1                  ##<< inorganic P pool
  , B_amend = 0                     ##<< microbial biomass
  , R_amend = 0                 ##<< N rich substrate
  , L_amend = 0
  , RN_amend = 0
  , LN_amend = 0
  , I_amend =  0
  , RP_amend = 0
  , LP_amend = 0
  , IP_amend =  0
  , alpha = 0.5              ##<< initial community composition
)
x <- x0
getX0Sum <- function(x0){
  x <- parms0$multiPoolFractions
  x <- x$setX(x, x0)
  x$tot
}
x0Sum <- getX0Sum(x0)
test_that("equal initial sum of pools",{
  expect_equivalent(x0Sum["alpha"], x0["alpha"])
  expect_equivalent(x0Sum[1:9], x0[1:9])
})


x0Nlim <- x0
.R <- 1000
x0Nlim["R_SOM"] = .R
x0Nlim["RN_SOM"] = .R/parms0$cnIR    ##<< N rich substrate N pool
x0Nlim["RP_SOM"] = .R/parms0$cpIR
x0NlimSum <- getX0Sum(x0Nlim)

#parmsInit <- parms0
testParmsScen <- function(parmsInit){
  ans0 <- derivSesam4F(0, x0, parms = parmsInit)
  ans0E <- derivSesam4a(0, x0Sum, parms = parmsInit)
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
  parmsInitPlim <- within(parmsInit, cpIL <- 160)
  #times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0Plim), times, derivSesam3s
    , parms = within(parmsInitPlim, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Plim, times, derivSesam4a
    , parms = within(parmsInitPlim,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  expect_true( xETest["alpha"] < xEExp["alpha"])
  # interestingly this leads to slightly higher biomass under colimitation
  # expect_true( xETest["B"] < xEExp["B"])
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-2)
  #expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  #
  # Microbial starvation
  x0Starv <- x0; x0Starv["B"] <- 80
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0Starv), times, derivSesam3s
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Starv, times, derivSesam4a
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  xETest <- unlist(tail(resTest,1))
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-6)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  xEExp2 <- xETest[2:11]; xEExp2[names(xEExp[c(2:8)])] <- xEExp[c(2:8)]
}

.tmp.f <- function(){
  # plots of results
  library(dplyr)
  library(ggplot2)
  # res <- suppressWarnings(bind_rows(
  #   cbind(resTest, scen = "Test"), cbind(resTest2, scen = "Test2"), cbind(resExp, scen = "Exp")))
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
  parmsFixedS <- within(parms0,{
    isFixedS <- TRUE
  })
  testParmsScen(parmsFixedS)
})


test_that("mass balance also with feedback to DOM", {
  parmsInit <- within(parms0, {isFixedI <- TRUE})
  ans0 <- derivSesam4a(0, x0, parms = within(parmsInit, {cW <- 1; B0 <- 0}))
  ans0 <- derivSesam4a(0, x0, parms = within(parmsInit, {cW <- 0; B0 <- 1e20}))
  ans0 <- derivSesam4a(0, x0, parms = parmsInit)
})


test_that("same as sesam2 with substrate feedbacks", {
  parmsInit <- within(parms0, {isFixedI <- TRUE})
  testParmsScen(parmsInit)
})

.tmp.f <- function(){ test_that("balanceAlphaBetweenElementLimitationsMin", {
  alphas <- c(0.2,0.6,0.8)
  CsynBEs <- c(40,20,20)
  ans <- balanceAlphaBetweenElementLimitationsMin(alphas, CsynBEs)
  expect_equal( ans$alpha, alphas[2L] )
  expect_equal( ans$wELim, c(0,1,0))
})
}

.tmp.f <- function(){ test_that("compare balanceAlphaBetweenElementLimitationsMin", {
#  Warning messages: 1: In lsoda(x0CNLim, times, derivSesam4a, parms =
#  within(parmsInit,  : an excessive amount of work (> maxsteps ) was done, but
#  integration was not successful - increase maxsteps 2: In lsoda(x0CNLim,
#  times, derivSesam4a, parms = within(parmsInit,  :
  parmsInit <- within(parms0, {isFixedI <- TRUE})
  # from C to N limitation
  x0CNLim <- x0; x0CNLim["I"] <- 0
  #times <- seq(0,2100, length.out = 2)
  times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0CNLim), times, derivSesam3s
    , parms = within(parmsInit, {epsTvr <- epsPred})))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0CNLim, times, derivSesam4a
    , parms = within(parmsInit,{tauP <- tau/xEExp["B"]; tau <- 0})))
  # yields timeout when going near col-limitation
  xETest <- unlist(tail(resTest,1))
  resTest2 <- as.data.frame(lsoda(
    x0CNLim, times, derivSesam4a
    , parms = within(parmsInit,{
      tauP <- tau/xEExp["B"]; tau <- 0; isBalanceAlphaMin <- TRUE})))
  xETest2 <- unlist(tail(resTest2,1))
  expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-6)
  expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
})
}

.tmp.f <- function(){
  ans1 <- derivSesam4a(0, x0, parms0)
  times <- seq(0,500, length.out = 101)
  res <- res1 <- as.data.frame(lsoda( x0, times, derivSesam4a, parms = parms0))
}



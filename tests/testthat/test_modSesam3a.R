#require(testthat)
#test_file("tests/testthat/test_modSesam3a.R")
context("modSesam3a")

parms0 <- list(
  cnB = 7.16
  ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
  #,cnE = 7.16
  ,cnIR = 4.5     ##<< between micr and enzyme signal
  ,cnIL = 30      ##<< N poor substrate
  #,kN = 0.05     ##<< (per day) enzyme turnover
  ##,kN = 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
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
  ,epsTvr = 0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
  ,iR = 0         ##<< input modelled explicitely
  ,iL = 300       ##<< g/m2 input per year (half NPP)
  ,useFixedAlloc = FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
  ,kINPlant = 10.57 #0.0289652*365         ##<< plant uptake iP IN
  ,iBN = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iBN IN
  ,iIN = 0         ##<< input of mineral N
  ,lN = 0.96 #0.00262647*365       ##<< leaching rate of mineralN lN IN
  ,nuN = 0.9       ##<< microbial N use efficiency
  , isEnzymeMassFlux = TRUE  ##<< set to FALSE to neglect C and N mass fluxes by enzymes
  , kmN = 0.3 * 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
)
parms0 <- within(parms0,{
  #kmR <- kmL <- km
  #eps1 <- eps2 <- eps
  #cnER <- cnEL <- cnE
  #kNR <- kNL <- kN
  plantNUpAbs <- iL / cnIL	# same litter input as plant uptake
  kINPlant <- plantNUpAbs <- 0			# no plant uptake
})

parms <- parms0

x0 <- x0Orig <- c( #aE = 0.001*365
  B = 20                     ##<< microbial biomass
  , R = 7000                 ##<< N rich substrate
  , RN = 7000/parms0$cnIR    ##<< N rich substrate N pool
  , L = 200                  ##<< N poor substrate
  , LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  , IN =  1                   ##<< inorganic pool
  , alpha = 0.5              ##<< initial community composition
)
x <- x0

x0Nlim <- c( #aE = 0.001*365
  B = 20                     ##<< microbial biomass
  , R = 1000                 ##<< N rich substrate
  , RN = 1000/parms0$cnIR    ##<< N rich substrate N pool
  , L = 200                  ##<< N poor substrate
  , LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  , IN =  0                   ##<< inorganic pool
  , alpha = 0.5              ##<< initial community composition
)

test_that("same as seam for fixed substrates", {
  #x0["L"] <- calcLSteady5( parms0$iL, x0["R"], parms=parms0 )
  #x0["B"] <- calcBSteady5( parms0$iL, x0["R"], L=x0["L"], parms=parms0 )
  #x0[c("ER","EL")] <- calcESteady5( x0["B"], parms=parms0 )

  parmsFixedS <- within(parms0,{
    isFixedS <- TRUE
  })
  parmsTvr0 <- within(parms0,{
    isTvrNil <- TRUE
    iR <- 160
  })
  parmsInit <- parms0
  #parmsInit <- within(parms0, {isFixedS <- TRUE})
  #parmsInit <- within(parms0, {isAlphaFix <- TRUE})
  #parmsInit <- within(parms0, {isAlphaMatch <- TRUE})
  res <- derivSesam3a(0, x0, parmsInit)
  resExp <- derivSeam3a(
    0, getSeamStateFromSesam(x0, getKmKNParsFromKmN(parmsInit))
    , getKmKNParsFromKmN(parmsInit))
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  resTest <- as.data.frame(lsoda( x0, times, derivSesam3a, parms = parmsFixedS))
  resTestNoEnz <- as.data.frame(lsoda(
    x0, times, derivSesam3a
    , parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  #resExpOld <- as.data.frame(lsoda( x0Seam3, times, derivSeam3a, parms = parmsFixedS))
  resExp <- as.data.frame(lsoda(
    getSeamStateFromSesam(x0, getKmKNParsFromKmN(parmsFixedS))
    , times, derivSeam3a, parms = getKmKNParsFromKmN(parmsFixedS)))
  xETest <- unlist(tail(resTest,1))
  xETestNoEnz <- unlist(tail(resTestNoEnz,1))
  xEExp <- unlist(tail(resExp,1)); xEExp[-(3:4)]
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-6)
  expect_equal( xETestNoEnz["alphaC"], xEExp["alphaC"], tolerance = 1e-3)
  expect_equal( xETest[2:7], xEExp[c(2,5:9)], tolerance = 1e-6)
  expect_equal( xETestNoEnz[2:7], xEExp[c(2,5:9)], tolerance = 0.1)
  .tmp.f <- function(){
    derivSeam3a(0, xEExp[2:9], within(parmsInit, isRecover <- TRUE))
    derivSesam3a(0, xEExp[c(2,5:9)], within(parmsInit, isRecover <- TRUE))
  }
  #
  # N limitation
  #times <- seq(0,2100, length.out = 101)
  resTest <- as.data.frame(lsoda( x0Nlim, times, derivSesam3a, parms = parmsFixedS))
  resTestNoEnz <- as.data.frame(lsoda(
    x0Nlim, times, derivSesam3a
    , parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  #resExpOld <- as.data.frame(lsoda(
  #  x0NlimSeam3, times, derivSeam3a, parms = parmsFixedS))
  resExp <- as.data.frame(lsoda(
    getSeamStateFromSesam(x0Nlim, getKmKNParsFromKmN(parmsFixedS))
    , times, derivSeam3a, parms = getKmKNParsFromKmN(parmsFixedS)))
  xETest <- unlist(tail(resTest,1))
  xETestNoEnz <- unlist(tail(resTestNoEnz,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( xETest[2:7], xEExp[c(2,5:9)], tolerance = 1e-6)
  expect_equal( xETestNoEnz[2:7], xEExp[c(2,5:9)], tolerance = 0.1)
})

test_that("same as seam with substrate feedbacks", {
  parmsInit <- within(parms0, {isFixedI <- TRUE})
  times <- seq(0,800, length.out = 2)
  #times <- seq(0,800, length.out = 101)
  #times <- c(0,148:151)
  #times <- seq(0,2100, by = 2)
  #times <- seq(0,10000, length.out = 101)
  resTest <- as.data.frame(lsoda( x0, times, derivSesam3a, parms = parmsInit))
  resTestNoEnz <- as.data.frame(lsoda(
    x0, times, derivSesam3a
    , parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  resExp <- as.data.frame(lsoda(
    getSeamStateFromSesam(x0, getKmKNParsFromKmN(parmsInit))
    , times, derivSeam3a, parms = getKmKNParsFromKmN(parmsInit)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1));
  xETestNoEnz <- unlist(tail(resTestNoEnz,1))
  xpESteady <- unlist(head(tail(resTest,2),1))	# the previous before end
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-4)
  expect_equal( xETestNoEnz["alphaC"], xEExp["alphaC"], tolerance = 1e-2)
  expect_equal( xETest[2:7], xEExp[c(2,5:9)], tolerance = 1e-6)
  expect_equal( xETestNoEnz[2:7], xEExp[c(2,5:9)], tolerance = 0.1)
  rbind(xEExp[names(xETest)], xETest, xETestNoEnz)
  .tmp.f <- function(){
    derivSesam3a(0, xETest[2:8], within(parmsInit, isRecover <- TRUE))
    derivSeam3a(0, xEExp[2:10], within(parmsInit, isRecover <- TRUE))
    derivSesam3a(0, xETest, within(parmsInit, isRecover <- TRUE))
    derivSesam3a(0, xpESteady, within(parmsInit, isRecover <- TRUE))
  }
  #
  # N limitation
  resTest <- as.data.frame(lsoda( x0Nlim, times, derivSesam3a, parms = parmsInit))
  resTestNoEnz <- as.data.frame(lsoda(
    x0Nlim, times, derivSesam3a
    , parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  #resExp <- as.data.frame(lsoda( x0NlimSeam3, times, derivSeam3a, parms = parmsInit))
  resExp <- as.data.frame(lsoda(
    getSeamStateFromSesam(x0Nlim, getKmKNParsFromKmN(parmsInit))
    , times, derivSeam3a, parms = getKmKNParsFromKmN(parmsInit)))
  xETest <- unlist(tail(resTest,1))
  xETestNoEnz <- unlist(tail(resTestNoEnz,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-6)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( xETestNoEnz["B"], xEExp["B"], tolerance = 0.1)
  expect_equal( xETestNoEnz["alphaN"], xEExp["alphaN"], tolerance = 1e-2)
  # here alphaC will be different, because Crevenue is computed
  # with alphaC in sesam but with current alpha in seam
  expect_equal( xETest[2:7], xEExp[c(2,5:9)], tolerance = 1e-5)
  expect_equal( xETestNoEnz[2:7], xEExp[c(2,5:9)], tolerance = 0.1)
  rbind(xEExp[c(2,5:9)], xETest[2:7], xETestNoEnz[2:7])
  rbind(xEExp[names(xETest)], xETest, xETestNoEnz)
  #
  # from C to N limitation
  x0CNLim <- x0; x0CNLim["IN"] <- 0
  #x0CNLimSeam3 <- x0Seam3; x0CNLimSeam3["IN"] <- 0
  times <- seq(0,1200, length.out = 2)
  #times <- seq(0,1200, length.out = 101)
  #times <- c(0,seq(140,220, length.out = 101))
  resTest <- as.data.frame(lsoda( x0CNLim, times, derivSesam3a, parms = parmsInit))
  resTestNoEnz <- as.data.frame(lsoda(
    x0CNLim, times, derivSesam3a
    , parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  #resExp <- as.data.frame(lsoda( x0CNLimSeam3, times, derivSeam3a, parms = parmsInit))
  resExp <- as.data.frame(lsoda(
    getSeamStateFromSesam(x0CNLim, getKmKNParsFromKmN(parmsInit))
    , times, derivSeam3a, parms = getKmKNParsFromKmN(parmsInit)))
  xETest <- unlist(tail(resTest,1))
  xETestNoEnz <- unlist(tail(resTestNoEnz,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-5)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-5)
  expect_equal( xETestNoEnz["B"], xEExp["B"], tolerance = 0.1)
  expect_equal( xETestNoEnz["alphaN"], xEExp["alphaN"], tolerance = 0.1)
  expect_equal( xETest[2:7], xEExp[c(2,5:9)], tolerance = 1e-3)
  expect_equal( xETestNoEnz[2:7], xEExp[c(2,5:9)], tolerance = 0.1)
  rbind(xEExp[c(2,5:9)], xETest[2:7], xETestNoEnz[2:7])
  rbind(xEExp[names(xETest)], xETest, xETestNoEnz)
})

.tmp.f <- function(){
  # plots of results
  library(dplyr)
  library(ggplot2)
  res <- suppressWarnings(bind_rows(
    cbind(resTest, scen = "Test"), cbind(resTestNoEnz, scen = "TestNoEnz"), cbind(resExp, scen = "Exp")))
  ggplot(filter(res, time > 1), aes(time, B, color = scen)) + geom_line()
  ggplot(filter(res, time < 500 & time > 0), aes(time, alpha, color = scen)) + geom_point()
  ggplot(filter(res, time < 500), aes(time, B, color = scen)) + geom_line()
  ggplot(filter(res, time > 1), aes(time, L, color = scen)) + geom_line()
  ggplot(filter(res, time > 1), aes(time, R, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, respO, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, ER, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, EL, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, alphaC, color = scen)) + geom_line()
  ggplot(filter(res, time < 5000), aes(time, alphaN, color = scen)) + geom_line()
  ggplot(filter(res, time > 01), aes(time, PhiNB, color = scen, linetype = scen)) + geom_line()
}





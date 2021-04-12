#require(testthat)
context("modSeam3a")

parms0 <- list(
  cnB = 7.16
  ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
  #,cnE = 7.16
  ,cnIR = 4.5    ##<< between micr and enzyme signal
  ,cnIL = 30     ##<< N poor substrate
  #,kN = 0.05    ##<< (per day) enzyme turnover
  ,kN = 0.01*365 ##<< /yr enzyme turnover 1% turning over each day
  ,kNB = 0.8     ##<< amount of recycling enzyme turnover by biomass (added to uptake instead of R)
  #,kR = 0.2     ##<< substrate decomposition rate N-rich (here simulate large N stock)
  #,kL = 1       ##<< substrate decomposition rate N-poor
  #,kR = 5e-3    ##<< substrate decomposition rate N-rich (here simulate large N stock)
  #,kL = 10e-3   ##<< substrate decomposition rate N-poor
  #,kL = 5e-2    ##<< substrate decomposition rate N-poor
  #,aE = 0.05    ##<< C-uptake allocated to enzymes
  #,kR = 1/(5*365)    ##<< 5 years
  #,kL = 1/(0.5*365)  ##<< 1/2 years
  ,kR = 1/(50)        ##<< 1/(x years)
  ,kL = 1/(1)        ##<< 1/(x years)
  #,aE = 0.003*365   ##<< C biomass allocated to enzymes gC/day /microbial biomass
  ,aE = 0.001*365    ##<< C biomass allocated to enzymes gC/day /microbial biomass
  ,km = 0.3     ##<< enzyme half-saturation constant
  #,km = 0.03   ##<< enzyme half-saturation constant
  #,km = 14     ##<< enzyme half-saturation constant
  ,m = 0.02*365    ##<< maintenance respiration rate   gC/day /microbial biomass
  ,tau = 1/60*365  ##<< biomass turnover rate (12 days)
  ,eps = 0.5       ##<< carbon use efficiency
  ,epsTvr = 0.3    ##<< carbon use efficiency of microbial tvr (predators respire)
  ,iR = 0          ##<< input modelled explicitely
  ,iL = 300        ##<< g/m2 input per year (half NPP)
  #,plantNUp = 300/70*1/4  # plant N uptake balancing N inputs
  ,plantNUp = 0
  ,useFixedAlloc = FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
  ,kINPlant = 10.57 #0.0289652*365          ##<< plant uptake iP IN
  ,iBN = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iBN IN
  ,iIN = 0     ##<< input of mineral N
  ,lN = 0.96   #0.00262647*365       ##<< leaching rate of mineralN lN IN
	,nuN = 0.9     # microbial N use efficiency
)
parms0 <- within(parms0,{
  kmkN <- km*kN
  kmR <- kmL <- km
  eps1 <- eps2 <- eps
  cnER <- cnEL <- cnE
  kNR <- kNL <- kN
  plantNUpOrg <- 0 # no organic N uptake
  plantNUpAbs <- iL / cnIL	# inorganic plant N uptake same litter input
  kINPlant <- 1000  # high number -> constrain by plantNUpAbs
  kINPlant <- plantNUpAbs <- 0			# no plant uptake
})

parms <- parms0

x0 <- x0Orig <- c( #aE = 0.001*365
  B = 20                    ##<< microbial biomass
  ,ER  = 1.5*parms0$km    ##<< total enzyme pool
  ,EL  = 1.5*parms0$km    ##<< total enzyme pool
  ,R = 7000                 ##<< N rich substrate
  ,RN = 7000/parms0$cnIR    ##<< N rich substrate N pool
  ,L = 200                  ##<< N poor substrate
  ,LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  ,IN =  1                   ##<< inorganic pool
  ,alpha = 0.5              ##<< initial microbial community partitioning
)
x <- x0

x0Seam2 <- c(x0[setdiff(names(x0),"alpha")]
)[c("B","ER","EL","R","RN","L","LN","IN")]	##<< make sure to have same order as derivative of Seam2
#x0Seam2

x0Nlim <- c( #aE = 0.001*365
  B = 20                     ##<< microbial biomass
  ,ER  = 1.5*parms0$km     ##<< total enzyme pool
  ,EL  = 1.5*parms0$km     ##<< total enzyme pool
  ,R = 1000                 ##<< N rich substrate
  ,RN = 1000/parms0$cnIR    ##<< N rich substrate N pool
  ,L = 200                  ##<< N poor substrate
  ,LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  ,IN =  0                   ##<< inorganic pool
  ,alpha = 0.5              ##<< initial microbial community partitioning
)
x0NlimSeam2 <- c(x0Nlim[setdiff(names(x0),"alpha")]
)[c("B","ER","EL","R","RN","L","LN","IN")]	##<< make sure to have same order as derivative of Seam2


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
  res <- derivSeam3a(0, x0, parmsInit)
  # TODO test
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,5, length.out=101)
  resTest <- as.data.frame(lsoda(
    x0, times, derivSeam3a, parms = parmsFixedS
    , fBalanceAlpha = balanceAlphaBetweenCNLimitations))
  resExp <- as.data.frame(lsoda( x0Seam2, times, derivSeam2, parms = parmsFixedS))
  resExp$alphaTarget = resExp$alpha
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  rbind( xETest, xEExp[names(xETest)] )
  #sort( (xETest-xEExp[names(xETest)])/pmax(1e-10,xEExp[names(xETest)]))
  # relaxed precision because enzyme mass fluxes are different
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-2)
  .tmp.f <- function(){
    resComp <- rbind( cbind(scen = "Test", resTest)
                      , cbind(scen = "Exp", resExp))
    ggplot( resComp, aes(time, alpha, color = scen)) + geom_line()
    ggplot( resComp, aes(time, B, color = scen)) + geom_line()
  }
  # relaxed precision because enzyme mass fluxes are different
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-6)
  rbind( xETest, xEExp[names(xETest)] )
  expect_true( all(abs(xETest - xEExp[names(xETest)])/pmax(1e-10,xEExp[names(xETest)]) < 0.01))
  .tmp.f <- function(){
    derivSeam2(0, xEExp[2:9], within(parmsInit, isRecover <- TRUE))
    derivSeam3a(0, xEExp[c(2,5:9)], within(parmsInit, isRecover <- TRUE))
  }
  #
  # N limitation
  #times <- seq(0,2100, length.out = 101)
  resTest <- as.data.frame(lsoda(
    x0Nlim, times, derivSeam3a, parms = parmsFixedS
    , fBalanceAlpha = balanceAlphaBetweenCNLimitations))
  resExp <- as.data.frame(lsoda( x0NlimSeam2, times, derivSeam2, parms = parmsFixedS))
  resExp$alphaTarget <- resExp$alpha
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  rbind(xETest, xEExp[names(xETest)])
  #sort( (xETest-xEExp[names(xETest)])/pmax(1e-10,xEExp[names(xETest)]))
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-2)
  expect_true( all(abs(xETest - xEExp[names(xETest)])/pmax(1e-10,xEExp[names(xETest)]) < 0.01))
})

test_that("same as seam with substrate feedbacks", {
  parmsInit <- within(parms0, {isFixedI <- TRUE})
  times <- seq(0,800, length.out = 2)
  #times <- seq(0,800, length.out = 101)
  #times <- c(0,148:151)
  #times <- seq(0,2100, by = 2)
  #times <- seq(0,10000, length.out = 101)
  resTest <- as.data.frame(lsoda(
    x0, times, derivSeam3a, parms = parmsInit
    , fBalanceAlpha = balanceAlphaBetweenCNLimitations))
  resExp <- as.data.frame(lsoda( x0Seam2, times, derivSeam2, parms = parmsInit))
  resExp$alphaTarget <- resExp$alpha
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1));
  #xpESteady <- unlist(head(tail(resTest,2),1))	# the previous before end
  rbind(xETest, xEExp[names(xETest)])
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-6)
  # smaller then 10% relative derivation
  expect_true( all(abs(xETest - xEExp[names(xETest)])/pmax(1e-5,abs(xEExp[names(xETest)])) < 0.2))
  #
  # N limitation
  times <- seq(0,800, length.out = 8)
  #times <- seq(0,800, length.out = 101)
  resTest <- as.data.frame(lsoda(
    x0Nlim, times, derivSeam3a, parms = parmsInit
    , fBalanceAlpha = balanceAlphaBetweenCNLimitations))
  resExp <- as.data.frame(lsoda( x0NlimSeam2, times, derivSeam2, parms = parmsInit))
  resExp$alphaTarget <- resExp$alpha
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  rbind(xETest, xEExp[names(xETest)])
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  # relative differences get large for small changes
  namesCompare <- setdiff( names(xETest), c("dL","dR","dB"))
  expect_true( all(abs(xETest[namesCompare] - xEExp[namesCompare])/pmax(1e-5,abs(xEExp[namesCompare])) < 1e-2))
  # from C to N limitation
  x0CNLim <- x0; x0CNLim["IN"] <- 0
  x0CNLimSeam2 <- x0Seam2; x0CNLimSeam2["IN"] <- 0
  times <- seq(0,1200, length.out = 8)
  #times <- seq(0,1200, length.out = 101)
  #times <- c(0,seq(140,220, length.out = 101))
  resTest <- as.data.frame(lsoda( x0CNLim, times, derivSeam3a, parms = parmsInit))
  resExp <- as.data.frame(lsoda( x0CNLimSeam2, times, derivSeam2, parms = parmsInit))
  resExp$alphaTarget <- resExp$alpha
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  rbind(xETest, xEExp[names(xETest)])
  expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 5e-2)
  # larger relative differences in dR
  namesCompare <- setdiff( names(xETest), c("dL","dR","dB"))
  expect_true( all(abs(xETest[namesCompare] - xEExp[namesCompare])/pmax(1e-5,abs(xEExp[namesCompare])) < 1e-1))
})

.tmp.f <- function(){
  library(dplyr)
  library(ggplot2)
  res <- cbind(resTest,scen = "Test")
  res <- suppressWarnings(bind_rows(
    cbind(resTest, scen = "Test"), cbind(resExp,scen = "Exp")))
  ggplot(filter(res, time > 1), aes(time, B, color = scen)) + geom_line()
  ggplot(filter(res, time < 500 & time > 0), aes(time, alpha, color = scen)) + geom_point()
  ggplot(filter(res, time > 1), aes(time, L, color = scen)) + geom_line()
  ggplot(filter(res, time > 1), aes(time, R, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, respO, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, ER, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, EL, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, alphaC, color = scen)) + geom_line()
  ggplot(filter(res, time < 5000), aes(time, alphaN, color = scen)) + geom_line()
  ggplot(filter(res, time > 01), aes(time, PhiNB, color = scen, linetype = scen)) + geom_line()
  ggplot(filter(res, time > 01), aes(time, dR, color = scen, linetype = scen)) + geom_line()
  ggplot(filter(res, time > 01), aes(time, dL, color = scen, linetype = scen)) + geom_line()
}

test_that("mass balance with plant N uptake", {
  # mass balance chekced within derivSeam3a
  parmsInit <- within(parms0, {kINPlant = 1000; plantNUpAbs  <- iL/cnIL})
  times <- seq(0,800, length.out = 2)
  times <- seq(0,800, length.out = 101)
  #times <- c(0,148:151)
  #times <- seq(0,2100, by = 2)
  #times <- seq(0,10000, length.out = 101)
  resTest <- as.data.frame(lsoda(
    x0, times, derivSeam3a, parms = parmsInit
    , fBalanceAlpha = balanceAlphaBetweenCNLimitations))
  #
  # N limitation
  times <- seq(0,800, length.out = 8)
  #times <- seq(0,800, length.out = 101)
  resTest <- as.data.frame(lsoda(
    x0Nlim, times, derivSeam3a, parms = parmsInit
    , fBalanceAlpha = balanceAlphaBetweenCNLimitations))
  #
  # from C to N limitation
  x0CNLim <- x0; x0CNLim["IN"] <- 0
  x0CNLimSeam2 <- x0Seam2; x0CNLimSeam2["IN"] <- 0
  times <- seq(0,1200, length.out = 8)
  #times <- seq(0,1200, length.out = 101)
  #times <- c(0,seq(140,220, length.out = 101))
  resTest <- as.data.frame(lsoda( x0CNLim, times, derivSeam3a, parms = parmsInit))
  expect_true(resTest$B[8] > 0)
})



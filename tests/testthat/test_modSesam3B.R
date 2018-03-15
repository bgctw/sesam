#require(testthat)
#test_file("tests/testthat/test_modSesam3B.R")
context("modSesam3B")

parms0 <- list(
  cnB = 7.16
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
  ,epsTvr = 0.3   ##<< carbon use efficiency of microbial tvr (predators respire)
  ,iR = 0         ##<< input modelled explicitely
  ,iL = 300       ##<< g/m2 input per year (half NPP)
  #,plantNUp = 300/70*1/4  ##<< plant N uptake balancing N inputs
  ,plantNUp = 0   ##<< plant N uptake balancing N inputs
  ,useFixedAlloc = FALSE    ##<< set to true to use fixed enzyme allocation (alpha = 0.5)
  ,kIPlant = 10.57 #0.0289652*365         ##<< plant uptake iP I
  ,iB = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iB I
  ,iI = 0         ##<< input of mineral N
  ,l = 0.96 #0.00262647*365       ##<< leaching rate of mineralN l I
  ,nu = 0.9       ##<< microbial N use efficiency
  #, isEnzymeMassFlux = FALSE  ##<< steady state B solution neglects enyzme mass fluxes
  , isEnzymeMassFlux = TRUE  ##<< steady state B solution accounts for enyzme mass fluxes
)
parms0 <- within(parms0,{
  kmR <- kmL <- km
  eps1 <- eps2 <- eps
  cnER <- cnEL <- cnE
  kNR <- kNL <- kN
  kIPlant <- iL / cnIL	# same litter input as plant uptake
  kIPlant <- 0			# no plant uptake
})

parms <- parms0

x0 <- x0Orig <- c( #aE = 0.001*365
  B = 20                     ##<< microbial biomass
  , R = 7000                 ##<< N rich substrate
  , RN = 7000/parms0$cnIR    ##<< N rich substrate N pool
  , L = 200                  ##<< N poor substrate
  , LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  , I =  1                   ##<< inorganic pool
  , alpha = 0.5              ##<< initial community composition
)
x <- x0

# for SEAM model add enzyme state variables
x0Seam3 <- c(x0
             , ER  = 1.5*parms0$km                  ##<< total enzyme pool
             , EL  = 1.5*parms0$km                   ##<< total enzyme pool
             # make sure to have same order as derivative of Seam3
)[c("B","ER","EL","R","RN","L","LN","I","alpha")]
#x0Seam3

x0Nlim <- c( #aE = 0.001*365
  B = 20                     ##<< microbial biomass
  , R = 1000                 ##<< N rich substrate
  , RN = 1000/parms0$cnIR    ##<< N rich substrate N pool
  , L = 200                  ##<< N poor substrate
  , LN = 200/parms0$cnIL     ##<< N poor substrate N pool
  , I =  0                   ##<< inorganic pool
  , alpha = 0.5              ##<< initial community composition
)
x0NlimSeam3 <- c(x0Nlim
                 , ER  = 1.5*parms0$km                  ##<< total enzyme pool
                 , EL  = 1.5*parms0$km                   ##<< total enzyme pool
                 # make sure to have same order as derivative of Seam3
)[c("B","ER","EL","R","RN","L","LN","I","alpha")]

test_that("sesam3BSteadyClim",{
  parms <- parms0
  dL <- 500
  dR <- 200
  alpha <- 0.7
  ans <- sesam3BSteadyClim(dL, dR, alpha, parms = parms)
  B <- ans #max(ans)
  #B <- sesam:::.depr.sesam3BSteadyClim(dL, dR, alpha, parms = parms)
  aE = parms[["aE"]]
  kmkN = parms[["km"]] * parms[["kN"]]
  tau = parms[["tau"]]
  m = parms[["m"]]
  kappaE = parms[["kNB"]]
  eps = parms[["eps"]]
  decL <- dL*(1 - alpha)*aE*B/(kmkN + (1 - alpha)*aE*B)
  decR <- dR*(alpha)*aE*B/(kmkN + (alpha)*aE*B)
  decE <- kappaE*aE*B
  maint <- m*B
  synE <- aE*B/eps
  tvr <- tau*B
  #c(decLN, decRN, immNPot, decLN + decRN + immNPot, tvrBN)
  expect_equal(eps*(decL + decR + decE - synE - maint) - tvr, 0.0, tolerance = 1e-6)
  .tmp.f <- function(){
    Bs <- seq(0,50,length.out = 11)
    ans <- sapply(Bs, function(B){
      decL <- dL*(1 - alpha)*aE*B/(kmkN + (1 - alpha)*aE*B)
      decR <- dR*(alpha)*aE*B/(kmkN + (alpha)*aE*B)
      maint <- m*B
      tvr <- tau*B
      c(eps*(decL + decR + decE - maint), tvr)
      #ansi <- sesam3BSteadyNlim(dLN, dRN, alpha, parms = parms, immNPot = immNPot)
    })
    matplot(Bs, t(ans), type = "l")
  }
  .tmp.f <- function(){
    decFacs <- seq(0.1,0.15,length.out = 51)
    ans <- sapply(decFacs, function(decFac){
      dL2 <- dL*decFac
      dR2 <- dL*decFac
      B <- sesam3BSteadyClim(dL2, dR2, alpha, parms = parms)
    })
    plot(ans ~ decFacs, type = "l")
  }
})

test_that("sesam3BSteadyClim0",{
  parms <- parms0
  # if there is not enough substrate, turnover is larger than decomposition
  # already for very small biomass -> zero biomass
  dL <- 5
  dR <- 2
  alpha <- 0.7
  ans <- sesam3BSteadyClim(dL, dR, alpha, parms = parms)
  expect_equal(ans, 0)
})



test_that("sesam3BSteadyNlim",{
  parms <- parms0
  dLN <- 5 / parms$cnIL
  dRN <- 2 / parms$cnIR
  alpha <- 0.7
  immNPot <- 0.3
  ans <- sesam3BSteadyNlim(dLN, dRN, alpha, parms = parms, immNPot = immNPot)
  B <- ans #Re(ans[2])
  aE = parms[["aE"]]
  kmkN = parms[["km"]] * parms[["kN"]]
  tau = parms[["tau"]]
  betaB <- parms[["cnB"]]
  betaE <- parms[["cnE"]]
  kappaE = parms[["kNB"]]
  nu = parms[["nu"]]  # N efficiency during DON uptake
  decLN <- dLN*(1 - alpha)*aE*B/(kmkN + (1 - alpha)*aE*B)
  decRN <- dRN*(alpha)*aE*B/(kmkN + (alpha)*aE*B)
  synEN <- aE*B/betaE
  decEN <- kappaE*synEN
  tvrBN <- tau*B/betaB
  #c(decLN, decRN, immNPot, nu*(decLN + decRN) + immNPot, tvrBN)
  expect_equal(nu*(decLN + decRN + decEN) + immNPot - synEN - tvrBN, 0.0, tolerance = 1e-6)
  .tmp.f <- function(){
    nPots <- seq(0,10,length.out = 11)
    #nPots <- seq(0,.1,length.out = 11)
    ans <- sapply(nPots, function(immNPot){
      ansi <- sesam3BSteadyNlim(dLN, dRN, alpha, parms = parms, immNPot = immNPot)
    })
    #matplot(nPots, t(Re(ans)), type = "l"); abline(h = 0)
    # only second root seesm reasonably increase with immNPot
    plot(ans ~ nPots, type = "l"); abline(h = 0)
  }
})


test_that("same as sesam3s for fixed substrates", {
  parmsFixedS <- within(parms0,{
    isFixedS <- TRUE
  })
  parmsTvr0 <- within(parms0,{
    isTvrNil <- TRUE
    iR <- 160
  })
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  resTest <- as.data.frame(lsoda( x0[-1], times, derivSesam3B, parms = parmsFixedS))
  resExp <- as.data.frame(lsoda(
    x0, times, derivSesam3s
    , parms = parmsFixedS))
    #, parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-6)
  expect_equal( xETest[2:7], xEExp[3:8], tolerance = 1e-6)
  .tmp.f <- function(){
    derivSesam3B(0, xETest[c(2:7)], within(parmsFixedS, isRecover <- TRUE))
    derivSesam3B(0, xEExp[c(3:8)], within(parmsFixedS, isRecover <- TRUE))
    derivSesam3s(0, xEExp[c(2:8)], within(parmsFixedS, isRecover <- TRUE))
  }
  #
  # N limitation
  #times <- seq(0,2100, length.out = 101)
  resTest <- as.data.frame(lsoda( x0Nlim[-1], times, derivSesam3B, parms = parmsFixedS))
  resExp <- as.data.frame(lsoda(
    x0Nlim, times, derivSesam3s
    , parms = parmsFixedS))
    #, parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-6)
  expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( xETest[2:7], xEExp[3:8], tolerance = 1e-6)
})

test_that("same as sesam3s with substrate feedbacks", {
  parmsInit <- within(parms0, {isFixedI <- TRUE})
  times <- seq(0,800, length.out = 2)
  #times <- seq(0,800, length.out = 101)
  #times <- c(0,148:151)
  #times <- seq(0,2100, by = 2)
  #times <- seq(0,10000, length.out = 101)
  #ans1 <- derivSesam3B(0, x0[-1], within(parmsInit, isRecover <- TRUE) )
  ans1 <- derivSesam3B(0, x0[-1], parmsInit)
  expect_equal(ans1[[2]]["dB"], c(dB = 0), tolerance = 1e-6)
  #
  resTest <- as.data.frame(lsoda( x0[-1], times, derivSesam3B, parms = parmsInit))
  resExp <- as.data.frame(lsoda(
    x0, times, derivSesam3s
    , parms = parmsInit))
  #    , parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1));
  xpESteady <- unlist(head(tail(resTest,2),1))	# the previous before end
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-4)
  expect_equal( xETest[2:7], xEExp[c(3:8)], tolerance = 1e-6)
  #rbind(xEExp[names(xETest)], xETest)
  .tmp.f <- function(){
    derivSesam3B(0, x0[-1], within(parmsInit, isRecover <- TRUE))
  }
  #
  # N limitation
  resTest <- as.data.frame(lsoda( x0Nlim[-1], times, derivSesam3B, parms = parmsInit))
  resExp <- as.data.frame(lsoda(
    x0Nlim, times, derivSesam3s
    , parms = parmsInit))
    #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( xETest[2:7], xEExp[c(3:8)], tolerance = 5e-3)
  #rbind(xEExp[names(xETest)], xETest)
  #
  # from C to N limitation
  x0CNLim <- x0; x0CNLim["I"] <- 0
  x0CNLimSeam3 <- x0Seam3; x0CNLimSeam3["I"] <- 0
  times <- seq(0,1200, length.out = 2)
  #times <- seq(0,1200, length.out = 101)
  #times <- c(0,seq(140,220, length.out = 101))
  resTest <- as.data.frame(lsoda( x0CNLim[-1], times, derivSesam3B, parms = parmsInit))
  resExp <- as.data.frame(lsoda(
    x0CNLim, times, derivSesam3s
    , parms = parmsInit))
    #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-4)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-5)
  expect_equal( xETest[2:7], xEExp[c(3:8)], tolerance = 1e-4)
  rbind(xEExp[names(xETest)], xETest)
})

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





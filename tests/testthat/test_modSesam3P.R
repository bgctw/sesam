#require(testthat)
#test_file("tests/testthat/test_modSesam3P.R")
context("modSesam3P")

parms0 <- list(
  cnB = 7.16
  , cnBW = 10    ##<< C/N ratio of cell walls (that go to R, )
  ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
  #,cnE = 7.16
  ,cnIR = 4.5     ##<< between micr and enzyme signal
  ,cnIL = 30      ##<< N poor substrate
  #,kN = 0.05     ##<< (per day) enzyme turnover
  #,kN = 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
  ,kmN = 0.3*0.01*365  ##<< /yr enzyme turnover 1% turning over each day
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
  #,km = 0.3       ##<< enzyme half-saturation constant
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
)
parms0 <- within(parms0,{
#  kmR <- kmL <- km
 # eps1 <- eps2 <- eps
  #cnER <- cnEL <- cnE
  #kNR <- kNL <- kN
  kIPlant <- iL / cnIL	# same litter input as plant uptake
  kIPlant <- 0			# no plant uptake
  lP <- l       # leaching rate of inorganic P equals that of N
  nuP <- nu     # mineralization of P during decomposiition equals that of N
  kIPPlant <- kIPlant  # plant uptake rate of P equals that of N
  iIP <- l      # assume no P inputs compensate for leaching
})
# for compatibility set C:N:P ratio of cell walls to that of biomass
parms0 <- within(parms0,{
  cnBW <- cnB
  cpBW <- cpB
})

parms <- parms0

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
getX0NoP <- function(x0){x0[setdiff(names(x0),c("RP","LP","IP"))]}
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

test_that("balanceAlphaBetweenElementLimitations",{
  #balanceAlphaBetweenCNLimitations()
  scen <- "C - limited"
  alpha <- c(alphaC = 0.2, alphaN = 0.6, alphaP = 0.8)
  CsynBE <- c(CsynBC = 40, CsynBN = 6000, CSynBP = 6000)
  ce <- c(ccB = NA, cnB = 8, cnP = 8*8)
  eps <- 0.5
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  expect_equal(alphaBalanced$alpha, as.numeric(alpha[1]), tolerance = 1e-2)
  #
  scen <- "N - limited"
  CsynBE <- c(CsynBC = 6000, CsynBN = 40, CSynBP = 6000)
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  expect_equal(alphaBalanced$alpha, as.numeric(alpha[2]), tolerance = 1e-2)
  #
  scen <- "P - limited"
  CsynBE <- c(CsynBC = 6000, CsynBN = 6000, CSynBP = 40)
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  expect_equal(alphaBalanced$alpha, as.numeric(alpha[3]), tolerance = 1e-2)
  #
  scen <- "CN - limited"
  CsynBE <- c(CsynBC = 40, CsynBN = 40, CSynBP = 6000)
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  expect_equal(alphaBalanced$alpha, as.numeric(mean(alpha[1:2])), tolerance = 1e-2)
  #
  scen <- "CN - limited, slightly more N limted"
  CsynBE <- c(CsynBC = 40, CsynBN = 40 - 1e-1, CSynBP = 6000)
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  alphaBalanced
  expect_true(alphaBalanced$alpha > alpha[1] & alphaBalanced$alpha < alpha[2])
  #
  scen <- "NP - limited"
  CsynBE <- c(CsynBC = 6000, CsynBN = 40, CSynBP = 40)
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  expect_equal(alphaBalanced$alpha, as.numeric(mean(alpha[2:3])), tolerance = 1e-2)
  #
  scen <- "NP - limited"
  CsynBE <- c(CsynBC = 6000, CsynBN = 40, CSynBP = 40)
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  expect_equal(alphaBalanced$alpha, as.numeric(mean(alpha[2:3])), tolerance = 1e-2)
  #
  scen <- "equal co-limitation"
  CsynBE <- c(CsynBC = 40, CsynBN = 40, CSynBP = 40)
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  expect_equal(alphaBalanced$alpha, as.numeric(mean(alpha)), tolerance = 1e-2)
  #
  .tmp.f <- function(){
     epsNs <- seq(-0.4,+0.4, length.out = 51)
     alphaBs <- sapply(epsNs, function(epsN){
       CsynBE <- c(CsynBC = 40, CsynBN = 40 - epsN, CSynBP = 6000)
       alphaBalanced <- balanceAlphaBetweenElementLimitations(
         alpha, CsynBE, tauB = 1, delta = 5)$alpha
         #alpha, CsynBE, tauB = 1)#, ce, eps  )
     })
     plot(alphaBs ~ epsNs, type = "l")
     alphaBs10 <- sapply(epsNs, function(epsN){
       CsynBE <- c(CsynBC = 40, CsynBN = 40 - epsN, CSynBP = 6000)
       alphaBalanced <- balanceAlphaBetweenElementLimitations(
         alpha, CsynBE, tauB = 1, delta = 10)$alpha
     })
     lines(alphaBs10 ~ epsNs, lty = "dashed")
     alphaBs20 <- sapply(epsNs, function(epsN){
       CsynBE <- c(CsynBC = 40, CsynBN = 40 - epsN, CSynBP = 6000)
       alphaBalanced <- balanceAlphaBetweenElementLimitations(
         alpha, CsynBE, tauB = 1, delta = 20)$alpha
     })
     lines(alphaBs20 ~ epsNs, lty = "dotted")
  }
  #
  scen <- "strongly C - limited (negative synthesis)"
  # negative carbon balance
  alpha <- c(alphaC = 0.2, alphaN = 0.6, alphaP = 0.8)
  CsynBE <- c(CsynBC = -40, CsynBN = 4, CSynBP = 4)
  ce <- c(ccB = NA, cnB = 8, cnP = 8*8)
  eps <- 0.5
  alphaBalanced <- balanceAlphaBetweenElementLimitations(
    alpha, CsynBE, tauB = 1)#, ce, eps  )
  alphaBalanced
  # think of accounting for negative biomass balance
  expect_equal(alphaBalanced$alpha, as.numeric(alpha[1]), tolerance = 1e-2)
})

test_that("same as sesam3a for fixed substrates", {
  parmsFixedS <- within(parms0,{
    isFixedS <- TRUE
  })
  parmsTvr0 <- within(parms0,{
    isTvrNil <- TRUE
    iR <- 160
  })
  ans0 <- derivSesam3P(0, x0, parms = parmsFixedS)
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  resTest <- as.data.frame(lsoda( x0, times, derivSesam3P, parms = parmsFixedS))
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0), times, derivSesam3a
    , parms = parmsFixedS))
    #, parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-6)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  x0EExp2 <- xETest[2:11]; x0EExp2[names(xEExp[2:8])] <- xEExp[2:8]
  .tmp.f <- function(){
    derivSesam3P(0, xETest[c(2:11)], within(parmsFixedS, isRecover <- TRUE))
    derivSesam3P(0, x0EExp2, within(parmsFixedS, isRecover <- TRUE))
    derivSesam3a(0, xEExp[c(2:8)], within(parmsFixedS, isRecover <- TRUE))
  }
  #
  # N limitation
  #times <- seq(0,2100, length.out = 101)
  times <- seq(0, 8100, length.out = 2)
  resTest <- as.data.frame(lsoda( x0Nlim, times, derivSesam3P, parms = parmsFixedS))
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0Nlim), times, derivSesam3a
    , parms = parmsFixedS))
    #, parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-6)
  expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  #
  # NP col-limitation
  x0Plim <- x0Nlim; x0Plim["IP"] <- 0
  #times <- seq(0,2100, length.out = 101)
  resTest <- as.data.frame(lsoda( x0Plim, times, derivSesam3P, parms = parmsFixedS))
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0Plim), times, derivSesam3a
    , parms = parmsFixedS))
  #, parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_true( xETest["alpha"] < xEExp["alpha"])
  # interestingly this leads to slightly higher biomass under colimitation
  # expect_true( xETest["B"] < xEExp["B"])
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-2)
  #expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
})

test_that("same as sesam3a with substrate feedbacks", {
  parmsInit <- within(parms0, {isFixedI <- TRUE})
  ans0 <- derivSesam3P(0, x0, parms = parmsInit)
  times <- seq(0,800, length.out = 2)
  #times <- seq(0,800, length.out = 101)
  #times <- c(0,148:151)
  #times <- seq(0,2100, by = 2)
  #times <- seq(0,10000, length.out = 101)
  #ans1 <- derivSesam3P(0, x0, within(parmsInit, isRecover <- TRUE) )
  #
  resTest <- as.data.frame(lsoda( x0, times, derivSesam3P, parms = parmsInit))
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0), times, derivSesam3a
    , parms = parmsInit))
  #    , parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1));
  xpESteady <- unlist(head(tail(resTest,2),1))	# the previous before end
  expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-4)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  #rbind(xEExp[names(xETest)], xETest)
  .tmp.f <- function(){
    derivSesam3P(0, x0, within(parmsInit, isRecover <- TRUE))
  }
  #
  # N limitation
  resTest <- as.data.frame(lsoda( x0Nlim, times, derivSesam3P, parms = parmsInit))
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0Nlim), times, derivSesam3a
    , parms = parmsInit))
    #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  #rbind(xEExp[names(xETest)], xETest)
  #
  # from C to N limitation
  x0CNLim <- x0; x0CNLim["I"] <- 0
  times <- seq(0,1200, length.out = 2)
  #times <- seq(0,1200, length.out = 101)
  #times <- c(0,seq(140,220, length.out = 101))
  resTest <- as.data.frame(lsoda( x0CNLim, times, derivSesam3P, parms = parmsInit))
  resExp <- as.data.frame(lsoda(
    getX0NoP(x0CNLim), times, derivSesam3a
    , parms = parmsInit))
    #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xETest <- unlist(tail(resTest,1))
  xEExp <- unlist(tail(resExp,1))
  expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-4)
  expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-5)
  expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
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





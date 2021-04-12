.tmp.f <- function(){
  require(testthat)
}
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
  ,kINPlant = 10.57 #0.0289652*365         ##<< plant uptake iP IN
  ,iBN = 0.38 * 10.57 #0.0110068*365   ##<< immobilization flux iBN IN
  ,iIN = 0         ##<< input of mineral N
  ,lN = 0.96 #0.00262647*365       ##<< leaching rate of mineralN lN IN
  , nuN = 0.9       ##<< microbial N use efficiency
  #, isEnzymeMassFlux = FALSE  ##<< steady state B solution neglects enyzme mass fluxes
  , isEnzymeMassFlux = TRUE  ##<< steady state B solution accounts for enyzme mass fluxes
  , nuP = 0.3      ##<< microbial uptake of depolymerized P, (1-nuP) is mineralized
  , cpE = 50
  , cpB = 40
  , cpBW = 50
  , cpIR = 40
  , cpIL = 40*3
  , iBP = 0.38 * 10.57 # start with same as N
  , kSP = 1/30 # 1/x years
)
parms0 <- within(parms0,{
#  kmR <- kmL <- km
 # eps1 <- eps2 <- eps
  #cnER <- cnEL <- cnE
  #kNR <- kNL <- kN
  lP <- lN       # leaching rate of inorganic P equals that of N
  nuP <- nuN     # mineralization of P during decomposiition equals that of N
  kIPPlant <- kINPlant  # plant uptake rate of P equals that of N
  iIP <- lN      # assume no P inputs compensate for leaching
  plantNUpAbs <- iL / cnIL	# same litter input as plant uptake
  kINPlant <- plantNUpAbs <- 0			# no plant uptake
  kLP <- kRP <- kSP
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
  , IN =  1                   ##<< inorganic N pool
  , IP =  1                  ##<< inorganic P pool
  , alphaL = 0.4              ##<< initial community composition
  , alphaR = 0.5              ##<< initial community composition
  , alphaRP = 0.1              ##<< initial community composition
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
  , IN =  0                   ##<< inorganic N pool
  , IP =  10                  ##<< inorganic P pool
  , alpha = 0.5              ##<< initial community composition
)
# fixed substrate ---------------------------------------------------------
test_that("same as sesam3a for fixed substrates", {
  parmsFixedS <- within(parms0,{
    isFixedS <- TRUE
  })
  parmsTvr0 <- within(parms0,{
    isTvrNil <- TRUE
    iR <- 160
  })
  ans0 <- derivSesam3P(0, x0, parms = parms0)
  expect_true(length(ans0) > 0)
  # ans0 <- derivSesam3P(0, x0, parms = parmsFixedS)
  # times <- seq(0, 2100, length.out = 2)
  # #times <- seq(0,2100, length.out = 101)
  # resTest <- as.data.frame(lsoda( x0, times, derivSesam3P, parms = parmsFixedS))
  # resExp <- as.data.frame(lsoda(
  #   getX0NoP(x0), times, derivSesam3a
  #   , parms = parmsFixedS))
  #   #, parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  # xETest <- unlist(tail(resTest,1))
  # xEExp <- unlist(tail(resExp,1))
  # expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-6)
  # expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  # x0EExp2 <- xETest[2:11]; x0EExp2[names(xEExp[2:8])] <- xEExp[2:8]
  # .tmp.f <- function(){
  #   derivSesam3P(0, xETest[c(2:11)], within(parmsFixedS, isRecover <- TRUE))
  #   derivSesam3P(0, x0EExp2, within(parmsFixedS, isRecover <- TRUE))
  #   derivSesam3a(0, xEExp[c(2:8)], within(parmsFixedS, isRecover <- TRUE))
  # }
  # #
  # # N limitation
  # #times <- seq(0,2100, length.out = 101)
  # times <- seq(0, 8100, length.out = 2)
  # resTest <- as.data.frame(lsoda( x0Nlim, times, derivSesam3P, parms = parmsFixedS))
  # resExp <- as.data.frame(lsoda(
  #   getX0NoP(x0Nlim), times, derivSesam3a
  #   , parms = parmsFixedS))
  #   #, parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  # xETest <- unlist(tail(resTest,1))
  # xEExp <- unlist(tail(resExp,1))
  # expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-6)
  # expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
  # expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
  # expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
  # #
  # # NP col-limitation
  # x0Plim <- x0Nlim; x0Plim["IP"] <- 0
  # #times <- seq(0,2100, length.out = 101)
  # resTest <- as.data.frame(lsoda( x0Plim, times, derivSesam3P, parms = parmsFixedS))
  # resExp <- as.data.frame(lsoda(
  #   getX0NoP(x0Plim), times, derivSesam3a
  #   , parms = parmsFixedS))
  # #, parms = within(parmsFixedS, isEnzymeMassFlux <- FALSE)))
  # xETest <- unlist(tail(resTest,1))
  # xEExp <- unlist(tail(resExp,1))
  # expect_true( xETest["alpha"] < xEExp["alpha"])
  # # interestingly this leads to slightly higher biomass under colimitation
  # # expect_true( xETest["B"] < xEExp["B"])
  # expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-2)
  # #expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
})

# # substrate feedbacks ---------------------------------------------------------
#
#
# test_that("same as sesam3a with substrate feedbacks", {
#   parmsInit <- within(parms0, {isFixedI <- TRUE})
#   ans0 <- derivSesam3P(0, x0, parms = parmsInit)
#   times <- seq(0,800, length.out = 2)
#   #times <- seq(0,800, length.out = 101)
#   #times <- c(0,148:151)
#   #times <- seq(0,2100, by = 2)
#   #times <- seq(0,10000, length.out = 101)
#   #ans1 <- derivSesam3P(0, x0, within(parmsInit, isRecover <- TRUE) )
#   #
#   resTest <- as.data.frame(lsoda( x0, times, derivSesam3P, parms = parmsInit))
#   resExp <- as.data.frame(lsoda(
#     getX0NoP(x0), times, derivSesam3a
#     , parms = parmsInit))
#   #    , parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
#   xETest <- unlist(tail(resTest,1))
#   xEExp <- unlist(tail(resExp,1));
#   xpESteady <- unlist(head(tail(resTest,2),1))	# the previous before end
#   expect_equal( xETest["alphaC"], xEExp["alphaC"], tolerance = 1e-4)
#   expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
#   #rbind(xEExp[names(xETest)], xETest)
#   .tmp.f <- function(){
#     derivSesam3P(0, x0, within(parmsInit, isRecover <- TRUE))
#   }
#   #
#   # N limitation
#   resTest <- as.data.frame(lsoda( x0Nlim, times, derivSesam3P, parms = parmsInit))
#   resExp <- as.data.frame(lsoda(
#     getX0NoP(x0Nlim), times, derivSesam3a
#     , parms = parmsInit))
#     #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
#   xETest <- unlist(tail(resTest,1))
#   xEExp <- unlist(tail(resExp,1))
#   expect_equal( xETest["alpha"], xEExp["alpha"], tolerance = 1e-6)
#   expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-6)
#   expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
#   #rbind(xEExp[names(xETest)], xETest)
#   #
#   # from C to N limitation
#   x0CNLim <- x0; x0CNLim["IN"] <- 0
#   times <- seq(0,1200, length.out = 2)
#   #times <- seq(0,1200, length.out = 101)
#   #times <- c(0,seq(140,220, length.out = 101))
#   resTest <- as.data.frame(lsoda( x0CNLim, times, derivSesam3P, parms = parmsInit))
#   resExp <- as.data.frame(lsoda(
#     getX0NoP(x0CNLim), times, derivSesam3a
#     , parms = parmsInit))
#     #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
#   xETest <- unlist(tail(resTest,1))
#   xEExp <- unlist(tail(resExp,1))
#   expect_equal( xETest["B"], xEExp["B"], tolerance = 1e-4)
#   expect_equal( xETest["alphaN"], xEExp["alphaN"], tolerance = 1e-5)
#   expect_equal( getX0NoP(xETest[2:11]), xEExp[2:8], tolerance = 1e-6)
#   rbind(xEExp[names(xETest)], xETest)
# })
#
# .tmp.f <- function(){
#   # plots of results
#   library(dplyr)
#   library(ggplot2)
#   res <- suppressWarnings(bind_rows(
#     cbind(resTest, scen = "Test"), cbind(resExp, scen = "Exp")))
#   ggplot(filter(res, time > 1), aes(time, B, color = scen)) + geom_line(alpha = 0.5)
#   ggplot(filter(res, time > 0), aes(time, alpha, color = scen)) + geom_point(alpha = 0.5)
#   ggplot(filter(res, time < 500 & time > 0), aes(time, alpha, color = scen)) + geom_point()
#   ggplot(filter(res, time < 500), aes(time, B, color = scen)) + geom_line()
#   ggplot(filter(res, time >= 0), aes(time, L, color = scen)) + geom_line()
#   ggplot(filter(res, time >= 0), aes(time, R, color = scen)) + geom_line()
#   ggplot(filter(res, time < 500), aes(time, respO, color = scen)) + geom_line()
#   ggplot(filter(res, time > 10 & time < 500), aes(time, ER, color = scen)) + geom_line()
#   ggplot(filter(res, time > 10 & time < 500), aes(time, EL, color = scen)) + geom_line()
#   ggplot(filter(res, time < 500), aes(time, alphaC, color = scen)) + geom_line()
#   ggplot(filter(res, time < 5000), aes(time, alphaN, color = scen)) + geom_line()
#   ggplot(filter(res, time > 01), aes(time, PhiNB, color = scen, linetype = scen)) + geom_line()
# }
#
#
#
#

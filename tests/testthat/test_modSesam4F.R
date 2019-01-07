#require(testthat)
#test_file("tests/testthat/test_modSesam4F.R")
context("modSesam4F")

parms0 <- list(
  cnB = 7.16
  , cnBW = 10    ##<< C/N ratio of cell walls (that go to R, )
  ,cnE = 3.1     # Sterner02: Protein (Fig. 2.2.), high N investment (low P)
  #,cnE = 7.16
  ,cnIR = 4.5     ##<< between micr and enzyme signal
  ,cnIL = 30      ##<< N poor substrate
  #,kN = 0.05     ##<< (per day) enzyme turnover
  ##,kN = 0.01*365  ##<< /yr enzyme turnover 1% turning over each day
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
  ## ,km = 0.3       ##<< enzyme half-saturation constant
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
  relIL <- list( C = c(SOM = 0, amend = 1))
  relIL$N <- relIL$P <- relIL$C
  relIR <- list(C = c(SOM = 0, amend = 1))
  relIR$N <- relIR$P <- relIR$C
  relII <- c(SOM = 0, amend = 1)
  relIIP <- c(SOM = 0, amend = 1)
})
parms0 <- within(parms0,{
  #kmR <- kmL <- km
  #eps1 <- eps2 <- eps
  #cnER <- cnEL <- cnE
  #kNR <- kNL <- kN
  kIPlant <- iL / cnIL	# same litter input as plant uptake
  kIPlant <- 0			# no plant uptake
  lP <- l       # leaching rate of inorganic P equals that of N
  nuP <- nu     # mineralization of P during decomposiition equals that of N
  kIPPlant <- kIPlant  # plant uptake rate of P equals that of N
  iIP <- l      # assume no P inputs compensate for leaching
})

x0 <- x0Orig <- c( #aE = 0.001*365
  BC_SOM = 20                     ##<< microbial biomass
  , RC_SOM = 7000                 ##<< N rich substrate
  , LC_SOM = 200                  ##<< N poor substrate
  , resp_SOM = 0
  , BN_SOM = 20/parms0$cnB
  , RN_SOM = 7000/parms0$cnIR    ##<< N rich substrate N pool
  , LN_SOM = 200/parms0$cnIL     ##<< N poor substrate N pool
  , I_SOM =  1                   ##<< inorganic N pool
  , leachN_SOM = 0
  , BP_SOM = 20/parms0$cpB
  , RP_SOM = 7000/parms0$cpIR
  , LP_SOM = 200/parms0$cpIL     ##<< N poor substrate P pool
  , IP_SOM =  1                  ##<< inorganic P pool
  , leachP_SOM = 0
  , BC_amend = 0                     ##<< microbial biomass
  , RC_amend = 0                 ##<< N rich substrate
  , LC_amend = 0
  , resp_amend = 0
  , BN_amend = 0
  , RN_amend = 0
  , LN_amend = 0
  , I_amend =  0
  , leachN_amend = 0
  , BP_amend = 0
  , RP_amend = 0
  , LP_amend = 0
  , IP_amend =  0
  , leachP_amend = 0
  , alpha = 0.5              ##<< initial community composition
)
# TODO: need to resort so that all fractions of one pool are toegehter
x0F <- parms0$multiPoolFractions$setX(parms0$multiPoolFractions, x0)
x0 <- x0F$stateVec(x0F)
getX0Sum <- function(x0, parms = parms0){
  x <- parms$multiPoolFractions
  x <- x$setX(x, x0)
  x$tot
}
x0Sum <- getX0Sum(x0)
test_that("equal initial sum of pools",{
  expect_equivalent(x0Sum["alpha"], x0["alpha"])
  expect_equivalent(x0Sum[1:14], x0Orig[1:14])
})


x0Nlim <- x0
.RC <- 1000
x0Nlim["RC_SOM"] = .RC
x0Nlim["RN_SOM"] = .RC/parms0$cnIR    ##<< N rich substrate N pool
x0Nlim["RP_SOM"] = .RC/parms0$cpIR
x0NlimSum <- getX0Sum(x0Nlim)


.tmp.f <- function(){
  parms1 <- within(parms0, {isFixedI <- TRUE; l <- 0; iL <- 0})
  ans1 <- derivSesam4F(0, x0, parms1)
  x0A <- x0; x0A["LC_amend"] <- x0["LC_SOM"]; x0A["LC_SOM"] <- x0["LC_amend"]
  x0A["LN_amend"] <- x0["LN_SOM"]; x0A["LN_SOM"] <- x0["LN_amend"]
  x0A["LP_amend"] <- x0["LP_SOM"]; x0A["LP_SOM"] <- x0["LP_amend"]
  ans1 <- derivSesam4F(0, x0A, parms1)
  times <- seq(0,500, length.out = 101)
  resF <- res1 <- as.data.frame(lsoda(
    #x0
    x0A
    , times, derivSesam4F
    ,parms = parms1
    #, parms = within(parms0, {isFixedI <- TRUE})
    #, parms = within(parms0, {isFixedL <- TRUE; iI <- 1}) # works since no immobiization
    #, parms = within(parms0, {isFixedI <- TRUE; isFixedL <- TRUE}) #works
    #, parms = within(parms0, {isFixedI <- TRUE; isFixedR <- TRUE; isFixedL <- TRUE}) # works
    ))
  resF <- resF %>% mutate(scen = "tmpf1")
  #
  res <- sumMultiPoolFractions(x0F, resF)
}

namesF = c(       ##<< names in totals of sesam4F
  "BC", "RC", "LC", "RN", "LN", "I", "RP", "LP", "IP", "alpha")
namesA = c(       ##<< names in resA of the same order as namesF
  "BC", "RC", "LC", "RN", "LN", "I", "RP", "LP", "IP", "alpha")


testFState <- function(
  ### test numeric equality of results of sesam4a and sesam4F
  resF                ##<< data.frame or vector of state vectors of sesam4F
  , resA              ##<< data.frame or vector of state vectors of sesam4a
  , parms             ##<< list with entries cnB, cpB, and multiPoolFractions
  , tolerance = 1e-4  ##<< tolerance of differences
){
  resFs <- sumMultiPoolFractions(parms$multiPoolFractions, resF)
  expect_equal( resFs$BC, parms$cnB*resFs$BN, tolerance = tolerance)
  expect_equal( resFs$BC, parms$cpB*resFs$BP, tolerance = tolerance)
  if (!(is.matrix(resA) || is.data.frame(resA))) resA <- as.data.frame(t(resA))
  expect_equal( resFs[namesF], structure(resA[namesA], names = namesF)
                , tolerance = tolerance )
}

#parmsInit <- parms0
testParmsScen <- function(parmsInit){
  ans0 <- derivSesam4F(0, x0, parms = parmsInit)
  xF <- x0F$setX(x0F, x0)
  ans0E <- derivSesam4a(0, xF$tot[namesF], parms = parmsInit)
  testFState(ans0[[1]], ans0E[[1]], parmsInit)
  #
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  #times <- seq(0,0.2, length.out = 101)
  xF <- x0F$setX(x0F, x0)
  resExp <- as.data.frame(lsoda(
    xF$tot[namesF], times, derivSesam4a
    , parms = parmsInit))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0, times, derivSesam4F
    , parms = parmsInit))
  xETest <- unlist(tail(resTest,1))
  xFE <- x0F$setX(x0F, xETest)
  testFState( xETest, xEExp, parmsInit)
  #
  # N limitation
  #times <- seq(0,2100, length.out = 101)
  xF <- xF$setX(xF, x0Nlim)
  resExp <- as.data.frame(lsoda(
    xF$tot[namesF], times, derivSesam4a
    , parms = parmsInit))
  #, parms = within(parmsInit, isEnzymeMassFlux <- FALSE)))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Nlim, times, derivSesam4F
    , parms = parmsInit))
  xETest <- unlist(tail(resTest,1))
  testFState( xETest, xEExp, parmsInit)
  #
  # NP col-limitation
  x0Plim <- x0Nlim; x0Plim["IP_SOM"] <- x0Plim["IP_amend"] <- 1e-12
  xF <- xF$setX(xF, x0Plim)
  parmsInitPlim <- within(parmsInit, cpIL <- 160)
  #times <- seq(0,2100, length.out = 101)
  resExp <- as.data.frame(lsoda(
    xF$tot[namesF], times, derivSesam4a
    , parms = parmsInit))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Plim, times, derivSesam4F
    , parms = parmsInit))
  xETest <- unlist(tail(resTest,1))
  testFState( xETest, xEExp, parmsInit)
  #
  # Microbial starvation
  x0Starv <- setMultiPoolFractionsElements(xF, x0, "B", 160)
  .tmp.f <- function() {
    #track differences between N and P by setting all P to N dynamics
    x0Starv <- setMultiPoolFractionsPool(xF, x0Starv, "LP", x0Starv["LN_SOM"])
    x0Starv <- setMultiPoolFractionsPool(xF, x0Starv, "RP", x0Starv["RN_SOM"])
    x0Starv <- setMultiPoolFractionsPool(xF, x0Starv, "BP", x0Starv["BN_SOM"])
    x0Starv <- setMultiPoolFractionsPool(xF, x0Starv, "IP", x0Starv["I_SOM"])
    parmsInit <- within(parmsInit, {cpB <- cnB; cpE <- cnE; cpBW <- cnBW; nuP <- nu; iIP <- iI; cpIL <- cnIL})
  }
  xF <- xF$setX(xF, x0Starv)
  times <- seq(0, 2100, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  #times <- seq(0,2, length.out = 101)
  resExp <- as.data.frame(lsoda(
    xF$tot[namesF], times, derivSesam4a
    , parms = parmsInit))
  xEExp <- unlist(tail(resExp,1))
  resTest <- as.data.frame(lsoda(
    x0Starv, times, derivSesam4F
    , parms = parmsInit))
  xETest <- unlist(tail(resTest,1))
  testFState( xETest, xEExp, parmsInit)
}

.tmp.f <- function(){
  # debug error by running both models from xETest
  ans0 <- derivSesam4F(0, xETest, parms = parmsInit)
  xF <- x0F$setX(x0F, xETest)
  ans0E <- derivSesam4a(0, xF$tot[namesF], parms = parmsInit)
  testFState(ans0[[1]], ans0E[[1]], parmsInit)

  sort( abs(sumMultiPoolFractions(xF,ans0[[1]])[namesF] -
  ans0E[[1]][namesF]), decreasing = TRUE)

  x0B <- x0
  x0 <- xETest[2:24]
}

.tmp.f <- function(){
  # plots of results
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  # res <- suppressWarnings(bind_rows(
  #   cbind(resTest, scen = "Test"), cbind(resTest2, scen = "Test2"), cbind(resExp, scen = "Exp")))
  res <- suppressWarnings(bind_rows(
    cbind(sumMultiPoolFractions(xF, resTest), scen = "Test"), cbind(resExp, scen = "Exp")))
  ggplot(filter(res, time > 1), aes(time, BC, color = scen)) + geom_line(alpha = 0.5)
  ggplot(filter(res, time > 0), aes(time, alpha, color = scen)) + geom_point(alpha = 0.5)
  ggplot(filter(res, time < 500 & time > 0), aes(time, alpha, color = scen)) + geom_point()
  ggplot(filter(res, time < 500), aes(time, BC, color = scen)) + geom_line()
  ggplot(filter(res, time <= 2), aes(time, BC, color = scen, linetype = scen)) + geom_line()
  ggplot(filter(res, time <= 2), aes(time, starvB, color = scen, linetype = scen)) + geom_line()
  ggplot(filter(res, time >= 0), aes(time, LC, color = scen)) + geom_line()
  ggplot(filter(res, time >= 0), aes(time, RC, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, respO, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, ER, color = scen)) + geom_line()
  ggplot(filter(res, time > 10 & time < 500), aes(time, EL, color = scen)) + geom_line()
  ggplot(filter(res, time < 500), aes(time, alphaC, color = scen)) + geom_line()
  ggplot(filter(res, time < 5000), aes(time, alphaN, color = scen)) + geom_line()
  ggplot(filter(res, time > 01), aes(time, PhiB, color = scen, linetype = scen)) + geom_line()
}

.tmp.f <- function(){
  resTs <- resTest %>%
    select(time, contains("_")) %>%
    gather("var","value", -time) %>%
    separate(var,c("pool","frac")) %>%
    spread(pool, value)
  ggplot(filter(resTs, time <= 2), aes(time, BN, color = frac)) + geom_line(alpha = 0.5)
  ggplot(filter(resTs, time <= 2), aes(time, LN, color = frac)) + geom_line(alpha = 0.5)
  #ggplot(filter(resTs, time <= 2), aes(time, RN, color = frac)) + geom_line(alpha = 0.5)

}

test_that("same as sesam4a for fixed substrates", {
  parmsInit <- parmsFixedS <- within(parms0,{
    isFixedS <- TRUE
  })
  testParmsScen(parmsInit)
})


test_that("mass balance also with feedback to DOM", {
  parmsInit <- within(parms0, {isFixedI <- isFixedIP <- TRUE})
  # no turnover to DOM impossible with cnBW != cnB
  expect_error(
    ans0 <- derivSesam4F(0, x0, parms = within(parmsInit, {cW <- 1; B0 <- 0}))
  )
  # with no turnover by predation
  ans0 <- derivSesam4F(0, x0, parms = within(parmsInit, {B0 <- 1e20}))
  # with only turnover by predation
  ans0 <- derivSesam4F(0, x0, parms = within(parmsInit, {tau <- 0; B0 <- 0}))
  # similar to sesam3 but with tau increasing
  ans0 <- derivSesam4F(0, x0, parms = within(parmsInit, {
    cnBW = cnB; cpBW = cpB; cW <- 1; tau <- 0; B0 <- 0}))
})

test_that("same as sesam4a with substrate feedbacks", {
  parmsInit <- within(parms0, {isFixedI <- isFixedIP <- TRUE})
  testParmsScen(parmsInit)
})

test_that("mass balance in output pools", {
  parmsInit <- within(parms0,{})
  endTime <- 5
  times <- seq(0, endTime, length.out = 2)
  #times <- seq(0,2100, length.out = 101)
  #times <- seq(0,0.2, length.out = 101)
  x0F <- x0F$setX(x0F, x0)
  resTest <- as.data.frame(lsoda(
    x0, times, derivSesam4F
    , parms = parmsInit))
  xETest <- unlist(tail(resTest,1))
  # totals
  xF <- x0F$setX(x0F, xETest)
  expect_equal(
    as.numeric(sum(xF$tot[c("BC","RC","LC","resp")]))
    , as.numeric(sum(x0F$tot[c("BC","RC","LC","resp")]) + xETest["time"] * parmsInit$iL)
    , tolerance = 1e-8)
  expect_equal(
    as.numeric(sum(xF$tot[c("BN","RN","LN","I","leachN")]))
    , as.numeric(sum(x0F$tot[c("BN","RN","LN","I","leachN")]) + xETest["time"]*(
      parmsInit$iL/parmsInit$cnIL + parmsInit$iI))
    , tolerance = 1e-8)
  expect_equal(
    as.numeric(sum(xF$tot[c("BP","RP","LP","IP","leachP")]))
    , as.numeric(sum(x0F$tot[c("BP","RP","LP","IP","leachP")]) + xETest["time"]*(
      parmsInit$iL/parmsInit$cpIL + parmsInit$iIP))
    , tolerance = 1e-8)
  # fractions
  expect_equal(
    as.numeric(rowSums(elementMultiPoolFractions(xF,"C", c("BC","RC","LC","resp"))))
    , as.numeric(rowSums(elementMultiPoolFractions(x0F,"C", c("BC","RC","LC","resp"))) +
      xETest["time"]*(parmsInit$iL*parmsInit$relIL$C))
    , tolerance = 1e-8)
  expect_equal(
    as.numeric(rowSums(elementMultiPoolFractions(xF,"N", c("BN","RN","LN","I","leachN"))))
    , as.numeric(rowSums(elementMultiPoolFractions(x0F,"N", c("BN","RN","LN","I","leachN"))) +
                   xETest["time"]*(parmsInit$iL/parmsInit$cnIL*parmsInit$relIL$N +
                                     parmsInit$iI*parmsInit$relII))
    , tolerance = 1e-8)
  expect_equal(
    as.numeric(rowSums(elementMultiPoolFractions(xF,"P", c("BP","RP","LP","IP","leachP"))))
    , as.numeric(rowSums(elementMultiPoolFractions(x0F,"P", c("BP","RP","LP","IP","leachP"))) +
                   xETest["time"]*(parmsInit$iL/parmsInit$cpIL*parmsInit$relIL$P +
                                     parmsInit$iIP*parmsInit$relIIP))
    , tolerance = 1e-8)
})


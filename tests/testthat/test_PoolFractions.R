#require(testthat)
#test_file("tests/testthat/test_PoolFractions.R")
context("MultiPoolFractions")

test_that("MultiPoolFractions amend through all fractions", {
  elements <- c("C","N","P")
  units <- structure(lapply(elements, function(el) c(SOM = 1, amend = 1))
                     , names = elements)
  x <- createMultiPoolFractions(units, setX = createSesam4setX(units))
  expect_equal( x$units, units)
  expect_equal( names(x$tot), c("BC","RC","LC","BN","RN","LN","I","BP","RP","LP","IP","alpha"))
  expect_true( all( c("RC_SOM","RC_amend","RP_SOM","alpha") %in% names(x$stateVec(x)) ))
  #.self <- x  # rm(.self)  # poolName <- "RP"
  x0 <- structure(seq_along(x$stateVec(x)), names = names(x$stateVec(x)))
  x <- x$setX(x,x0)
  expect_equivalent( x$frac[["alpha"]], x0[c("alpha")])
  expect_equivalent( x$frac[["RP"]], x0[c("RP_SOM","RP_amend")])
  expect_equivalent( x$frac[["RC"]], x0[c("RC_SOM","RC_amend")])
  expect_equivalent( x$tot["alpha"], x0[c("alpha")])
  expect_equivalent( x$tot["RC"], sum(x0[c("RC_SOM","RC_amend")]*units$C))
  expect_equivalent( x$tot["RN"], sum(x0[c("RN_SOM","RN_amend")]*units$N))
  expect_equivalent( x$rel[["alpha"]], 1)
  expect_equivalent( x$rel[["RC"]], x$frac[["RC"]]/x$tot["RC"])
  expect_equivalent( x$rel[["RN"]], x$frac[["RN"]]/x$tot["RN"])
  # reli <- x$rel[["BC"]]
  expect_true( all(sapply(x$rel[c("BC","RC","LC")], function(reli) sum(reli * units$C)) - 1 < 1e-12) )
  expect_true( all(sapply(x$rel[c("RN","LN","I")], function(reli) sum(reli * units$N)) - 1 < 1e-12) )
  expect_true( all(sapply(x$rel[c("RP","LP","IP")], function(reli) sum(reli * units$P)) - 1 < 1e-12) )
  expect_true( all(unlist(x$rel[["alpha"]]) - 1 < 1e-12) )
})

test_that("MultiPoolFractions 13C 14N, noPIso", {
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppTrillion
    , N = c(N14 = 1, N15 = 0.01) # 15N in percent
  )
  x <- createMultiPoolFractions(units, setX = createSesam4CNsetX(units))
  expect_equal( names(x$tot), c("BC","RC","LC","BN","RN","LN","I","BP","RP","LP","IP","alpha"))
  expect_true( all( c("RC_C12","RC_C13","RN_N14","alpha","RP") %in% names(x$stateVec(x)) ))
  #.self <- x  # rm(.self)  # poolName <- "RP"
  x0 <- structure(seq_along(x$stateVec(x)), names = names(x$stateVec(x)))
  x <- x$setX(x,x0)
  expect_equivalent( x$frac[["alpha"]], x0[c("alpha")])
  expect_equivalent( x$frac[["RP"]], x0[c("RP")])
  expect_equivalent( x$frac[["RC"]], x0[c("RC_C12","RC_C13","RC_C14")])
  expect_equivalent( x$tot["alpha"], x0[c("alpha")])
  expect_equivalent( x$tot["RC"], sum(x0[c("RC_C12","RC_C13","RC_C14")]*units$C))
  expect_equivalent( x$tot["RN"], sum(x0[c("RN_N14","RN_N15")]*units$N))
  expect_equivalent( x$rel[["alpha"]], 1)
  expect_equivalent( x$rel[["RC"]], x$frac[["RC"]]/x$tot["RC"])
  expect_equivalent( x$rel[["RN"]], x$frac[["RN"]]/x$tot["RN"])
  expect_true( all(sapply(x$rel[c("BC","RC","LC")], function(reli) sum(reli * units$C)) - 1 < 1e-12) )
  expect_true( all(sapply(x$rel[c("RN","LN","I")], function(reli) sum(reli * units$N)) - 1 < 1e-12) )
  expect_true( all(sapply(x$rel[c("RP","LP","IP")], function(reli) sum(reli * units$P)) - 1 < 1e-12) )
  expect_true( all(unlist(x$rel[["alpha"]]) - 1 < 1e-12) )
})

test_that("sumMultiPoolFractions", {
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppTrillion
    , N = c(N14 = 1, N15 = 0.01) # 15N in percent
  )
  x <- createMultiPoolFractions(units, setX = createSesam4CNsetX(units))
  x0 <- structure(seq_along(x$stateVec(x)), names = names(x$stateVec(x)))
  ds <- do.call(rbind, lapply(1:5, function(time){ data.frame(
    time = time, t(x0), note = "bla")}))
  #.self <- x
  dsSum <- sumMultiPoolFractions(x, ds)
  expect_equal( nrow(dsSum), nrow(ds) )
  expect_true( all(names(x$tot) %in% names(dsSum)) )
  expect_equal( dsSum$RC, rowSums(ds[,c("RC_C12","RC_C13","RC_C14")]*matrix(
    x$units$C, nrow = nrow(ds), ncol = length(x$units$C), byrow = TRUE)))
})

test_that("sumMultiPoolFractionsVars", {
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppTrillion
    , N = c(N14 = 1, N15 = 0.01) # 15N in percent
  )
  x <- createMultiPoolFractions(units, setX = createSesam4CNsetX(units))
  x0 <- structure(seq_along(x$stateVec(x)), names = names(x$stateVec(x)))
  ds <- do.call(rbind, lapply(1:5, function(time){ data.frame(
    time = time, t(x0), decNL_N14 = 1:5, decNL_N15 = 6:10)}))
  #.self <- x
  dsSum <- sumMultiPoolFractionsVars(x, ds, list(N = c("decNL")))
  expect_equal( nrow(dsSum), nrow(ds) )
  expect_equal( dsSum$decNL, rowSums(ds[,c("decNL_N14","decNL_N15")]*matrix(
    x$units$N, nrow = nrow(ds), ncol = length(x$units$N), byrow = TRUE)))
})

test_that("setMultiPoolFractionsElements", {
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppTrillion
    , N = c(N14 = 1, N15 = 0.01) # 15N in percent
    , P = c(tot = 1) # P only one pool
  )
  x <- createMultiPoolFractions(units, setX = createSesam4setX(units))
  #.self <- x
  x0 <- structure(seq_along(x$stateVec(x)), names = names(x$stateVec(x)))
  C12 <- 100; cnB <- 8.3; cpB <- 47.3
  x1 <- setMultiPoolFractionsElements(x, x0, "B", C12
                            , c(N = cnB,  P = cpB)
                            , list(C = c(C13 = 2, C14 = 4), N = c(N15 = 8)))
  x1F <- x$setX(x, x1)
  expect_equal(as.vector(x1["BC_C12"]), C12)
  expect_equal(as.vector(x1["BC_C13"]), C12*2)
  expect_equal(as.vector(x1["BC_C14"]), C12*4)
  expect_equal(as.vector(x1F$tot["BC"]), sum(C12*c(1,2,4)*units$C))
  # other fractions
  Ntot <- x1F$tot["BN"]
  expect_equal(as.vector(Ntot*cnB), as.vector(x1F$tot["BC"])) # cn-ratio
  N14 <- x1["BN_N14"]
  expect_equal(as.vector(x1["BN_N15"]), as.vector(N14*8))
  # scalar elements
  expect_equal(as.vector(x1["BP_tot"]*cpB), as.vector(x1F$tot["BC"]))
  #
  # test setting without specifying fractions
  C12b <- 200
  x1b <- setMultiPoolFractionsElements(x, x1, "B", C12b
                                   , c(N = cnB*2,  P = cpB*2))
  x1bF <- x$setX(x, x1b)
  expect_equal(as.vector(x1bF$tot["BN"]*cnB*2), as.vector(x1bF$tot["BC"])) # cn-ratio
  expect_equal(as.vector(x1bF$tot["BP"]*cpB*2), as.vector(x1bF$tot["BC"])) # cp-ratio
  expect_equal(as.vector(x1bF$frac[["BC"]][-1]/x1bF$frac[["BC"]][1]), c(2,4))
  #
  # test setting without specifying elemental ratios
  C12c <- 400
  x1c <- setMultiPoolFractionsElements(x, x1, "B", C12c)
  x1cF <- x$setX(x, x1c)
  expect_equal(as.vector(x1cF$tot["BN"]*cnB), as.vector(x1cF$tot["BC"])) # cn-ratio
  expect_equal(as.vector(x1cF$tot["BP"]*cpB), as.vector(x1cF$tot["BC"])) # cp-ratio
  expect_equal(as.vector(x1cF$frac[["BC"]][-1]/x1cF$frac[["BC"]][1]), c(2,4))
})




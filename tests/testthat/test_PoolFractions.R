#require(testthat)
#test_file("tests/testthat/test_PoolFractions.R")
context("MultiPoolFractions")

test_that("MultiPoolFractions amend through all fractions", {
  elements <- c("C","N","P")
  units <- structure(lapply(elements, function(el) c(SOM = 1, amend = 1))
                     , names = elements)
  x <- createMultiPoolFractions(units, setX = createSesam4setX(units))
  expect_equal( x$units, units)
  x <- x$setX(x, numeric())
  expect_equal( x$poolNames(x), c("B","R","L","RN","LN","I","RP","LP","IP","alpha"))
  expect_true( all( c("R_SOM","R_amend","RP_SOM","alpha") %in% x$stateNames(x) ))
  #.self <- x  # rm(.self)  # poolName <- "RP"
  x0 <- structure(seq_along(x$stateNames(x)), names = x$stateNames(x))
  x <- x$setX(x,x0)
  expect_equivalent( x$frac[["alpha"]], x0[c("alpha")])
  expect_equivalent( x$frac[["RP"]], x0[c("RP_SOM","RP_amend")])
  expect_equivalent( x$frac[["R"]], x0[c("R_SOM","R_amend")])
  expect_equivalent( x$tot["alpha"], x0[c("alpha")])
  expect_equivalent( x$tot["R"], sum(x0[c("R_SOM","R_amend")]*units$C))
  expect_equivalent( x$tot["RN"], sum(x0[c("RN_SOM","RN_amend")]*units$N))
  expect_equivalent( x$rel[["alpha"]], 1)
  expect_equivalent( x$rel[["R"]], x$frac[["R"]]*units$C/x$tot["R"])
  expect_equivalent( x$rel[["RN"]], x$frac[["RN"]]*units$N/x$tot["RN"])
  expect_true( all(sapply(x$rel, sum) - 1 < 1e-12) )
})

test_that("MultiPoolFractions 13C 14N, noPIso", {
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppTrillion
    , N = c(N14 = 1, N15 = 0.01) # 15N in percent
  )
  x <- createMultiPoolFractions(units, setX = createSesam4CNsetX(units))
  x <- x$setX(x, numeric())
  expect_equal( x$poolNames(x), c("B","R","L","RN","LN","I","RP","LP","IP","alpha"))
  expect_true( all( c("R_C12","R_C13","RN_N14","alpha","RP") %in% x$stateNames(x) ))
  #.self <- x  # rm(.self)  # poolName <- "RP"
  x0 <- structure(seq_along(x$stateNames(x)), names = x$stateNames(x))
  x <- x$setX(x,x0)
  expect_equivalent( x$frac[["alpha"]], x0[c("alpha")])
  expect_equivalent( x$frac[["RP"]], x0[c("RP")])
  expect_equivalent( x$frac[["R"]], x0[c("R_C12","R_C13","R_C14")])
  expect_equivalent( x$tot["alpha"], x0[c("alpha")])
  expect_equivalent( x$tot["R"], sum(x0[c("R_C12","R_C13","R_C14")]*units$C))
  expect_equivalent( x$tot["RN"], sum(x0[c("RN_N14","RN_N15")]*units$N))
  expect_equivalent( x$rel[["alpha"]], 1)
  expect_equivalent( x$rel[["R"]], x$frac[["R"]]*units$C/x$tot["R"])
  expect_equivalent( x$rel[["RN"]], x$frac[["RN"]]*units$N/x$tot["RN"])
  expect_true( all(sapply(x$rel, sum) - 1 < 1e-12) )
})




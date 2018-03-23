#require(testthat)
#test_file("tests/testthat/test_PoolFractions.R")
context("MultiPoolFractions")

test_that("MultiPoolFractions amend through all fractions", {
  elements <- c("C","N","P")
  units <- structure(lapply(elements, function(el) c(SOM = 1, amend = 1))
                     , names = elements)
  x <- createMultiPoolFractions(units, setX = createSesam4setX(units))
  expect_equal( x$units, units)
  expect_equal( names(x$tot), c("B","R","L","BN","RN","LN","I","BP","RP","LP","IP","alpha"))
  expect_true( all( c("R_SOM","R_amend","RP_SOM","alpha") %in% names(x$stateVec(x)) ))
  #.self <- x  # rm(.self)  # poolName <- "RP"
  x0 <- structure(seq_along(x$stateVec(x)), names = names(x$stateVec(x)))
  x <- x$setX(x,x0)
  expect_equivalent( x$frac[["alpha"]], x0[c("alpha")])
  expect_equivalent( x$frac[["RP"]], x0[c("RP_SOM","RP_amend")])
  expect_equivalent( x$frac[["R"]], x0[c("R_SOM","R_amend")])
  expect_equivalent( x$tot["alpha"], x0[c("alpha")])
  expect_equivalent( x$tot["R"], sum(x0[c("R_SOM","R_amend")]*units$C))
  expect_equivalent( x$tot["RN"], sum(x0[c("RN_SOM","RN_amend")]*units$N))
  expect_equivalent( x$rel[["alpha"]], 1)
  expect_equivalent( x$rel[["R"]], x$frac[["R"]]/x$tot["R"])
  expect_equivalent( x$rel[["RN"]], x$frac[["RN"]]/x$tot["RN"])
  # reli <- x$rel[["B"]]
  expect_true( all(sapply(x$rel[c("B","R","L")], function(reli) sum(reli * units$C)) - 1 < 1e-12) )
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
  expect_equal( names(x$tot), c("B","R","L","BN","RN","LN","I","BP","RP","LP","IP","alpha"))
  expect_true( all( c("R_C12","R_C13","RN_N14","alpha","RP") %in% names(x$stateVec(x)) ))
  #.self <- x  # rm(.self)  # poolName <- "RP"
  x0 <- structure(seq_along(x$stateVec(x)), names = names(x$stateVec(x)))
  x <- x$setX(x,x0)
  expect_equivalent( x$frac[["alpha"]], x0[c("alpha")])
  expect_equivalent( x$frac[["RP"]], x0[c("RP")])
  expect_equivalent( x$frac[["R"]], x0[c("R_C12","R_C13","R_C14")])
  expect_equivalent( x$tot["alpha"], x0[c("alpha")])
  expect_equivalent( x$tot["R"], sum(x0[c("R_C12","R_C13","R_C14")]*units$C))
  expect_equivalent( x$tot["RN"], sum(x0[c("RN_N14","RN_N15")]*units$N))
  expect_equivalent( x$rel[["alpha"]], 1)
  expect_equivalent( x$rel[["R"]], x$frac[["R"]]/x$tot["R"])
  expect_equivalent( x$rel[["RN"]], x$frac[["RN"]]/x$tot["RN"])
  expect_true( all(sapply(x$rel[c("B","R","L")], function(reli) sum(reli * units$C)) - 1 < 1e-12) )
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
  ds <- do.call(rbind, lapply(1:5, function(time){ data.frame(time = time, t(x0), note = "bla")}))
  #.self <- x
  dsSum <- sumMultiPoolFractions(x, ds)
  expect_equal( nrow(dsSum), nrow(ds) )
  expect_true( all(names(x$tot) %in% names(dsSum)) )
  expect_equal( dsSum$R, rowSums(ds[,c("R_C12","R_C13","R_C14")]*matrix(
    x$units$C, nrow = nrow(ds), ncol = length(x$units$C), byrow = TRUE)))
})

.tmp.f <- function(){
test_that("setMultiPoolFractionsPool", {
  units <- list(
    C = c(C12 = 1, C13 = 0.01, C14 = 1e-12) # 13C in percent, 14C ppTrillion
    , N = c(N14 = 1, N15 = 0.01) # 15N in percent
  )
  x <- createMultiPoolFractions(units, setX = createSesam4CNsetX(units))
  #.self <- x
  x0 <- structure(seq_along(x$stateVec(x)), names = names(x$stateVec(x)))
  C12 <- 100; cnB <- 8.3; cpB <- 47.3
  x1 <- setMultiPoolFractionsPool(x, x0, "B", C12
                            , c(N = cnB,  P = cpP)
                            , list(C = c(C13 = 2, C14 = 4), N = c(N15 = 8)))
  x1F <- x$setX(x, x1)
  expect_equivalent(x1["B_C12"], C12)
  expect_equivalent(x1["B_C13"], C12 + C12*2*units$C$C13)
  expect_equivalent(x1["B_C14"], C12 + C12*4*units$C$C14)
  expect_equivalent(x1F$tot["B"], sum(C12*c(1,2,4)*units$C))
  Ntot <- x1F$tot["BN"]
  expect_equivalent(Ntot*cnB, x1F$tot["B"])
  expect_equivalent(x1["B_N14"], C12)

})
}




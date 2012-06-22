#require(testthat)
context("twDev")

tmp.f <- function(){
test_that("R and C version coinside",{
		outcw2 <-  solve.modSpaceTvr(x0c,times,parmsw, useRImpl=FALSE)
		#head(outcw2)
		expect_true( max(abs(outcw2-outcw)) < 1e-8 )
		# apply( abs(outcw2-outcw), 2, max )
		# head(outcw2)
	})
}


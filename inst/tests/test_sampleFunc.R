#require(testthat)
context("sampleFunc")

test_that("two vectors are added correctly",{
		a = 1:10
		b = 1
		c = a+b
		expect_that( sampleFunc(a,b), equals(c) )
		expect_that( sampleFunc(a), throws_error() )
	})


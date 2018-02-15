## unit tests will not be done if RUnit is not available
if (require("testthat", quietly = TRUE)) {
	pkg <- "sesam" # <-- Change to package name!
	test_check(pkg)
} else {
	warning("cannot run unit tests -- package testthat is not available")
}
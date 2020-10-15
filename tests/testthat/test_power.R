context("RnaSeqSampleSize Power")
library(RnaSeqSampleSize)

test_that("est_power is working correctly", {
	expect_equal(est_power(n=63, rho=2, lambda0=5, phi0=0.5,alpha=0.05), 0.995,tolerance = 1.0e-3)
	expect_equal(est_power(n=63, rho=2, lambda0=5, phi0=0.5,f=0.01), 0.8,tolerance = 1.0e-3)
})

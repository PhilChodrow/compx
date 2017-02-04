library(compx, warn.conflicts = FALSE, quietly = TRUE)
library(tidyverse, warn.conflicts = FALSE, quietly = TRUE)

context('Divergence Functions')

test_that("Divergence equal distributions is 0", {
	expect_equal(DKL(p = rep(1/10, 10), q = rep(1/10,10)),0)
	expect_equal(DKL(p = rep(1/5, 5), q = rep(1/5,5)),0)
})

test_that("Can't compute divergence of distributions over different lengths", {
	expect_error(DKL(p = c(1/2, 1/2, 0), q = c(1/2, 1/2)))
})

p <- c(.1, .4, .5)
q <- c(.5, .2, .3)


test_that("Quadratic approximation of euclidean distance is perfect",{
	expect_equal(euc(p, q), quad(p, q, h = euc_))
})

test_that("NA_multiply implements expected logic",{

	# single dimension cases
	delta_1 <- as.matrix(c(0, NA, 1))
	H_1     <- DKL_(c(0, 1, 0))
	result  <- NA_multiply(delta_1, H_1) %>% as.numeric()
	expect(is.na(result), 'This example should generate NA')

	H_2    <- DKL_(c(.5, .5, 0))
	result <- NA_multiply(delta_1, H_2) %>% as.numeric()
	expect(is.na(result), 'This example should generate NA')

	delta_3 <- as.matrix(c(0, 0, 1))
	H_3     <- DKL_(as.numeric(delta_3))
	result  <- NA_multiply(delta_3, H_3) %>% as.numeric()
	expect(!is.na(result), 'This example should be ok')

	# Multiple dimensions

	delta_4 <- matrix(c(1, 0, 3, 3, 0, 0), 3, 2)
	H_4     <- diag(c(1, Inf, 1))
	result  <- NA_multiply(delta_4, H_4)
	expect(!is.na(result), 'This example should be ok')

	# this one should work too

	delta_5 <- matrix(c(NA, 1, 1, NA, 2, 2), 3, 2)
	H_5 <- diag(c(0, 1, 1))
	result  <- NA_multiply(delta_5, H_5)
	expected_result <- matrix(c(2, 4, 4, 8), 2, 2)
	expect_equal(result, expected_result,tolerance = 1e-10, scale = NULL, message = 'This example should be ok')

})





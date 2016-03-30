library(compx)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(mvtnorm)

context('Estimation Functions')


n <- sample(5:10, 1)
I <- sample(4:10, 1)
J <- sample(2:20, 1)
K <- sample(10:20, 1)

dims <- list(n = n, I = I, J = J, K = K)

test_that("With one representative function, all estimates are equal", {

	dims$K <- 1

	data <- random_data(dims)
	pars <- random_params(dims)

	est <- estimate(data, pars)

	compare <- do.call(cbind, replicate(I, pars$Q, simplify = FALSE))
	expect_equal(est, compare)

})

test_that("With multiple rep. functions, all results are valid",{

	dims$K <- sample(2:100, 1) # more than one

	data <- random_data(dims)
	pars <- random_params(dims)

	ests <- estimate(data, pars)
	expect_equal(min(apply(ests, MARGIN = 2, simplex_check)), 1)
})


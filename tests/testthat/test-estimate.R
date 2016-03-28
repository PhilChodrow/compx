library(compx)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(mvtnorm)

context('Estimation Functions')

test_that("With one representative function, all estimates are equal", {

	n <- sample(2:10, 1)
	K <- 1
	J <- sample(2:20, 1)

	random_data(n, K, J) %>% attach(warn.conflicts = FALSE)
	random_params(n, K, J) %>% attach(warn.conflicts = FALSE)

	est <- estimate(X = X, Q = Q, Mu = Mu, Sigma = Sigma, C = C)
	expect_equal(min(est == Q), 1)
})

test_that("With multiple rep. functions, all results are valid",{
	n     <- sample(2:10, 1)
	K     <- sample(2:100, 1) # more than one
	J     <- sample(2:100,1)

	random_data(n, K, J) %>% attach(warn.conflicts = FALSE)
	random_params(n, K, J) %>% attach(warn.conflicts = FALSE)

	ests <- estimate(X = X, Q = Q, Mu = Mu, Sigma = Sigma, C = C)
	expect_equal(min(apply(ests, MARGIN = 2, simplex_check)), 1)
})


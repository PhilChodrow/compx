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

	est <- est(data, pars)
	est
	compare <- do.call(rbind, replicate(I, pars$Q, simplify = FALSE))
	expect_equal(est, compare)
})

test_that("With multiple rep. functions, all results are valid",{

	dims$K <- sample(2:100, 1) # more than one

	data <- random_data(dims)
	pars <- random_params(dims)

	ests <- est(data, pars)

	expect_equal(dim(ests), dim(data$P))
	expect_equal(min(apply(ests, MARGIN = 1, simplex_check)), 1)
})

test_that("Spatial influence function always adds to 1",{
	data <- random_data(dims)
	pars <- random_params(dims)

	spatial_influence <- spatial_influence_constructor(pars)
	influence <- spatial_influence(data$X)
	expect_equal(rowSums(influence),rep(1, I))
})


library(compx)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(matrixcalc)

context('Random data and param generators')


n <- sample(5:10, 1)
I <- sample(4:10, 1)
J <- sample(2:20, 1)
K <- sample(10:20, 1)

dims <- list(n = n, I = I, J = J, K = K)

test_that("Random compositional data lies on simplex", {
	P <- random_distributions(dims)

	expect_equal(min(apply(P, MARGIN = 1, simplex_check)), 1)
})

test_that('Random data generator gives as many locs as distributions',{

	data <- random_data(dims)
	expect_equal(dim(data$X)[1], dim(data$P)[1])
})

test_that('Random PD matrix generator does indeed return a PD matrix',{

	M <- random_PD_matrix(n)
	eigs <- eigen(M)
	expect_equal(min(eigs$values > 0), 1)
})

test_that('Random param generator gives valid params',{

	pars <- random_params(dims)
	Q <- pars$Q
	Mu <- pars$Mu
	Sigma <- pars$Sigma
	C <- pars$C
	expect_equal(length(C), K)
	expect_equal(length(Sigma), K)
	expect_equal(all.equal(dim(Q), c(K, J)), TRUE)
	expect_equal(length(Mu), K)
	expect_equal(unlist(lapply(Mu, length)), rep(n, K))

})

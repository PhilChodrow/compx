library(compx)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(matrixcalc)

context('Check properties of random data generators')

test_that("Random compositional data lies on simplex", {
	I <- sample(1:10, 1)
	J <- sample(1:20, 1)
	P <- random_distributions(I, J)

	expect_equal(min(apply(P, MARGIN = 2, simplex_check)), 1)
})

test_that('Random data generator gives as many locs as distributions',{
	n <- sample(1:10, 1)
	I <- sample(1:10, 1)
	J <- sample(1:20, 1)

	data <- random_data(n, I, J)
	expect_equal(length(data$locs), dim(data$distributions)[2])
})

test_that('Random PD matrix generator does indeed return a PD matrix',{
	n <- sample(5:10,1)
	M <- random_PD_matrix(n)
	eigs <- eigen(M)
	expect_equal(min(eigs$values > 0), 1)
})

test_that('Random param generator gives valid aprams',{

	n <- sample(5:10,1)
	K <- sample(10:20, 1)
	J <- sample(5:10, 1)
	params <- random_params(n, K, J)

	Q = params$Q
	Mu = params$Mu
	Sigma = params$Sigma
	C = params$C

	# Shape tests
	expect_equal(length(C), K)
	expect_equal(length(Sigma), K)
	expect_equal(all.equal(dim(Q), c(J, K)), TRUE)
	expect_equal(length(Mu), K)
	expect_equal(unlist(lapply(Mu, length)), rep(n, K))

})

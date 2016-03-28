library(compx)


context("Argument conversions")

test_that("matrix ravel and unravel are inverses.",{

	n <- sample(2:10, 1)
	M <- random_PD_matrix(n)
	N <- M %>% UT_unravel() %>%
		UT_ravel()

	expect_that(M, equals(N))

	l <- n*(n+1) / 2
	v <- runif(l, 0, 1)
	u <- v %>% UT_ravel() %>%
		UT_unravel()

	expect_that(v, equals(u))
})

test_that("parameter conversion to matrix is invertible", {
	n <- sample(2:10, 1)
	I <- sample(2:40, 1)
	J <- sample(2:10, 1)
	K <- sample(10:20, 1)

	# random_params(n = n, K = K, J = J) %>% attach
	params <- random_params(n = n, K = K, J = J)
	params %>% attach

	M <- to_matrix(Q, Mu, Sigma, C)
	p <- from_matrix(M, n = n, J = J)

	expect_true(all.equal(p, params, check.attributes = FALSE))
})

test_that('conversion to vector is invertible',{
	n <- sample(2:10, 1)
	I <- sample(2:40, 1)
	J <- sample(2:10, 1)
	K <- sample(10:20, 1)

	# random_params(n = n, K = K, J = J) %>% attach
	params <- random_params(n = n, K = K, J = J)
	params %>% attach

	V <- to_vector(Q, Mu, Sigma, C)
	p <- from_vector(V, n = n, J = J, K = K)

	expect_true(all.equal(p, params, check.attributes = FALSE))
})





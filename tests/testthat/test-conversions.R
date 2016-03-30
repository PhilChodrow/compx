library(compx)

context("Argument conversions")

n <- sample(5:10, 1)
I <- sample(4:10, 1)
J <- sample(2:20, 1)
K <- sample(10:20, 1)

dims <- list(n = n, I = I, J = J, K = K)
pars <- random_params(dims)

test_that("matrix ravel and unravel are inverses.",{

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

	M <- to_matrix(pars)
	p <- from_matrix(M, dims)

	expect_true(all.equal(p, pars, check.attributes = FALSE))
})

test_that('conversion to vector is invertible',{

	V <- to_vector(pars)
	p <- from_vector(V, dims)
	v <- to_vector(p)

	expect_true(all.equal(p, pars, check.attributes = FALSE))
	expect_true(all.equal(v, V, check.attributes = FALSE))
})

test_that('matrix square root is correct',{
	A <- random_PD_matrix(10)
	B <- square_root(A)
	expect_true(all.equal(B%*%B, A))
})




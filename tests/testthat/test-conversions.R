library(compx)

context("Argument conversions")

n <- sample(2:5, 1)
I <- sample(4:10, 1)
J <- sample(2:20, 1)
K <- sample(10:20, 1)

dims <- list(n = n, I = I, J = J, K = K)

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


test_that('to and from vector are inverses',{

	b <- runif(K, 0, 1)
	V <- rnorm(n*(n+3) / 2 * K)
	Q <- rnorm(K * J) %>% matrix(K, J)
	par <- list(b = b, V = V, Q = Q)

	converted <- to_vec(par) %>% from_vec(dims)

	expect_true(all.equal(converted$b, b))
	expect_true(all.equal(converted$V, V))
	expect_true(all.equal(converted$Q, Q))
})







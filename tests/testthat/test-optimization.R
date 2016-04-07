library(compx)

context('optimization infrastructure')


I <- sample(2:10, 1)
n <- sample(2:5, 1)
K <- sample(2:5, 1)
J <- 4
dims <- list(n = n, I = I, J = J, K = K)
data <- random_data(dims)


x <- rnorm(n)
V <- rnorm(n*(n+3) / 2 * K)
b <- runif(K, 0, 1)
Q <- random_representatives(dims)
par <- list(b = b, V = V, Q = Q)
vec <- to_vec(par)

n_params <- K * (n*(n+3)/2 + 1 + J)

problem <- make_problem(data, dims)

test_that("Objective, constraints, and gradients are dimensionally correct",{
	expect_true(problem$objective(vec) > 0)
	expect_equal(dim(problem$gradient(vec)), c(1, n_params))
	expect_equal(problem$heq(vec), rep(0, K))
	expect_equal(dim(problem$grad_heq(vec)), c(K, n_params))
	expect_true(min(problem$hin(vec) > 0) == 1)
	expect_equal(dim(problem$grad_hin(vec)), c(K*(J+1), n_params))
})



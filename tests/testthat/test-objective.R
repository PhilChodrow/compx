library(compx)

context('Objective function')


test_that('Objective function returns a nonnegative scalar', {
	n <- sample(2:10, 1)
	I <- sample(2:40, 1)
	J <- sample(2:10, 1)
	K <- sample(10:20, 1)

	# generate random data and params, and place them in the namespace
	random_data(n, I = I, J = J) %>% attach(warn.conflicts = FALSE)
	random_params(n, K, J)       %>% attach(warn.conflicts = FALSE)

	objective <- make_objective(P = P, X = X)
	expect_true(objective(Q = Q, Mu = Mu, Sigma = Sigma, C = C) > 0)
})

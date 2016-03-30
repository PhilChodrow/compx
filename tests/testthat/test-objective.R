library(compx)

context('Objective function')

n <- sample(5:10, 1)
I <- sample(4:10, 1)
J <- sample(2:20, 1)
K <- sample(10:20, 1)

dims <- list(n = n, I = I, J = J, K = K)

data <- random_data(dims)
pars <- random_params(dims)

test_that("Divergence equal distributions is 0", {
	expect_equal(DKL(p = rep(1/10, 10), q = rep(1/10,10)),0)
	expect_equal(DKL(p = rep(1/5, 5), q = rep(1/5,5)),0)
})

test_that("Can't compute divergence of distributions over different lengths", {
	expect_error(DKL(p = c(1/2, 1/2, 0), q = c(1/2, 1/2)))
})

test_that("Divergence warns on nonnormalized vectors",{
	expect_warning(DKL(p = c(1,1), q = c(1,1)))
})

test_that('Total DKL between equal matrices is zero',{

	expect_equal(total_DKL(pars$Q, pars$Q), 0)
})

test_that('Objective function returns a nonnegative scalar', {

	# generate random data and params, and place them in the namespace

	objective <- make_objective(data)
	expect_true(objective(pars) > 0)
})

test_that('Problem construction is conformable',{
	v <- to_vector(pars)
	problem <- make_problem(data, dims)
	expect_true(problem$objective(v) > 0)
	expect_true(min(problem$heq(v) == 0) == 1)
	expect_true(min(problem$hin(v) > 0) == 1)
})

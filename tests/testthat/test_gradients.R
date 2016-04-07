
library(compx)
library(dplyr, warn.conflicts = FALSE)
library(numDeriv)
library(Matrix, warn.conflicts = FALSE)
library(mvtnorm)

context("gradient math")

I <- sample(2:10, 1)
n <- sample(2:5, 1)
K <- sample(2:5, 1)
J <- 4
dims <- list(n = n, I = I, J = J, K = K)

x <- rnorm(n)
V <- rnorm(n*(n+3) / 2 * K)
b <- runif(K, 0, 1)
Q <- random_representatives(dims)
par <- list(b = b, V = V, Q = Q)
vec <- to_vec(par)



test_that("normed hadamard gradient",{
	x <- rnorm(K)
	y <- rnorm(K)

	dphi <- d_phi(x,y)

	dx1 <- dphi[,1:K]
	dy1 <- dphi[,(K+1):(2*K)]

	dx2 <- jacobian(function(x) phi(x,y), x)
	dy2 <- jacobian(function(y) phi(x,y), y)

	expect_true(all.equal(dx1, dx2))
	expect_true(all.equal(dy1, dy2))
})


test_that("spatial gradient is correct",{
	d1 <- d_psi(x, vec, dims)
	d2 <- jacobian(function(vec) psi(x, vec, dims), vec)
	expect_true(all.equal(d1, d2))
})


test_that("analytic gradient equals numeric gradient of objective function",{

	data <- random_data(dims)
	obj <- obj_constructor(data = data, dims = dims)
	gradient <- grad_constructor(data, dims)
	par <- list(b = b, V = V, Q = Q)
	vec <- to_vec(par)

	grad1 <- gradient(vec)
	grad2 <- jacobian(obj, vec)

	expect_true(all.equal(grad1, grad2, check.attributes = FALSE))

})













library(compx)
library(dplyr)
library(numDeriv)

context("gradient math")

test_that("chain rule logic works for normal distribution",{
	# reference: mask for dmv norm that takes pars with square root of covariance matrix
	dmvnorm2 <- function(par, n, xx){
		mu <- par[1:n]
		sig <- par[(n+1):length(par)] %>% UT_ravel()
		Sigma <- t(sig) %*% sig
		dmvnorm(x = xx, mean = mu, sigma = Sigma)
	}

	par <- rnorm(5)
	n <- 2
	xx <- rnorm(2)

	# test derivatives
	expect_true(all.equal(d_lambda(par, n, xx),
						  jacobian(func=dmvnorm2, x=par, n = n, xx = xx),
						  check.attributes = FALSE))
})


















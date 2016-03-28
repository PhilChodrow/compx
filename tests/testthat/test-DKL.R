library(compx)

context('KL Divergence')

test_that("Divergence equal distributions is 0", {
	expect_equal(DKL(p = rep(1/10, 10), q = rep(1/10,10)),0)
	expect_equal(DKL(p = rep(1/5, 5), q = rep(1/5,5)),0)
})

test_that("Can't compute divergence of distributions over different lengths", {
	expect_error(DKL(p = c(1/2, 1/2, 0), q = c(1/2, 1/2)))
})

test_that("Error on nonnormalized vectors",{
	expect_error(DKL(p = c(1,1), q = c(1,1)))
})

test_that('Total DKL between equal matrices is zero',{
	K     <- 4 # number of representative distributions
	J     <- sample(1:100,1)
	Q     <- rep(1,K) %>%
		lapply(FUN = function(x) runif(J, 1, 100) %>% simplex_normalize()) %>%
		do.call(cbind, .)

	expect_equal(total_DKL(Q, Q), 0)
})






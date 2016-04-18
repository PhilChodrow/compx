library(compx)

context('Objective function')

n <- sample(5:10, 1)
I <- sample(4:10, 1)
J <- sample(2:20, 1)
K <- sample(10:20, 1)

dims <- list(n = n, I = I, J = J, K = K)


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

test_that("Divergence works if p has zeros",{
	p <- runif(10, 0,1)
	p[1] <- 0
	p <- p / sum(p)

	q <- runif(10, 0, 1)
	q <- q/sum(q)

	expect_true(!is.na(DKL(p,q)))
	expect_true(is.infinite(DKL(q,p)))
})

test_that("Divergence across two matrices is zero only if they are the same",{
	P <- runif(K * J, 0,1) %>%
		matrix(K,J) %>%
		apply(MARGIN = 1, FUN = simplex_normalize) %>%
		t

	Q <- runif(K * J, 0,1) %>%
		matrix(K,J) %>%
		apply(MARGIN = 1, FUN = simplex_normalize) %>%
		t

	obj1 <- 1:dim(P)[1] %>%
		matrix %>%
		apply(MARGIN = 1, FUN = function(i) DKL(P[i,], Q[i,])) %>%
		sum

	obj2 <- 1:dim(P)[1] %>%
		matrix %>%
		apply(MARGIN = 1, FUN = function(i) DKL(P[i,], P[i,])) %>%
		sum

	expect_true(obj1 > 0)
	expect_true(obj2 == 0)
	# Ok, this appears to work correctly, so what we should do now is use this same idiom to calculate the objective function (i.e. first calculate the estimate matrix and then iterate over it.)

})


test_that("Entropy is divergence from uniform",{
	p <- runif(10, 0,1)
	p <- p / sum(p)
	u <- rep(1/length(p), length(p))
	expect_equal(H(p), log(length(p)) - DKL(p,u))
})

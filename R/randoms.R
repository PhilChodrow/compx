#' randoms
#' @name randoms
#' @docType package
#' @import dplyr

NULL

#' Create a random set of params that are valid for the objective function
#' @param n the spatial dimension
#' @param K the number of representative distributions
#' @param J the number of categories
#' @export

random_params <- function(dims){
	J <- dims$J
	K <- dims$K
	n <- dims$n

	Q     <- random_representatives(dims) #
	Mu    <- lapply(rep(1,K), function(x) runif(n, -5, 5))
	Sigma <- lapply(rep(1,K), function(x) random_PD_matrix(n))
	C     <- c(abs(rnorm(K))) # weight vector
	return(list(Q = Q, Mu = Mu, Sigma = Sigma, C = C))
}

#' Create a random PD matrix
#' https://stat.ethz.ch/pipermail/r-help/2008-February/153708.html
#' @param n the dimension
#' @param ev specified eigenvalues
#' @return a random positive-definite matrix
#' @export

random_PD_matrix <- function (n, ev = runif(n, 1, 10))
{
	Z <- matrix(ncol=n, rnorm(n^2))
	decomp <- qr(Z)
	Q <- qr.Q(decomp)
	R <- qr.R(decomp)
	d <- diag(R)
	ph <- d / abs(d)
	O <- Q %*% diag(ph)
	Z <- t(O) %*% diag(ev) %*% O
	return(Z)
}



#' Make a matrix whose columns are random valid probability distributions
#' @param I the number of distributions (columns)
#' @param J the number of bins (rows)
#' @export

random_representatives <- function(dims){
	K <- dims$K
	J <- dims$J
	rep(1, K) %>%
		lapply(FUN = function(x) runif(J, 1, 100) %>% simplex_normalize()) %>%
		do.call(rbind, .)

}



#' Make a matrix whose columns are random valid probability distributions
#' @param I the number of distributions (columns)
#' @param J the number of bins (rows)
#' @export

random_distributions <- function(dims){
	I <- dims$I
	J <- dims$J

	rep(1, I) %>%
		lapply(FUN = function(x) runif(J, 1, 100) %>% simplex_normalize()) %>%
		do.call(rbind, .)

}

#' Generate a list of random locations
#' @param n the dimension of space in which we are working
#' @param I the number of locations to generate
#' @export

random_locs <- function(dims){
	n <- dims$n
	I <- dims$I
	X <- matrix(runif(n*I, -1, 1), I, n)
	X
}

#' Make a random, evenly-spaced, spatially structured data set in one dimension.
#' @export

random_1d_data <- function(dims){

	q <- runif(dims$K*dims$J,0,1) %>%
		matrix(dims$K,dims$J) %>%
		apply(MARGIN = 1, FUN = simplex_normalize) %>%
		t
	mu <- 1:dims$K %>% matrix

	sig <- .5
	b <- runif(dims$K, 0, 1)
	v <- cbind(mu, sig) %>% t %>% c
	x <- seq(0,dims$K+1,by = .1) %>% matrix

	pars <- list(Q = q, b = b, V = v)
	gen_vec <- to_vec(pars)
	p <- x %>% apply(MARGIN = 1,
					 FUN = psi,
					 vec = gen_vec,
					 dims = dims) %>% t

	list(P = p, X = x)
}


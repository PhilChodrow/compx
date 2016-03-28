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

random_params <- function(n, K, J){
	Q     <- random_distributions(K, J) #
	Mu    <- lapply(rep(1,K), function(x) runif(n))
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

random_distributions <- function(I, J){
	rep(1, I) %>%
		lapply(FUN = function(x) runif(J, 1, 100) %>% simplex_normalize()) %>%
		do.call(cbind, .)

}

#' Generate a list of random locations
#' @param n the dimension of space in which we are working
#' @param I the number of locations to generate
#' @export

random_locs <- function(n, I){
	X <- lapply(rep(1,I), function(x) runif(n))
	X
}

#' Generate random data for analysis
#' @param n the spatial dimension
#' @param I the number of locations
#' @param J the number of bins
#' @return data a list of [1] locations and [2] corresponding distributions
#' @export

random_data <- function(n, I, J){
	list(X = random_locs(n, I), P = random_distributions(I,J))
}

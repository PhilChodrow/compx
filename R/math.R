#' math
#' @name math
#' @docType package
#' @import dplyr mvtnorm

NULL

#' Find the Kullback-Leibler divergence of two empirical distributions.
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @return numeric the KL divergence between p and q
#' @export

DKL <- function(p,q){
	return(p %*% log(p/q))
}

#' Find the gradient of the part of the Kullback-Leibler divergence that depends on f
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @return p dot log(q)
#' @export

grad_f <- function(p,q){
	return( - p %*% log(q))
}

#' Generate some random data with correct structure for testing functions
#' @param the number of locations x
#' @param the spatial dimension
#' @param the number of characteristic distributions to use in modeling
#' @return list a list of locations x, means Mu, covariance matrices Sigma, and representatives Q
#' @export
create_test_data <- function(n, d, K){
	print("Not implemented")
}

#' Compute a vector of values of different normal densities at a point x
#' @param x vector the location
#' @param Mu a list of means
#' @param Sigma a list of covariance matrices
#' @return a vector of densities at x corresponding to the different densities
#' @export
normal_vec <- function(x, Mu, Sigma){
	mapply(mvtnorm::dmvnorm, x = x, mean = Mu, sigma = Sigma)
}

#' Find the spatially-structured estimate at a location x
#' @param x vector the location
#' @param Q matrix the matrix of representative distributions
#' @param Mu list a list of means
#' @param Sigma list a list of covariance matrices
#' @param C vector the vector of scale coefficients
#' @return vector the estimates at x
#' @export

estimate <- function(x, Q, Mu, Sigma, C = 1){
	densities <- normal_vec(x, Mu, Sigma)
	vec <- C * densities / C %*% densities
	return(Q %*% vec)
}

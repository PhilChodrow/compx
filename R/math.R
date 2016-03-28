#' math
#' @name math
#' @docType package
#' @import dplyr mvtnorm matrixcalc

NULL

#' Find the Kullback-Leibler divergence of two empirical distributions.
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @return numeric the KL divergence between p and q
#' @export

DKL <- function(p,q){
	if(!simplex_check(p) | !simplex_check(q)){
		stop('p or q are not on the simplex; try simplex_normalize(p)')
	}
	if(length(p) != length(q)){
		stop('Distribution alphabets are different size')
		}
	return(as.numeric(p %*% log(p/q)))
}

#' Check whether a vector is an element of the probability simplex
#' @param p vector the vector to check
#' @param allow_zero boolean whether to allow elements with zero entries (boundary of simplex)
#' @return boolean whether p is a valid probability distribution
#' @export

simplex_check <- function(p, allow_zero = TRUE){
	nonneg <- ifelse(allow_zero, min(p >= 0), min(p > 0))
	normed <- abs(sum(p) - 1) < .00000001 # numerical tolerance
	nonneg & normed
}


#' Find the gradient of the part of the Kullback-Leibler divergence that depends on f
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @return p dot log(q)
#' @export

grad_f <- function(p,q){
	return( - p %*% log(q))
}


#' Compute a vector of values of different normal densities at a point x
#' @param x vector the location
#' @param Mu a list of means
#' @param Sigma a list of covariance matrices
#' @return a vector of densities at x corresponding to the different densities
#' @export

normal_vec <- function(x, Mu, Sigma){
#
# 	mapply(mvtnorm::dmvnorm, x = x, mean = Mu, sigma = Sigma)
#
	f <- function(Mu, Sigma){
		mvtnorm::dmvnorm(x, Mu, Sigma)
	}
	mapply(f, Mu, Sigma)

}

#' Find the spatially-structured estimate at a location x
#' @param x vector the location
#' @param Q matrix the matrix of representative distributions
#' @param Mu list a list of means
#' @param Sigma list a list of covariance matrices
#' @param C vector the vector of scale coefficients
#' @return vector the estimates at x

single_estimate <- function(x, Q, Mu, Sigma, C = 1){
	densities <- normal_vec(x, Mu, Sigma)
	vec <- C * densities / C %*% densities
	Q %*% vec
}

#' Find the spatially-structured estimate at multiple locations X
#' @param X list of location vectors
#' @param Q matrix the matrix of representative distributions
#' @param Mu list a list of means
#' @param Sigma list a list of covariance matrices
#' @param C vector the vector of scale coefficients
#' @return vector the estimates at x
#' @export

estimate <- function(X, Q, Mu, Sigma, C = 1){
	ests <- lapply(X, FUN = single_estimate, Q = Q, Mu = Mu, Sigma = Sigma, C = C) %>%
		do.call(cbind, .)
	return(ests)
}

#' Normalize a nonnegative vector so that it lies in the probability simplex
#' @param p vector the vector
#' @export

simplex_normalize <- function(p){
	if (min(p >= 0) == 0){
		stop("Vector has negative entries")
	}
	normed <- p / sum(p)
	if (min(p > 0) == 0){
		print('Distribution lies on simplex boundary')
	}
	normed
}

#' Compute the total KL divergence between a set of estimates and the data
#' @param matrix P a matrix in which each column is a distribution
#' @param matrix Q a matrix of estimates
#' @return numeric div the total divergence of the columns of Q from the columns of P
#' @export

total_DKL <- function(P, Q){
	1:dim(P)[2] %>%
		matrix %>%
		apply(1, function(i) DKL(P[,i], Q[,i])) %>%
		sum
}

#' Create the objective function for subsequent optimization
#' @param matrix P the matrix of true distributions
#' @export

make_objective <- function(P, X){
	objective <- function(Q, Mu, Sigma, C){
		ests <- estimate(X = X, Q = Q, Mu = Mu, Sigma = Sigma, C = C)
		total_DKL(P, ests)
	}
	objective
}

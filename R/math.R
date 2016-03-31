#' math
#' @name math
#' @docType package
#' @import dplyr mvtnorm matrixcalc

NULL


#' Check whether a vector is an element of the probability simplex
#' @param p vector the vector to check
#' @param allow_zero boolean whether to allow elements with zero entries (boundary of simplex)
#' @return boolean whether p is a valid probability distribution
#' @export

simplex_check <- function(p, allow_zero = TRUE){
	nonneg <- ifelse(allow_zero, min(p >= 0), min(p > 0))
	normed <- abs(sum(p) - 1) < .1 # numerical tolerance
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

normal_vec <- function(x, pars){

	f <- function(Mu, Sigma){
		mvtnorm::dmvnorm(x, Mu, Sigma)
	}
	mapply(f, pars$Mu, pars$Sigma)
}

#' Spatial influence function
#' @param x
#' @param pars
#' @export

spatial_influence_constructor <- function(pars){
	single <- function(x){
		densities <- normal_vec(x, pars)
		pars$C * densities / as.numeric(pars$C %*% densities)
	}
	influence <- function(X){
		out <- apply(X, MARGIN = 1, FUN = single)
		if(dim(pars$Q)[1] == 1){
			return(out %>% as.matrix)
		}
		return(out %>% t)
	}
	return(influence)
}


#' Find the spatially-structured estimate at multiple locations X
#' @param X matrix of location vectors
#' @param Q matrix the matrix of representative distributions
#' @param Mu list a list of means
#' @param Sigma list a list of covariance matrices
#' @param C vector the vector of scale coefficients
#' @return vector the estimates at x
#' @export

est <- function(data, pars){
	spatial_influence <- spatial_influence_constructor(pars)
	spatial <- spatial_influence(data$X)
	return(spatial %*% pars$Q)
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
		warning('p lies on simplex boundary')
	}
	normed
}



#' Find the Kullback-Leibler divergence of two empirical distributions.
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @return numeric the KL divergence between p and q
#' @export

DKL <- function(p,q){

	if(!simplex_check(p) | !simplex_check(q)){
		message <- paste0('p or q are not on the simplex: sum(p) = ', sum(p), ' and sum(q) = ', sum(q))
		warning(message)
	}
	if(length(p) != length(q)){
		stop('Distribution alphabets are different size')
	}

	drop <- p < 10^(-10)
	p <- p[!drop]
	q <- q[!drop]

	return(as.numeric(p %*% log(p/q)))
}

#' Compute the total objective value between two matrices
#' @param matrix P a matrix in which each column is a distribution
#' @param matrix Q a matrix of estimates
#' @return numeric div the total divergence of the columns of Q from the columns of P
#' @export

total_obj_constructor <- function(fun = DKL){
	total_obj <- function(P,Q){
		1:dim(P)[1] %>%
			matrix %>%
			apply(1, function(i) fun(P[i,], Q[i,])) %>%
			sum
	}
}

#' Compute the sum of squares distance between two probability distributions
#' @param p
#' @param q
#' return numeric the distance
#' @export

sum_of_squares <- function(p,q){
	if(length(p) != length(q)){
		stop('Distribution alphabets are different size')
	}
	(p - q)^2 %>% sum
}

#' Compute the entropy of a vector
#' @param the vector to compute the entropy of
#' @export
H <- function(p){
	- p %*% log(p) %>% as.numeric
}





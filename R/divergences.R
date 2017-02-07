#' divergences
#' @name divergences
#' @docType package
#' @import tidyverse
NULL

#' Find the Kullback-Leibler divergence of two empirical distributions.
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @return numeric the divergence between p and q
#' @export

DKL <- function(n,m){
	if(length(n) != length(m)){
		stop('Distribution alphabets are different size')
	}

	p <- simplex_normalize(n)
	q <- simplex_normalize(m)

	return(sum(p * log(p/q), na.rm = T))
}

#' @inherit DKL
#' @export
DJS <- function(n, m){

	p <- simplex_normalize(n)
	q <- simplex_normalize(m)
	r <- simplex_normalize(m + n)

	w <- sum(n) / sum(n + m)

	w*DKL(p, r) + (1-w) * DKL(q, r)
}

#' @inherit DKL
#' @export
euc <- function(n,m){

	p <- simplex_normalize(n)
	q <- simplex_normalize(m)

	1/2 * sum((p-q)^2)
}

# Agglomeration Loss ----------------------------------------------------------

#' Compute the agglomeration loss associated with combining two sets of counts
#' under a user-specified divergence
#' @param n the first set of counts to agglomerate
#' @param m the second set of counts to agglomerate
#' @param marginal the global count set
#' @param agg_div the divergence function under which to compute loss
#' @param ... additional arguments passed to div
#' @return a nonnegative scalar
#' @export
agg_loss <- function(n, m, marginal = n + m, div = DKL){

	r <- n + m

	1/sum(marginal) *
		(sum(n)     * div(n, marginal) +
		 sum(m)     * div(m, marginal) -
		 sum(n + m) * div(r, marginal))
}


# Hessians --------------------------------------------------------------------
# The functions in this section generate hessian matrices corresponding to
# various Bregman divergences.

#' Compute the Hessian of a Bregman divergence generator at a point
#' @param base the base point at which to calculate the Hessian
#' @export
DKL_ <- function(base){
	p <- simplex_normalize(base)
	diag(1/p)
}

#' Euclidean
#' @inherit DKL_
euc_ <- function(base){

	diag(length(base))
}

#' @inherit DKL_
#' @export
cum_euc_ <- function(base){
	lower <- matrix(0, length(base), length(base))
	lower[lower.tri(lower, diag = T)] <- 1
	lower %*% t(lower)
}


# Quadratic Approximators -----------------------------------------------------

#' Quadratic form multiplication with NA and zero logic. The modified logic is:
#' 0 * Inf = 0
#' 0 * NA = 0
#' @param delta a vector or matrix
#' @param H a square matrix
#' @return a matrix, the product t(delta) %*% .H %*% delta
#' @export

NA_multiply <- function(delta, H){
	mat                       <- delta %*% t(delta)
	zero_ind                  <- mat == 0
	zero_ind[is.na(zero_ind)] <- FALSE
	na_ind                    <- is.na(mat)

	.H           <- H
	.H[zero_ind] <- 0

	if(sum(H[na_ind] == 0) == sum(na_ind)){
		delta[is.na(delta)] <- 0
	}

	t(delta) %*% .H %*% delta
}

#' Quadratic approximation to a Bregman divergence
#' @param n the first point
#' @param m the second point
#' @param base the point at which to compute the Hessian
#' @param h the hessian function. Will be passed base as argument
#' @return the quadratic approximation of the Bregman divergence of which h is
#' the hessian.
#' @export

quad <- function(n, m, base = m, h = euc_){

	p <- simplex_normalize(n)
	q <- simplex_normalize(m)

	delta <- as.matrix(p-q)
	(1/2 * NA_multiply(delta, h(base))) %>% as.numeric()
}

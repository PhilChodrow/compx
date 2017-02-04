#' divergences
#' @name divergences
#' @docType package
#' @import tidyverse
NULL

#' Find the Kullback-Leibler divergence of two empirical distributions.
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @param check whether to check that p and q are valid probability distributions
#' @return numeric the KL divergence between p and q
#' @export

DKL <- function(p,q){
	if(length(p) != length(q)){
		stop('Distribution alphabets are different size')
	}
	return(sum(p * log(p/q), na.rm = T))
}

#' @export
DJS <- function(p, q, w = NULL){
	if(is.null(w)){
		w <- sum(p) / (sum(p) + sum(q))
	}
	r <- w * p + (1-w) * q
	w*DKL(p, r) + (1-w) * DKL(q, r)
}

euc <- function(p,q){
	1/2 * sum((p-q)^2)
}

# Agglomeration Loss ----------------------------------------------------------

agg_loss <- function(p, q, w = 1/2, marginal = w*p + (1-w)*q, div = DKL, ...){
	r <- w*p + (1-w) * q
	w * div(p,r, ...) + (1-w)*div(q,r, ...) - div(r, marginal, ...)
}


# Hessians --------------------------------------------------------------------

DKL_ <- function(base){
	diag(1/base)
}

euc_ <- function(base){
	diag(length(base))
}

cum_euc_ <- function(base){
	lower <- matrix(0, length(base), length(base))
	lower[lower.tri(lower, diag = T)] <- 1
	lower %*% t(lower)
}


# Quadratic Approximators -----------------------------------------------------

#' should implement NA*0 = 0 and 0 * inf = 0 logic.
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

quad <- function(p, q, base = (p+q)/2, h = euc_){
	delta <- as.matrix(p-q)
	(1/2 * NA_multiply(delta, h(base))) %>% as.numeric()
}

# tomorrow morning:
# add this stuff into the compx source code to enable flexible, extensible analysis with a wide variety of functions. See if we can implement some interesting calculations for clustering this way, such as pairwise agglomeration losses as
# our distance measure.



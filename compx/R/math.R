#' math
#' @name math
#' @docType package
#' @import dplyr

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

#'

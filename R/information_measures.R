#' info
#' @name info
#' @docType package
#' @import tidyverse

NULL

# -----------------------------------------------------------------------------
# INFORMATION MEASURES
# -----------------------------------------------------------------------------

#' Check whether a vector is an element of the probability simplex
#' @param p vector the vector to check
#' @param allow_zero boolean whether to allow elements with zero entries (boundary of simplex)
#' @return boolean whether p is a valid probability distribution
#' @export

simplex_check <- function(p, allow_zero = TRUE){
	if(min(is.na(p)) == 1) return(FALSE)
	nonneg <- ifelse(allow_zero, min(p >= 0), min(p > 0))
	normed <- abs(sum(p) - 1) < .01 # numerical tolerance
	nonneg & normed
}

#' Normalize a nonnegative vector so that it lies in the probability simplex
#' @param p vector the vector
#' @export

simplex_normalize <- function(p){
	if (min(p >= 0) == 0){
		stop("Vector has negative entries")
	}

	normed <- p / sum(p)
	# if (min(p > 0) == 0){
	# 	warning('p lies on simplex boundary')
	# }
	normed
}

#' Find the Kullback-Leibler divergence of two empirical distributions.
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @param check whether to check that p and q are valid probability distributions
#' @return numeric the KL divergence between p and q
#' @export
DKL <- function(p,q, check = FALSE){
	# if(!is.numeric(p) | !is.numeric(q)){
	# 	return(NaN)
	# }
	if(check){
		if(!simplex_check(p) | !simplex_check(q)){
			message <- paste0('p or q are not on the simplex: sum(p) = ', sum(p), ' and sum(q) = ', sum(q))
			warning(message)
		}
	}
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



#' Find the entropy of a distribution
#' @param p a normalized vector of nonnegative probabilities.
#' @export
H <- function(p){
	if(!simplex_check(p)){
		message <- paste0('p is not on the simplex: sum(p) = ', sum(p))
		warning(message)
	}
	drop <- p < 10^(-10)
	p <- p[!drop]
	as.numeric(- p %*% log(p))
}
#' Compute the binary entropy function for a fixed proportion.
#' @param p the proportion
#' @export
H_B <- function(p){
	-(p * log(p) + (1-p) * log(1-p))
}

#' Compute the mutual information of a data frame or matrix of crosstabs
#' @param input a data frame or matrix representing a crosstab
#' @param drop_threshold numerical keyword to KL divergence avoiding 0*0 multiplication. Should not be modified.
#' @export
mutual_info <- function(input){
	if(class(input) != 'matrix'){
		input <- as.matrix(input)
	}
	ind <- as.matrix(rowSums(input) / sum(input)) %*%
		   t(as.matrix(colSums(input)/sum(input)))              # margin product
	DKL(input/sum(input), ind)
}

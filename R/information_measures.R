#' info
#' @name info
#' @docType package
#' @import dplyr tidyr

NULL

# -----------------------------------------------------------------------------
# INFORMATION GEOMETRY
# -----------------------------------------------------------------------------

#' Check whether a vector is an element of the probability simplex
#' @param p vector the vector to check
#' @param allow_zero boolean whether to allow elements with zero entries (boundary of simplex)
#' @return boolean whether p is a valid probability distribution
#' @export

simplex_check <- function(p, allow_zero = TRUE){
	if(min(is.na(p)) == 1) return(FALSE)
	nonneg <- ifelse(allow_zero, min(p >= 0), min(p > 0))
	normed <- abs(sum(p) - 1) < .1 # numerical tolerance
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

	drop <- p < 10^(-10)
	p <- p[!drop]
	q <- q[!drop]

	return(as.numeric(p %*% log(p/q)))
}

euc <- function(p,q){
	sum((p-q)^2)
}


#' Find the entropy of a distribution
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

#' @export

H_B <- function(proportion){
	-(proportion * log(proportion) + (1-proportion) * log(1-proportion))
}


# -----------------------------------------------------------------------------
# MUTUAL INFORMATION
# -----------------------------------------------------------------------------

#' Compute joint distribution for a grouped table.
#' @export

compute_distributions <- function(data, group_col = NULL){

	if(is.null(group_col)){
		data$group_col <- row.names(data)
		group_col = 'group_col'
	}

	category_cols <- names(data)[names(data) != group_col]

	gathered <- data %>% group_by_(group_col) %>%
		summarise_each(funs(sum)) %>%
		gather_( key_col = 'category', value_col = 'n', gather_cols = category_cols) %>%
		mutate(p = n / sum(n)) %>%
		select(-n)

	joint <- gathered %>%
		spread_(key_col = 'category', value_col = 'p', fill = 0) %>%
		select_(.dots = category_cols) %>%
		as.matrix()

	group <- gathered %>%
		group_by_(group_col) %>%
		summarise(p = sum(p)) %>%
		select(p) %>%
		as.matrix()

	category <- gathered %>%
		group_by_('category') %>%
		summarise(p = sum(p)) %>%
		spread(category, p) %>%
		select_(.dots = category_cols) %>%
		gather(key = category, value = p) %>%
		select(p) %>%
		as.matrix()

	list(joint = joint, group = group, category = category)
}

#' Compute mutual information for grouped data.
#' @export
mutual_info <- function(data, group_col = NULL){

	if(nrow(data) < 2) return(0) # would NA be better?

	distributions <- compute_distributions(data, group_col)

	independent <- distributions$group %*% t(distributions$category)
	(distributions$joint * log(distributions$joint / independent)) %>%
		sum(na.rm = T)

	# f <- function(i){
	# 	DKL(distributions$joint[i,] / sum(distributions$joint[i,]),
	# 		distributions$category)
	# }
	#
	# divs <- 1:length(distributions$group) %>%
	# 	as.matrix() %>%
	# 	apply(MARGIN = 1,
	# 		  FUN = f)
	# divs %*% distributions$group %>% c
}


#' @export
entropy <- function(data){
	if(nrow(data) < 2) return(0) # would NA be better?
	distributions <- compute_distributions(data)
	H(distributions$category)
}





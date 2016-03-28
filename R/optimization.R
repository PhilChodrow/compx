#' optimization
#' @name optimization
#' @docType package
#' @import dplyr matrixcalc

NULL

#' Unravel the upper triangular part of a matrix as a vector
#' @param M matrix
#' @return v a vector
#' @export

UT_unravel <- function(M){
	if(!isSymmetric(M)){
		warning("Warning: matrix M is not symmetric, information lost")
	}
	M[lower.tri(t(M), diag = TRUE)]
}

#' Find the dimensions of a matrix for a given size of the upper-triangular part
#' @param l, the triangular number to convert
#' @return m, the dimension of the corresponding matrix.
#' @export

find_nrow <- function(v){
	n <- (-1 + sqrt(1 + 8*v))/2
	if(n %% 1 != 0){
		stop("v is not a triangular number")
	}
	n
}

#' Inverse of UT_unravel
#' @export

UT_ravel <- function(v){
	l <- length(v)
	n <- find_nrow(l)
	x <- matrix(0, n, n)
	x[lower.tri(x, diag = TRUE)] <- v
	x <- t(x) + x - diag(diag(x))
	x
}

#' Convert a list of params to matrix form
#' @param Q matrix the matrix of representative distributions
#' @param Mu list a list of means
#' @param Sigma list a list of covariance matrices
#' @param C vector the vector of scale coefficients
#' @return M a matrix
#' @export
to_matrix <- function(Q, Mu, Sigma, C){
	mu <- do.call(cbind, Mu)
	sigma <- do.call(cbind, lapply(Sigma, UT_unravel))
	rbind(Q, mu, sigma, C)
}

#' Convert a matrix to a list of params
#' @param M the matrix to convert
#' @param n the spatial dimension
#' @param J the number of buckets
#' @return param a list of parameters
#' @export
from_matrix <- function(M, n, J){
	Q <- M[1:J,]
	Mu <- M[(J+1):(J+n),] %>% split(col(.))
	Sigma <- M[(J + n + 1):(J + n + n*(n+1)/2),] %>%
		plyr::alply(.margins = 2, .fun = UT_ravel)
	C <- M[nrow(M),]

	list(Q = Q, Mu = Mu, Sigma = Sigma, C = C)
}

#' Convert a list of params to vector form
#' @param Q matrix the matrix of representative distributions
#' @param Mu list a list of means
#' @param Sigma list a list of covariance matrices
#' @param C vector the vector of scale coefficients
#' @return v a vector
#' @export

to_vector <- function(Q, Mu, Sigma, C){
	to_matrix(Q = Q, Mu = Mu, Sigma = Sigma, C = C) %>%
		as.vector()
}

#' Convert a vector to a list of params
#' @param M the vector to convert
#' @param n the spatial dimension
#' @param J the number of buckets
#' @return param a list of parameters
#' @export
from_vector <- function(V, n, J, K){
	n_rows <- J + n + n*(n+1)/2 + 1
	matrix(V, n_rows, K) %>%
		from_matrix(n = n, J = J)
}


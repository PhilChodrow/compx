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
	if(!all.equal(M, t(M))){
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
to_matrix <- function(pars){
	mu <- do.call(rbind, pars$Mu)
	sigma <- lapply(pars$Sigma, square_root)
	sigma <- do.call(rbind, lapply(sigma, UT_unravel))
	cbind(pars$Q, mu, sigma, pars$C)
}

#' Convert a matrix to a list of params
#' @param M the matrix to convert
#' @param n the spatial dimension
#' @param J the number of buckets
#' @return param a list of parameters
#' @export
from_matrix <- function(M, dims){
	J <- dims$J
	n <- dims$n

	Q <- M[,1:J]
	Mu <- M[,(J+1):(J+n)] %>% split(row(.))
	Sigma <- M[,(J + n + 1):(J + n + n*(n+1)/2)] %>%
		plyr::alply(.margins = 1, .fun = UT_ravel)
	Sigma <- lapply(Sigma, function(S) S %*% S)
	C <- M[,ncol(M)]
	list(Q = Q, Mu = Mu, Sigma = Sigma, C = C)
}

#' Convert a list of params to vector form
#' @param Q matrix the matrix of representative distributions
#' @param Mu list a list of means
#' @param Sigma list a list of covariance matrices
#' @param C vector the vector of scale coefficients
#' @return v a vector
#' @export

to_vector <- function(pars){
	to_matrix(pars) %>%
		as.vector()
}

#' Convert a vector to a list of params
#' @param M the vector to convert
#' @param n the spatial dimension
#' @param J the number of buckets
#' @return param a list of parameters
#' @export

from_vector <- function(V, dims){
	n <- dims$n
	J <- dims$J
	K <- dims$K
	n_cols <- J + n + n*(n+1)/2 + 1
	matrix(V, K, n_cols) %>%
		from_matrix(dims)
}

#' Compute a symmetric square root of a positive-definite, symmetric matrix
#' @param A the matrix to factor
#' @return B a symmetric square root of A
#' @export
square_root <- function(A){
	e <- eigen(A)
	V <- e$vectors
	V %*% diag(sqrt(e$values)) %*% t(V)
}


#' Create the objective function for subsequent optimization
#' @param matrix P the matrix of true distributions
#' @export

make_objective <- function(data, obj_fun){
	objective <- function(pars){
		total_obj <- total_obj_constructor(fun = obj_fun)
		ests <- est(data, pars)
		total_obj(data$P, ests)
	}
	objective
}


#' Consider doing this with a single overall function that returns a list of functions that will then serve as the objective and constraints.

#' Create the objective and constraints for the problem
#' @export
make_problem <- function(data, dims, obj_fun){
	m_objective <- make_objective(data, obj_fun) # matrix objective

	objective <- function(V){
		pars <- from_vector(V, dims)
		m_objective(pars)
	}

	heq <- function(V){
		pars <- from_vector(V = V, dims)
		rowSums(pars$Q) - 1
	}

	hin <- function(V){
		pars <- from_vector(V = V, dims)
		cbind(pars$Q, pars$C) %>% as.numeric()
	}
	list(objective = objective, heq = heq, hin = hin)
}

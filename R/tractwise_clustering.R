#' tractwise_clustering
#' @name tractwise_clustering
#' @docType package
#' @import tidyverse sp maptools rgeos
NULL

#' Construct an affinity matrix for a graph g based on the dist edge attribute
#' @param g an igraph::graph with a dist edge attribute/
#' @param sigma the multiplicative constant in the quadratic exponential
#' affinity kernel
#' @return a matrix containing the affinities computed with the quadratic exponential
#' affinity kernel
#' @export
affinity_matrix <- function(g, sigma = 1){
	dists <- distances(g, weights = E(g)$dist^2)
	A <- exp(-dists * sigma)
	A[is.na(A)] <- 0
	A
}

#' Construct a generalized laplacian matrix
#' @param input either an affinity matrix or an igraph::graph. If
#' class(input) == 'igraph', then an affinity matrix is computed and used.
#' Otherwise, it is assumed that the input is an affinity matrix.
#' @param ... additional arguments passed to affinity_matrix if
#' class(input) == 'igraph'.
#' @return a generalized laplacian matrix suitable for spectral analysis.
#' @details This function involves a matrix inversion and is therefore fairly
#' expensive for large graphs.
#' @export
generalized_laplacian_matrix <- function(input, ...){

	if(class(input) == 'igraph'){
		A <- affinity_matrix(input, ...)
	}else{
		A <- input
	}

	D <- diag(rowSums(A))
	L <- solve(D) %*% (D - A)
}






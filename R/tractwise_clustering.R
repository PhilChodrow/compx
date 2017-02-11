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
	assertthat::assert_that(class(input) %in% c('matrix', 'igraph'))
	if(class(input) == 'igraph'){
		A <- affinity_matrix(input, ...)
	}else{
		A <- input
	}

	D <- diag(rowSums(A))
	L <- solve(D) %*% (D - A)
}

#' Cluster a graph or affinity matrix.
#' Output is a kmeans object giving the clusters.
#' @export
spectral_cluster <- function(g, sigma = 1, k = 2, nreps = 100){

	A <- affinity_matrix(g, sigma = sigma)

	L   <- generalized_laplacian_matrix(A, sigma)
	evL <- eigen(L, symmetric = T)
	Z <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]

	models <- data_frame(n = 1:nreps) %>%
		dplyr::mutate(model = map(n, ~ kmeans(Z, centers = k))) %>%
		mutate(perf = map_dbl(model, ~.$tot.withinss))

	model <- models %>%
		filter(perf == min(perf))

	df <- data_frame(tract = row.names(A), cluster = model$model[[1]]$cluster)
	g <- g %>% set_vertex_attr('cluster', index = df$tract, value = df$cluster)
	g
}


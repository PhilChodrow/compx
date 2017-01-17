#' tractwise_clustering
#' @name tractwise_clustering
#' @docType package
#' @import tidyverse sp maptools rgeos
NULL


#' Construct a list of edges for an adjacency graph from an spdf
#' @param tracts an spdf
#' @param ... additional arguments passed to id_lookup
#' @return a tibble with two columns, with each row corresponding to an edge
#' or adjacency. This tibble contains no duplicates (i.e. is suitable for
#' undirected) graphs.

list_edges <- function(tracts, ...){

	ids <- id_lookup(tracts, ...)

	gRelate(tracts, byid = TRUE, pattern = '****1****') %>%
		as.data.frame() %>%
		rownames_to_column('id_1') %>%
		tbl_df() %>%
		gather(key = id_2, value = adj, -id_1) %>%
		filter(adj) %>%
		select(-adj) %>%
		mutate(key = paste(pmax(id_1, id_2), pmin(id_1, id_2))) %>%
		filter(!duplicated(key),
			   id_1 != id_2) %>%
		select(-key) %>%
		left_join(ids, by = c('id_1' = 'id')) %>%
		left_join(ids, by = c('id_2' = 'id'), suffix = c('_1', '_2')) %>%
		select(-id_1, -id_2)
}


#' Construct an edge list with nested demographic and coordinate data
#' @param tracts an spdf
#' @param data the tibble containing demographic information
#' @param ... additional arguments passed to id_lookup
#' @return a tibble containing node-level demographic summary information.
edge_list_with_data <- function(tracts, data, ...){

	edge_list <- list_edges(tracts,...)

	join_data <- data %>%
		group_by(tract) %>%
		ungroup() %>%
		nest(-tract)

	edge_list %>%
		left_join(join_data, by = c('geoid_1' = 'tract')) %>%
		left_join(join_data, by = c('geoid_2' = 'tract'), suffix = c('_1', '_2')) %>%
		mutate(missing = map2_lgl(data_1, data_2, ~ is.null(.x) | is.null(.y))) %>%
		filter(!missing) %>%
		select(-missing) %>%
		mutate(n_1  = map(data_1, ~.$n),
			   n_2  = map(data_2, ~.$n),
			   n_12 = map2(n_1, n_2, ~.x + .y),
			   p_1  = map(n_1, ~ ./sum(.)),
			   p_2  = map(n_2, ~ ./sum(.)),
			   p_12 = map(n_12, ~./sum(.)))
}


#' Construct an edge list with information-geodesic distance between nodes.
#' @param tracts an spdf
#' @param data a tibble containing demographic information keyed to tracts
#' @param divergence one of c('KL', 'euc', 'cum_euc', 'js') giving the desired
#' divergence measure between probability distributions
#' @param ... additional arguments passed to id_lookup
#' @return an edge list as a tibble with geodesic distances between pairs
#' of nodes.

edge_list_with_dists <- function(tracts, data, divergence = 'KL', ...){
	edge_list_df <- edge_list_with_data(tracts, data, ...)

	jensen_shannon <- function(n_1,n_2, weighted = FALSE){
		tryCatch({

			n_12 = n_1 + n_2

			r = n_12 / sum(n_12)
			p = n_1 / sum(n_1)
			q = n_2 / sum(n_2)

			if(weighted){
				return((sum(n_1) / sum(n_12)) * compx::DKL(p, r) +
					   	(sum(n_2) / sum(n_12)) * compx::DKL(q,r))
			}else{
				return(compx::DKL(p, (p+q) / 2) + compx::DKL(q, (p+q) / 2))
			}

		}, error = function(e) NA)
	}


	cum_euc <- function(p){
		lower <- matrix(0, length(p), length(p))
		lower[lower.tri(lower, diag = T)] <- 1
		lower %*% t(lower)
	}

	hessians <- list(KL      = function(p) diag(1/p),
					 euc     = function(p) diag(length(p)),
					 cum_euc = cum_euc)


	out_df <- edge_list_df
	if(divergence %in% c('KL', 'euc', 'cum_euc')){
		out_df <- out_df %>%
			mutate(H = map(p_12, hessians[[divergence]]),
				   d_alpha = map2(p_1, p_2, ~ .x - .y),
				   sqdist = map2_dbl(d_alpha, H, NA_multiply),
				   dist = sqrt(sqdist))
	}else{
		out_df <- out_df %>%
			mutate(js = map2_dbl(n_1, n_2, jensen_shannon, weighted = FALSE),
				   dist = sqrt(js))
	}
	out_df %>%
		filter(!is.na(dist)) %>%
		select(geoid_1, geoid_2, dist)
}

#' Construct an igraph::graph from an spdf with demographic data
#' @param tracts an spdf
#' @param data a tibble containing demographics keyed to tracts
#' @param divergence one of c('KL', 'euc', 'cum_euc', 'js') giving the desired
#' divergence measure between probability distributions
#' @param ... additional arguments passed to id_lookup
#' @return g an igraph::graph with geodesic distance as the dist edge attribute
#' @export

graph_from_tracts <- function(tracts, data, divergence = 'KL', ...){

	edge_list <- edge_list_with_dists(tracts, data, divergence)
	coords <- coords_df(tracts, ...) %>%
		filter(geoid %in% edge_list$geoid_1 | geoid %in% edge_list$geoid_2)

	g <- edge_list %>%
		graph_from_data_frame(directed = FALSE) %>%
		set.vertex.attribute(name = 'x', index = coords$geoid, value = coords$x) %>%
		set.vertex.attribute(name = 'y', index = coords$geoid, value = coords$y)
}

#' Construct an affinity matrix for a graph g based on the dist edge attribute
#' @param g an igraph::graph with a dist edge attribute/
#' @param sigma the multiplicative constant in the quadratic exponential
#' affinity kernel
#' @return a matrix containing the affinities computed with the quadratic exponential
#' affinity kernel
#' @export

affinity_matrix <- function(g, sigma = 1){
	dists <- distances(g, weights = E(g)$dist^2)
	exp(-dists * sigma)
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






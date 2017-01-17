#' tractwise geometry
#' @name tractwise_geometry
#' @docType package
#' @import tidyverse sp maptools rgeos
NULL

#' Construct a lookup data frame from an SPDF object
#' @param tracts the object whose keys to retrieve
#' @param key_col the column of tracts@data containing the identifier.
#' @return df a two-column tibble in which the first column is the id and
#' the second is the lookup identifier.
#' @export
id_lookup <- function(tracts, key_col = 'GEOID'){
	data_frame(id = map_chr(tracts@polygons, ~.@ID),
			   geoid = tracts@data[,key_col])
}

#' Construct a tibble of geographic distances between pairs of tract centroids.
#' @param tracts an spdf
#' @param ... additional parameters passed to id_lookup
#' @return a tibble with three columns: two are tract names, and the third
#' is the geographic distance between them.
distance_df <- function(tracts, ...){
	ids <- id_lookup(tracts, ...)

	centroids <- gCentroid(tracts,byid=TRUE)

	dist(centroids@coords) %>%
		as.matrix() %>%
		as.data.frame() %>%
		rownames_to_column('id_1') %>%
		tbl_df() %>%
		gather(key = id_2, value = distance, -id_1) %>%
		left_join(ids, by = c('id_1' = 'id')) %>%
		left_join(ids, by = c('id_2' = 'id'), suffix = c('_1', '_2')) %>%
		select(geoid_1, geoid_2, distance)
}

#' Construct a tibble of centroid coordinates
#' @param tracts the spdf whose coordinates we wish to compute.
#' @param ... additional parameters passed to id_lookup
#' @export
coords_df <- function(tracts, ...){
	ids <- id_lookup(tracts, ...)

	centroids <- gCentroid(tracts,byid=TRUE)

	centroids %>%
		as.data.frame() %>%
		rownames_to_column('id') %>%
		tbl_df() %>%
		left_join(ids, by = c('id' = 'id')) %>%
		select(-id, x, y)
}

#' Construct a tibble of tracts that are adjacent to each other.
#' @param tracts an spdf
#' @param ... additional parameters passed to id_lookup
#' @return tibble with two columns; two tract ids appear next to each
#' other if they are geographically adjacent.
#' @export
adjacency_df <- function(tracts, ...){
	ids <- id_lookup(tracts, ...)

	dists <- distance_df(tracts, ...)

	gRelate(tracts, byid = TRUE, pattern = '****1****') %>%
		as.data.frame() %>%
		rownames_to_column('id_1') %>%
		tbl_df() %>%
		gather(key = id_2, value = adj, -id_1) %>%
		filter(adj) %>%
		select(-adj) %>%
		left_join(ids, by = c('id_1' = 'id')) %>%
		left_join(ids, by = c('id_2' = 'id'), suffix = c('_1', '_2')) %>%
		select(-id_1, -id_2) %>%
		left_join(dists, by = c('geoid_1' = 'geoid_1', 'geoid_2' = 'geoid_2'))
}

#' Construct a (nested) tibble with tract-level demographics and coordinates
#' @param tracts the spdf
#' @param data the tibble containing the demographics for the tracts.
#' @return a nested tibble containing tract-level demographics and coordinates
#' as vectors in list-columns.
node_df <- function(tracts, data, ...){

	coords <- coords_df(tracts, ...)

	data %>%
		select(-t) %>%
		group_by(tract) %>%
		mutate(p = n / sum(n)) %>%
		ungroup() %>%
		select(-group) %>%
		nest(-tract, .key = 'counts') %>%
		mutate(total = map_dbl(counts, ~ sum(.$n)),
			   p     = map(counts, ~ .$p)) %>%
		select(-counts) %>%
		left_join(coords, by = c('tract' = 'geoid')) %>%
		mutate(coords = map2(x, y, ~ c(.x,.y))) %>%
		select(-x, -y)
}

#' Compute numerical derivatives from data using weighted linear regression
#' on each tract's 1-ego network.
#' @param tracts, the spdf
#' @param data, a tibble containing the demographic data
#' @param sigma the multiplicative constant in the square exponential kernel.
#' @return a tract-level tibble with a column D_alpha giving the derivatives of
#' demographic frequencies with respect to spatial coordinates.
compute_derivatives <- function(tracts, data, sigma = 1, ...){

	do_regression <- function(X, Y, W){

		(solve((t(X) %*% W) %*% X) %*% t(X)) %*% (W %*% Y)
	}

	node_data      <- node_df(tracts, data, ...)

	adjacency_data <- adjacency_df(tracts, ...)

	adjacency_data %>%
		left_join(node_data, by = c('geoid_1' = 'tract')) %>%
		left_join(node_data, by = c('geoid_2' = 'tract'), suffix = c('_1', '_2')) %>%
		mutate(X = map2(coords_2, coords_1, ~ .x - .y),
			   Y = map2(p_1, p_2, ~.x - .y)) %>%
		mutate(regression_weight = exp(- sigma * distance^2)) %>%
		select(geoid_1, geoid_2, regression_weight, X, Y) %>%
		group_by(geoid_1) %>%
		filter(n() > 2) %>%
		do(X = reduce(.$X, rbind),
		   Y = reduce(.$Y, rbind),
		   W = reduce(.$regression_weight, .Primitive('c'))) %>%
		ungroup() %>%
		mutate(W = map(W, ~ diag(., nrow = length(.)))) %>%
		mutate(D_alpha = pmap(list(X, Y, W), do_regression))
}

#' Compute the hessian (or Riemannian metric) at each tract.
#' @param tracts an spdf
#' @param data a tibble containing demographic data
#' @param divergence one of c('KL', 'euc', 'cum_euc')
#' @param sigma multiplicative constant on the square exponential kernel weights
#' @return a tibble with a column local_H containing the Hessian
#' (or Riemannian metric) in local coordinates.
#' @export
compute_hessian <- function(tracts, data, divergence = 'KL', sigma = 1, ...){

	cum_euc <- function(p){
		lower <- matrix(0, length(p), length(p))
		lower[lower.tri(lower, diag = T)] <- 1
		lower %*% t(lower)
	}

	hessians <- list(KL      = function(p) diag(1/p),
					 euc     = function(p) diag(length(p)),
					 cum_euc = cum_euc)

	derivs_df <- compute_derivatives(tracts, data, sigma, ...)
	node_data <- node_df(tracts, data, ...)

	node_data %>%
		left_join(derivs_df, by = c('tract' = 'geoid_1')) %>%
		mutate(missing = map_lgl(D_alpha, is.null)) %>%
		filter(!missing) %>%
		mutate(H  = map(p,  hessians[[divergence]]),
			   local_H = map2(D_alpha, H, ~ .x %*% .y %*% t(.x))) %>%
		select(tract, total, local_H)
}

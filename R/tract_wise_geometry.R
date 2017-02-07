#' tractwise geometry
#' @name tractwise_geometry
#' @docType package
#' @import tidyverse sp maptools rgeos
NULL


#
# #' Construct a tibble of geographic distances between pairs of tract centroids.
# #' @param tracts an spdf
# #' @param ... additional parameters passed to id_lookup
# #' @return a tibble with three columns: two are tract names, and the third
# #' is the geographic distance between them.
# distance_df <- function(tracts, ...){
# 	ids <- id_lookup(tracts, ...)
#
# 	centroids <- gCentroid(tracts,byid=TRUE)
#
# 	dist(centroids@coords) %>%
# 		as.matrix() %>%
# 		as.data.frame() %>%
# 		rownames_to_column('id_1') %>%
# 		tbl_df() %>%
# 		gather(key = id_2, value = distance, -id_1) %>%
# 		left_join(ids, by = c('id_1' = 'id')) %>%
# 		left_join(ids, by = c('id_2' = 'id'), suffix = c('_1', '_2')) %>%
# 		select(geoid_1, geoid_2, distance)
# }
#
#' Construct a tibble of centroid coordinates
#' @param tracts the spdf whose coordinates we wish to compute.
#' @param ... additional parameters passed to id_lookup
#' @export
coords_df <- function(tracts, km = FALSE, ...){
	ids <- id_lookup(tracts, ...)

	centroids <- gCentroid(tracts,byid=TRUE)

	centroids <- centroids %>%
		as.data.frame() %>%
		rownames_to_column('id') %>%
		tbl_df() %>%
		left_join(ids, by = c('id' = 'id')) %>%
		select(-id, x, y)

	if(km){
		centroids <- centroids %>%
			mutate(x = x * cos(y / 360) * 111,
				   y = y * 111)
	}
	return(centroids)
}
#
# #' Construct a tibble of tracts that are adjacent to each other.
# #' @param tracts an spdf
# #' @param ... additional parameters passed to id_lookup
# #' @return tibble with two columns; two tract ids appear next to each
# #' other if they are geographically adjacent.
# #' @export
# adjacency_df <- function(tracts, ...){
# 	ids <- id_lookup(tracts, ...)
#
# 	dists <- distance_df(tracts, ...)
#
# 	gRelate(tracts, byid = TRUE, pattern = '****1****') %>%
# 		as.data.frame() %>%
# 		rownames_to_column('id_1') %>%
# 		tbl_df() %>%
# 		gather(key = id_2, value = adj, -id_1) %>%
# 		filter(adj) %>%
# 		select(-adj) %>%
# 		left_join(ids, by = c('id_1' = 'id')) %>%
# 		left_join(ids, by = c('id_2' = 'id'), suffix = c('_1', '_2')) %>%
# 		select(-id_1, -id_2) %>%
# 		left_join(dists, by = c('geoid_1' = 'geoid_1', 'geoid_2' = 'geoid_2'))
# }

#' Construct a (nested) tibble with tract-level demographics and coordinates
#' @param tracts the spdf
#' @param data the tibble containing the demographics for the tracts.
#' @return a nested tibble containing tract-level demographics and coordinates
#' as vectors in list-columns.
adj_with_coords <- function(adj, tracts, km = FALSE, ...){

	coords <- coords_df(tracts, km, ...)

	adj <- adj %>%
		left_join(coords, by = c('geoid_1' = 'geoid')) %>%
		left_join(coords, by = c('geoid_2' = 'geoid'), suffix = c('_1', '_2'))

	if('t_1' %in% names(adj)){
		adj <- adj %>%
			mutate(coords_1 = pmap(list(x_1, y_1, t_1), c),
				   coords_2 = pmap(list(x_2, y_2, t_2), c))
	}else{
		adj <- adj %>%
			mutate(coords_1 = pmap(list(x_1, y_1), c),
				   coords_2 = pmap(list(x_2, y_2), c))
	}
	return(adj)
}

adj_with_geo_distance <- function(adj, tracts,...){

	ids <- id_lookup(tracts, ...)

	centroids <- gCentroid(tracts,byid=TRUE)

	dist_df <- dist(centroids@coords) %>%
		as.matrix() %>%
		as.data.frame() %>%
		rownames_to_column('id_1') %>%
		tbl_df() %>%
		gather(key = id_2, value = distance, -id_1) %>%
		left_join(ids, by = c('id_1' = 'id')) %>%
		left_join(ids, by = c('id_2' = 'id'), suffix = c('_1', '_2')) %>%
		select(geoid_1, geoid_2, distance)

	adj %>%
		left_join(dist_df, by = c('geoid_1' = 'geoid_1', 'geoid_2' = 'geoid_2'))
}



#' Compute numerical derivatives from data using weighted linear regression
#' on each tract's 1-ego network.
#' @param tracts, the spdf
#' @param data, a tibble containing the demographic data
#' @param sigma the multiplicative constant in the square exponential kernel.
#' @return a tract-level tibble with a column D_alpha giving the derivatives of
#' demographic frequencies with respect to spatial coordinates.
compute_derivatives <- function(adj, sigma = 1, ...){

	do_regression <- function(X, Y, W = diag(dim(X)[2])){

		(solve((t(X) %*% W) %*% X) %*% t(X)) %*% (W %*% Y)
	}

	if('t_1' %in% names(adj)){
		group_vars <- c('geoid_1', 't_1')
	}else{
		group_vars <- c('geoid_1')
	}

	select_vars <- group_vars %>% append(c('geoid_2', 'regression_weight', 'X', 'Y'))

	adj <- adj %>%
		mutate(X = map2(coords_2, coords_1, ~ .x - .y),
			   p_1 = map(n_1, simplex_normalize),
			   p_2 = map(n_2, simplex_normalize),
			   Y = map2(p_1, p_2, ~.x - .y))

	out_df <- adj %>%
		mutate(regression_weight = exp(- sigma * distance^2)) %>%
		select_(.dots = select_vars) %>%
		group_by_(.dots = group_vars) %>%
		filter(n() > 2) %>%
		do(X = reduce(.$X, rbind),
		   Y = reduce(.$Y, rbind),
		   W = reduce(.$regression_weight, .Primitive('c'))) %>%
		ungroup() %>%
		mutate(W = map(W, ~ diag(., nrow = length(.)))) %>%
		mutate(D_alpha = pmap(list(X, Y, W), do_regression))

	smoother <- adj %>%
		mutate(regression_weight = exp(- sigma * distance^2)) %>%
		select_(.dots = group_vars %>% append(c('n_2', 'regression_weight'))) %>%
		group_by_(.dots = group_vars) %>%
		mutate(n_2 = map2(n_2, regression_weight, ~ .y * .x)) %>%
		do(smoothed_n = reduce(.$n_2, `+`)) %>%
		ungroup()

	lookup_df <- adj %>%
		select_(.dots = group_vars %>% append(c('p_1', 'n_1'))) %>%
		distinct_(.dots = group_vars, .keep_all = TRUE) %>%
		mutate(total = map_dbl(n_1, sum)) %>%
		select(-n_1)

	names(group_vars) <- group_vars
	out_df %>%
		left_join(lookup_df, by = group_vars) %>%
		left_join(smoother, by = group_vars)
}





#' Compute the hessian (or Riemannian metric) at each tract.
#' @param tracts an spdf
#' @param data a tibble containing demographic data
#' @param divergence one of c('KL', 'euc', 'cum_euc')
#' @param sigma multiplicative constant on the square exponential kernel weights
#' @return a tibble with a column local_H containing the Hessian
#' (or Riemannian metric) in local coordinates.
#' @export
compute_hessian <- function(adj, hessian = DKL_, smooth = F){

	select_vars <- c('geoid', 'total', 'local_g')
	adj <- adj %>% rename(geoid = geoid_1)

	if('t_1' %in% names(adj)){
		 select_vars <- select_vars %>% append('t')
		 adj <- adj %>% rename(t = t_1)
	}

	if(smooth){
		out <- adj %>%
			mutate(smoothed_p = map(smoothed_n, ~ . / sum(.)),
				   H  = map(smoothed_p,  hessian),
				   local_g = map2(D_alpha, H, ~ NA_multiply(t(.x), .y)))
	}else{
		out <- adj %>%
			mutate(H  = map(p_1,  hessian),
				   local_g = map2(D_alpha, H, ~ NA_multiply(t(.x), .y)))
	}
	out %>%
		select_(.dots = select_vars)
}

#'
#' @export
compute_metric <- function(tracts, data, km = T, sigma = 100, hessian = euc_, smooth = F){
	out <- tracts[tracts@data$GEOID %in% data$tract,] %>%
		make_adjacency()

	if('t' %in% names(data)){
		out <- out %>% add_temporal(unique(data$t))
	}
	adj <- out %>%
		add_data(data) %>%
		adj_with_coords(tracts, km = km) %>%
		adj_with_geo_distance(tracts) %>%
		compute_derivatives(sigma = sigma) %>%
		compute_hessian(hessian = hessian, smooth = smooth)
}

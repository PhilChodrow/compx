#' tractwise geometry
#' @name geometry
#' @docType package
#' @import tidyverse sp maptools rgeos units sf
NULL

#' Construct a tibble of centroid coordinates from an sf data.frame containing polygons in its geometry column.
#' @param tracts the sf data.frame whose coordinates we wish to compute.
#' @param km boolean should the centroids be converted so that the euclidean distances between them are approximately measurable in km?
#' @export
coords_df <- function(tracts, km = FALSE, ...){

	centroids <- st_centroid(tracts) %>%
		mutate(x = map_dbl(geometry, ~.[1]),
			   y = map_dbl(geometry, ~.[2])) %>%
		tbl_df() %>%
		select(GEOID, x, y) %>%
		rename(geoid = GEOID)

	if(km){
		centroids <- centroids %>%
			mutate(x = x * cos(y / 360) * 111,
				   y = y * 111)
	}
	return(centroids)
}

#' Construct a (nested) tibble with tract-level demographics and coordinates
#' @param tracts the sf data.frame supplying geographic information
#' @param data the tibble containing the demographics for the tracts.
#' @param km boolean should the centroids be converted so that the euclidean distances between them are approximately measurable in km?
#' @return a nested tibble containing tract-level demographics and coordinates
#' as vectors in list-columns.

adj_with_coords <- function(adj, tracts, km = FALSE){

	coords <- coords_df(tracts, km)

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

#' Compute numerical derivatives from data using weighted linear regression
#' on each tract's 1-ego network.
#' @param adj the adjacency matrix on which to compute. Must include centroids.
#' @param r_sigma geographic bandwidth of the regression weight decay
#' @param s_sigma if smooth = TRUE, the bandwidth of a geographic rbf smoother
#' @param smooth boolean should the derivatives be computed on a geographically smoothed demographic distribution?
#' @return a tract-level tibble with a column D_alpha giving the derivatives of
#' demographic frequencies with respect to spatial coordinates.
compute_derivatives <- function(adj, r_sigma = 1, s_sigma = 1, smooth = F){

	do_regression <- function(X, Y, W = diag(dim(X)[2])){
		tryCatch({
			(solve((t(X) %*% W) %*% X) %*% t(X)) %*% (W %*% Y)
			},
			error = function(e) matrix(NA, dim(X)[1], dim(Y)[2]) )
	}

	if('t_1' %in% names(adj)){
		group_vars <- c('geoid_1', 't_1')
		group_vars_renamed <- c('geoid', 't')
	}else{
		group_vars <- c('geoid_1')
		group_vars_renamed <- c('geoid')
	}

	select_vars <- group_vars %>% append(c('geoid_2', 'regression_weight', 'X', 'Y'))

	adj <- adj %>%
		mutate(X = map2(coords_2, coords_1, ~ .x - .y),
			   p_1 = map(n_1, simplex_normalize),
			   p_2 = map(n_2, simplex_normalize),
			   Y = map2(p_1, p_2, ~.x - .y))

	out_df <- adj %>%
		mutate(regression_weight = exp(- distance^2 / (2*r_sigma))) %>%
		select_(.dots = select_vars) %>%
		group_by_(.dots = group_vars) %>%
		filter(n() > 2) %>%
		do(X = reduce(.$X, rbind),
		   Y = reduce(.$Y, rbind),
		   W = reduce(.$regression_weight, .Primitive('c'))) %>%
		ungroup() %>%
		mutate(W = map(W, ~ diag(., nrow = length(.)))) %>%
		mutate(D_alpha = pmap(list(X, Y, W), do_regression))

	names(group_vars) <- group_vars
	lookup_df <- adj %>%
		select_(.dots = group_vars %>% append(c('p_1', 'n_1'))) %>%
		distinct_(.dots = group_vars, .keep_all = TRUE) %>%
		mutate(total = map_dbl(n_1, sum)) %>%
		select(-n_1)
	out_df <- out_df %>%
		left_join(lookup_df, by = group_vars)

	if(smooth){
		smoother <- adj %>%
			mutate(smoothing_weight = exp(- distance^2 / (2*s_sigma))) %>%
			select_(.dots = group_vars %>% append(c('n_2', 'smoothing_weight'))) %>%
			group_by_(.dots = group_vars) %>%
			mutate(n_2 = map2(n_2, smoothing_weight, ~ .y * .x)) %>%
			do(smoothed_n = reduce(.$n_2, `+`)) %>%
			ungroup()
		out_df <- out_df %>%
			left_join(smoother, by = group_vars) %>%
			mutate(p_1 = map(smoothed_n, ~ . / sum(.)))
	}
	if('t_1' %in% names(out_df)){
		out_df <- out_df %>%
			rename(t = t_1)
	}
	out_df %>%
		rename(geoid = geoid_1, p = p_1) %>%
		select_(.dots = group_vars_renamed %>% append(c('p', 'total', 'D_alpha')))
}

#' Compute the hessian (or Riemannian metric) at each tract.
#' @param adj an adjacency tibble in which derivatives have already been computed
#' @param hessian one of `c(DKL_, euc_, cum_euc_)`
#' @return a tibble with a column local_H containing the Hessian
#' (or Riemannian metric) in local coordinates.
compute_hessian <- function(adj, hessian = DKL_){

	select_vars <- c('geoid', 'total', 'g')

	if('t' %in% names(adj)){
		select_vars <- select_vars %>% append('t')
	}

	out <- adj %>%
		mutate(H  = map(p,  hessian),
			   g = map2(D_alpha, H, ~ NA_multiply(t(.x), .y)))
	out %>%
		select_(.dots = select_vars)
}

#' Compute the pullback metric tensor from an sf data frame and an accompanying set of
#' demographic data.
#'
#' @param tracts the sf data.frame supplying geographic information
#' @param data the tibble supplying demographic data
#' @param km whether to compute rectangular distances in km
#' @param r_sigma regression smoother for derivative computation
#' @param s_sigma geographic smoother for derivative computation
#' @param smooth boolean should a geographic smoother be used for derivatives?
#' @param hessian the hessian function to be used, one of c(DKL_, euc_, cum_euc_)
#'
#' @return a tibble in which each row has a list-column containing the pullback metric tensor stored as a matrix.
#' @export

compute_metric <- function(tracts, data, km = T, r_sigma = 100, s_sigma = 1, smooth = F, hessian = euc_){


	assertthat::assert_that(st_crs(tracts)$epsg == 4326, msg = "compx expects lon-lat projected data (crs = 4326). Please try tracts <- st_transform(tracts, 4326).")



	out <- tracts %>%
		filter(GEOID %in% data$tract) %>%
		make_adjacency()

	if('t' %in% names(data)){
		out <- out %>% add_temporal(unique(data$t))
	}
	adj <- out %>%
		add_data(data) %>%
		adj_with_coords(tracts, km = km) %>%
		mutate(distance = sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)) %>%
		compute_derivatives(r_sigma = r_sigma, s_sigma = s_sigma, smooth = smooth) %>%
		compute_hessian(hessian = hessian)
	adj
}

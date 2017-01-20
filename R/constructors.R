#' constructors
#' @name constructors
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
	data_frame(id = purrr::map_chr(tracts@polygons, ~.@ID),
			   geoid = tracts@data[,key_col]) %>%
		mutate(geoid = as.character(geoid))
}

#' make a thing
#' @param tracts the spdf
#' @export
make_adjacency <- function(tracts, ...){
	ids <- id_lookup(tracts, ...)

	gRelate(tracts, byid = TRUE, pattern = '****1****') %>%
		as.data.frame() %>%
		rownames_to_column('id_1') %>%
		tbl_df() %>%
		gather(key = id_2, value = adj, -id_1) %>%
		filter(adj) %>%
		select(-adj) %>%
		left_join(ids, by = c('id_1' = 'id')) %>%
		left_join(ids, by = c('id_2' = 'id'), suffix = c('_1', '_2')) %>%
		select(-id_1, -id_2)
}

# this is testable via simple counting
add_temporal <- function(adj, t_vec = NULL){

	# first construct layers for each value of t
	space_cxns <- t_vec %>%
		map(~ mutate(adj, t_1 = .)) %>%
		reduce(rbind) %>%
		mutate(t_2 = t_1)
	# then add connections through time
	time_cxns <- expand.grid(geoid_1 = unique(adj$geoid_1),
				t_1 = t_vec,
				stringsAsFactors = FALSE) %>%
		tbl_df() %>%
		group_by(geoid_1) %>%
		mutate(t_2a = lead(t_1),
			   t_2b = lag(t_1)) %>%
		ungroup() %>%
		gather(time, t_2, -geoid_1, -t_1) %>%
		select(-time) %>%
		mutate(geoid_2 = geoid_1) %>%
		select(geoid_1, geoid_2, t_1, t_2) %>%
		filter(!is.na(t_1), !is.na(t_2))

	return(rbind(space_cxns, time_cxns))

}

undirect <- function(adj, allow_self_loops = FALSE){

	# construct spatial keys
	out <- adj %>%
		mutate(key = paste(pmax(geoid_1, geoid_2),
						   pmin(geoid_1, geoid_2)))

	# construct temporal keys if there are time columns in adj
	if(('t_1' %in% names(adj))){
		out <- out %>%
			mutate(t_key = paste(pmax(t_1, t_2), pmin(t_1, t_2)),
				   key = paste(key, t_key)) %>%
			select(-t_key)
	}

	# remove spatial (and temporal) duplicates
	out <- out %>%
		filter(!duplicated(key))

	# remove self loops if option is toggled
	if(!allow_self_loops){
		out <- out %>%
			filter(!(geoid_1 == geoid_2))
		if('t_1' %in% names(adj)){
			out <- out %>%
				filter(t_1 == t_2)
		}
	}

	out %>% select(-key)
}

add_data <- function(adj, data){

	all_ids <- c(adj$geoid_1, adj$geoid_2) %>% unique()
	join_data <- data %>%
		filter(tract %in% all_ids) %>%
		nest_(key_col = 'data', nest_cols = c('group', 'n'))

	# join depending on whether there are temporal columns
	if('t_1' %in% names(adj)){
		join_keys_1   <- c('geoid_1' = 'tract', 't_1' = 't')
		join_keys_2   <- c('geoid_2' = 'tract', 't_2' = 't')
	}else{
		join_keys_1 <- c('geoid_1' = 'tract')
		join_keys_2 <- c('geoid_2' = 'tract')
	}

	# add information from data df, and construct some useful numerical columns
	# for later computations.
	adj %>%
		left_join(join_data, by = join_keys_1) %>%
		left_join(join_data, by = join_keys_2, suffix = c('_1', '_2')) %>%
		mutate(n_1  = map(data_1, ~.$n),
			   n_2  = map(data_2, ~.$n),
			   n_12 = map2(n_1, n_2, ~.x + .y),
			   p_1  = map(n_1, ~ ./sum(.)),
			   p_2  = map(n_2, ~ ./sum(.)),
			   p_12 = map(n_12, ~./sum(.)))
}

add_hessian <- function(adj, divergence = 'DKL', p = 'p_12'){
	assertthat::assert_that(divergence %in% c('DJS', 'DKL', 'euc'))
	cum_euc <- function(p){
		lower <- matrix(0, length(p), length(p))
		lower[lower.tri(lower, diag = T)] <- 1
		lower %*% t(lower)
	}

	hessians <- list(DKL      = function(p) diag(1/p),
					 euc     = function(p) diag(length(p)),
					 cum_euc = cum_euc)

	adj %>%
		mutate(H = map(p_12, hessians[[divergence]]))
}


#' Early stage: pass hessians as names, pass DJS as a function.
information_distances <- function(adj, divergence, ...){

	divergences <- list(DKL = NULL,
						euc = NULL,
						cum_euc = NULL,
						DJS = DJS)

	assertthat::assert_that(divergence %in% names(divergences))

	NA_multiply <- function(D_alpha, H){
		na_locs           <- is.na(t(D_alpha) %*% D_alpha)
		H[is.infinite(H)] <- 0
		result            <- t(D_alpha) %*% H %*% D_alpha
		result[na_locs]   <- NA
		result
	}

	if(divergence %in% c('DKL', 'euc', 'cum_euc')){
		adj <- adj %>%
			add_hessian(divergence, p = 'p_12') %>%
			mutate(d_alpha = map2(p_1, p_2, ~ .x - .y),
				   sqdist = map2_dbl(d_alpha, H, NA_multiply),
				   dist = sqrt(sqdist))
	}else{
		adj <- adj %>%
			mutate(dist = map2_dbl(p_1, p_2, ~ divergences[[divergence]](.x, .y, ...)))
	}
	return(adj)
}

#' @export
make_graph <- function(adj){

	if('t_1' %in% names(adj)){
		g <- adj %>%
			mutate(key_1 = paste(geoid_1, t_1, sep = '_'),
				   key_2 = paste(geoid_2, t_2, sep = '_')) %>%
			select(key_1, key_2, dist) %>%
			igraph::graph_from_data_frame(directed = FALSE)

		g <- g %>%
			set.vertex.attribute(name = 'geoid',
								 value = names(V(g)) %>% stringr::str_sub(end = -6)) %>%
			set.vertex.attribute(name = 't',
								 value = names(V(g)) %>% stringr::str_sub(start = -4) %>% as.numeric())
		return(g)
	}else{
		adj %>%
			select(geoid_1, geoid_2, dist) %>%
			igraph::graph_from_data_frame()

	}
}

add_coordinates <- function(g, tracts, ...){
	ids <- id_lookup(tracts, ...)

	centroids <- gCentroid(tracts,byid=TRUE)

	coords <- centroids %>%
		as.data.frame() %>%
		rownames_to_column('id') %>%
		tbl_df() %>%
		left_join(ids, by = c('id' = 'id')) %>%
		select(-id, x, y)

	if(is.null(V(g)$geoid)){
		g <- g %>%
			set.vertex.attribute(name = 'x', index = coords$geoid, value = coords$x) %>%
			set.vertex.attribute(name = 'y', index = coords$geoid, value = coords$y)
		return(g)
	}else{
	coords <- data_frame(geoid = V(g)$geoid, key = names(V(g))) %>%
		left_join(coords, by = c('geoid' = 'geoid'))

	g <- g %>%
		set.vertex.attribute(name = 'x', index = coords$key, value = coords$x) %>%
		set.vertex.attribute(name = 'y', index = coords$key, value = coords$y)
	return(g)

	}


}

#' @export
construct_information_graph <- function(tracts, data, divergence = 'euc'){
	adj <- tracts[tracts@data$GEOID %in% data$tract,] %>%
		make_adjacency()

	if('t' %in% names(data)){
		adj <- adj %>%
			add_temporal(unique(data$t))
	}
	adj %>%
		undirect() %>%
		add_data(data) %>%
		information_distances(divergence = 'DKL', w = .5) %>%
		make_graph() %>%
		add_coordinates(tracts)
}



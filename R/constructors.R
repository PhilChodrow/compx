#' constructors
#' @name constructors
#' @docType package
#' @import tidyverse sp maptools rgeos igraph units sf
NULL

#' Construct a lookup data frame from an SPDF object
#' @param tracts the object whose keys to retrieve
#' @param key_col the column of tracts@data containing the identifier.
#' @return df a two-column tibble in which the first column is the id and
#' the second is the lookup identifier.
#' @export
id_lookup <- function(tracts, key_col = 'GEOID'){
	tracts[[key_col]] %>%
		data_frame(row = as.character(1:length(.)), geoid = .)
}

#' Construct an adjacency tibble from an sf data.frame of polygons
#' @param tracts an sf data.frame in which the `geometry` column holds polygons
#' @return df a two-column tibble in which two ids appear in the same row iff
#' their corresponding polygons are related according to the '****1****' DE9-IM pattern. Learn more about these patterns at the [wikipedia page](https://en.wikipedia.org/wiki/DE-9IM).
#' @export
make_adjacency <- function(tracts){
	lookup  <- id_lookup(tracts)
	adj_mat <- st_relate(tracts, pattern = '****T****', sparse = TRUE) # as sparse list
	1:length(adj_mat) %>%
		map(~data_frame(from = as.character(.),
						to = as.character(adj_mat[[.]]))) %>%
		reduce(rbind) %>%
		left_join(lookup, by = c('from' = 'row')) %>%
		left_join(lookup, by = c('to' = 'row'), suffix = c('_1', '_2')) %>%
		select(-from, -to)
}

#' Add temporal data to a nontemporal adjacency df
#' @param adj the nontemporal df to which we add a temporal column
#' @param t_vec a list of the timestamps to add.
#' @export
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
		if('t_1' %in% names(out)){
			out <- out %>%
				filter(geoid_1 != geoid_2 | t_1 != t_2)
		}else{
			out <- out %>%
				filter(geoid_1 != geoid_2)
		}
	}

	out %>% select(-key)
}

#' @export
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
		mutate(n_1  = map(data_1,    ~.$n),
			   n_2  = map(data_2,    ~.$n))
}


#' We are assuming a divergence function, defined by the user,
#' that takes in n_1 and n_2 and returns a number.

information_distances <- function(adj, divergence){

	adj %>%
		mutate(dist = map2_dbl(n_1, n_2, divergence),
			   dist = ifelse(is.na(dist), Inf, dist)) # bit of a hack
}


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
			igraph::graph_from_data_frame(directed = FALSE)

	}
}


add_n <- function(h, data_df){

	lookup <- data_df %>%
		tidyr::nest_(key_col = 'n_col', nest_cols = c('group', 'n')) %>%
		dplyr::mutate(n = purrr::map(n_col, ~.$n))

	valid_names <- igraph::V(h)$name

	if('t' %in% names(lookup)){
		lookup <- lookup %>%
			unite(col = key, tract, t)
	}else{
		lookup <- lookup %>%
			mutate(key = tract)
	}

	lookup <- lookup %>%
		filter(key %in% valid_names)

	h %>% igraph::set.vertex.attribute('n', index = lookup$key, value = lookup$n)
}



add_coordinates <- function(g, tracts, ...){
	ids <- id_lookup(tracts, ...)

	coords <- coords_df(tracts, km = FALSE)

	if(is.null(V(g)$geoid)){
		coords <- coords %>%
			filter(geoid %in% names(V(g)))

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

#' Construct an igraph graph object from the centroids of tracts stored as an sf data.frame, where tracts are linked iff adjacent.
#' @param tracts an sf data.frame storing polygons in the geometry column
#' @param data_df a tibble containing demographic data
#' @param divergence the divergence function with which to compare demographic distributions, one of c(DKL, euc, cum_euc).
#' @export

construct_information_graph <- function(tracts, data_df, divergence){
	adj <- tracts %>%
		filter(GEOID %in% data_df$tract) %>%
		make_adjacency()

	if('t' %in% names(data_df)){
		adj <- adj %>%
			add_temporal(unique(data_df$t))
	}
	adj %>%
		undirect() %>%
		add_data(data_df) %>%
		information_distances(divergence) %>%
		make_graph() %>%
		add_n(data_df) %>%
		add_coordinates(tracts)
}

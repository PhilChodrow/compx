#' geodesics
#' @name geodesics
#' @docType package
#' @import tidyverse igraph
NULL



weighted_sum <- function(a, b, c, d){
	if(is.null(b)){
		b = NA
	}
	if(is.null(d)){
		d = NA
	}
	(a*c + b*d) / (c + d)
}

sqdist <- function(a, b){
	if(length(a) == 0){
		a = rep(NA, dim(b)[1])
	}
	na_locs           <- is.na(t(a) %*% a)
	result            <- t(a) %*% b %*% a
	result[na_locs]   <- NA
	result
}

#
# 	if(is.null(b)){
# 		b = as.matrix(NA, length(a), length(a))
# 	}else{
# 		b[is.na(b)] = 0
# 	}
# 	as.numeric((t(a) %*% b) %*% a)
# }


make_graph_deprecated <- function(grid_data, coords){

	adj <- make_adjacency(grid_data, c('x', 'y', 't')) %>%
		select(coord_key, neighbor)

	hessians <- grid_data %>%
		select(coord_key, total, local_H)

	all_coords <- grid_data %>%
		select_(.dots = coords %>% append('coord_key')) %>%
		nest(-coord_key, .key = 'coords') %>%
		left_join(hessians, by = c('coord_key' = 'coord_key'))

	adj %>%
		left_join(all_coords, by = c('coord_key' = 'coord_key')) %>%
		left_join(all_coords, by = c('neighbor' = 'coord_key'), suffix = c('_1', '_2')) %>%
		mutate(mean_H = pmap(list(local_H_1,
								  local_H_2,
								  total_1,
								  total_2),
							 weighted_sum),
			   diff   = map2(coords_1, coords_2,  ~ as.numeric(.y) - as.numeric(.x)),
			   sqdist = pmap_dbl(list(diff, mean_H), sqdist),
			   weight = sqrt(sqdist)) %>%
		select(coord_key,neighbor, weight, sqdist) %>%
		mutate(key = paste(pmax(coord_key, neighbor), pmin(coord_key ,neighbor))) %>%
		filter(!duplicated(key),
			   coord_key != neighbor) %>%
		select(-key) %>%
		na.omit() %>%
		igraph::graph_from_data_frame(directed = FALSE)
}



all_dists <- function(g){
	distances(g, mode = 'out', weight = E(g)$weight) %>%
		as.data.frame() %>%
		rownames_to_column(var = 'cell_key') %>%
		gather(key = cell_key_2, value = dist, -cell_key) %>%
		tbl_df()
}


geodesic_centroid <- function(dists){
	dists %>%
		group_by(cell_key) %>%
		summarise(max_dist = max(dist)) %>%
		filter(max_dist == min(max_dist))
}

geodesic_center_of_mass <- function(dists, grid_data){
	pops <- grid_data %>%
		mutate(cell_key = paste0(cell, '_', t)) %>%
		select(cell_key, total)

	dists %>%
		left_join(pops, by = c('cell_key_2' = 'cell_key')) %>%
		group_by(cell_key) %>%
		summarise(wmean = weighted.mean(dist, total, na.rm = T)) %>%
		filter(wmean == min(wmean, na.rm = T))
}










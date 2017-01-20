#' geometry
#' @name geometry
#' @docType package
#' @import tidyverse sp maptools
NULL

#' Compute numerical first derivatives on a grid_data data frame with a
#' user-specified difference method.
#' @param grid_data, a data frame giving counts on a grid.
#' @param difference_method, one of c('central', 'forward', 'backward')
#' @return grid_data with new columns for derivatives, one for each grid
#' coordinate
#' @export
compute_derivs <- function(grid_data, difference_methods = NULL){

	core_names <- c('x', 'y', 'coord_key', 'group', 'n')
	nonspatial_names <- names(grid_data) %>% setdiff(core_names) %>% setdiff('cell')
	coord_names <- c('x', 'y') %>% append(nonspatial_names)


	diff_methods <- list()
	for(coord in coord_names){
		diff_methods[coord] = 'central'
	}

	if(!is.null(difference_methods)){
		diff_methods <- diff_methods %>% modifyList(difference_methods)
	}

	arranged <- grid_data %>%
		arrange_(.dots = nonspatial_names %>% append(c('x', 'y', 'group'))) %>%
		group_by_(.dots = nonspatial_names %>% append('coord_key')) %>%
		mutate(total = sum(n),
			   p = n / total)

	methods <- list(
		central  = ~((lead(p) - lag(p)) / (lead(coord) - lag(coord))),
		forward  = ~((lead(p) - p) / (lead(coord) - coord)),
		backward = ~((p - lag(p)) / (coord - lag(coord)))
		)

	for(coord in coord_names){
		group_coords <- coord_names %>% setdiff(coord) %>% append('group')

		varname <- paste0('dp_d', coord)
		varval  <- lazyeval::interp(methods[[diff_methods[[coord]]]],
									coord = as.name(coord))

		arranged <- arranged %>%
			group_by_(.dots = group_coords) %>%
			mutate_(.dots = setNames(list(varval), varname))
	}
	arranged %>% ungroup()
}

#' Construct an adjacency data frame for a grid_data data frame.
#'
#' @param grid_data a grid_data data frame
#' @return a data frame giving the adjacencies between grid cells, suitable for
#' use as an edge-list to construct an igraph graph.

# make_adjacency <- function(grid_data, coords){
#
# 	adj <- grid_data[,coords %>% append('coord_key')] %>%
# 		arrange_(.dots = coords) %>%
# 		unique()
#
# 	for(coord in coords){
# 		group_coords <- coords %>%
# 			setdiff(coord)
#
# 		varname <- list(paste0('forward_', coord), paste0('backward_', coord))
# 		forward  <- lazyeval::interp(~lead(coord_key))
# 		backward <- lazyeval::interp(~lag(coord_key))
# 		varval <- list(forward, backward)
#
# 		adj <- adj %>%
# 			group_by_(.dots = group_coords) %>%
# 			mutate_(.dots = setNames(varval, varname)) %>%
# 			ungroup()
# 	}
# 	gather_cols <-names(adj)[grepl(pattern = 'ward', names(adj))]
# 	gather_cols <- gather_cols %>% append('self')
# 	adj %>%
# 		mutate(self = coord_key) %>%
# 		gather_(key_col = 'direction', value_col = 'neighbor', gather_cols = gather_cols)
# }

#' Smooth a grid_data data frame by constructing a new column of proportions
#' computed on the sum of a cell and each of its neighbors.
#' @param grid_data the grid_data data frame to smooth.
#' @return a new grid_data object with smoothed proportions
#' @export
adjacency_smoother <- function(grid_data, coords){

	adj <- make_adjacency(grid_data, coords)
	adj %>%
		select(coord_key, neighbor) %>%
		left_join(grid_data, by = c('coord_key' = 'coord_key')) %>%
		group_by(coord_key, group) %>%
		summarise(smoothed_n = sum(n)) %>%
		group_by(coord_key) %>%
		mutate(smoothed_p = smoothed_n / sum(smoothed_n))
}

#' Compute the Riemannian metric in local coordinates on a grid_data object
#' @param grid_data a grid_data object
#' @param divergence the divergence function chosen to induce a metric in
#' information space. One of c('KL', 'euc', 'cum_euc').
#' @param coord_names the names of the coordinate columns in grid_data
#' @return a new grid_data object with hessians (Riemannian metrics) computed
#' in local coordinates

compute_metric <- function(grid_data, divergence, coord_names){

	cum_euc <- function(p){
		lower <- matrix(0, length(p), length(p))
		lower[lower.tri(lower, diag = T)] <- 1
		lower %*% t(lower)
	}

	hessians <- list(KL      = function(p) return(diag(1/p)),
					 euc     = function(p) diag(length(p)),
					 cum_euc = cum_euc)

	NA_multiply <- function(D_alpha, H){
		na_locs           <- is.na(t(D_alpha) %*% D_alpha)
		H[is.infinite(H)] <- 0
		result            <- t(D_alpha) %*% H %*% D_alpha
		result[na_locs]   <- NA
		result
	}

	deriv_names <- paste0('dp_d', coord_names)
	no_nest <- coord_names %>% append(c('cell', 'coord_key', 'total'))
	to_nest <- names(grid_data) %>% setdiff(no_nest)

	grid_data %>%
		nest_(key_col = 'data', nest_cols = to_nest) %>%
		mutate(D_alpha = map(data,        ~ as.matrix(.[deriv_names])),
			   H       = map(data,        ~ hessians[[divergence]](.$smoothed_p)),
			   local_H = map2(D_alpha, H, ~ NA_multiply(.x, .y)))
}

#' Compute local Riemannian metrics from tracts and a demographics data frame.
#' The user may specify the grid resolution, the type of divergence used, and
#' the numerical difference method employed
#' @param tracts an SPDF
#' @param data a data frame consisting at least three columns:
#' - tract, the key relating the data frame to the tracts SPDF
#' - group, the group labels
#' - n a count for each group label in each tract
#' - (Optional) additional nonspatial coordinates, such as a time coordinate
#' @param resolution the spatial resolution of the grid on which to compute, in
#' km.
#' @param divergence one of c('KL', 'euc', 'cum_euc') specifying the divergence
#' used to compare distributions
#' @param difference_method one of c('central', 'forward', 'backward'),
#' specifying the computation of numerical first derivatives
#' @return grid_data data frame with the the Riemannian metric tensor in local
#' coordinates for each grid cell.
#' @export
compute_geometry <- function(tracts, data, resolution, divergence = 'KL', difference_methods = NULL){

	coord_names <- names(data) %>%
		setdiff(c('group', 'n', 'tract'))
	coord_names <- c('x', 'y') %>% append(coord_names)

	grid_info <- grid_aggregate(tracts,resolution, data)
	smoothed  <- grid_info$grid_data %>% adjacency_smoother(coord_names)
	grid_data <- grid_info$grid_data %>%
		compute_derivs(difference_methods) %>%
		left_join(smoothed) %>%
		compute_metric(divergence, coord_names)
	grid_data
}


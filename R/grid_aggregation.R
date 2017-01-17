
#' grid
#' @name grid
#' @docType package
#' @import tidyverse sp maptools rgeos
NULL

#' Construct a square grid over input tracts with specified
#' spatial resolution.
#' @param tracts a SpatialPolygonsDataFrame
#' @param resolution the length in km of the grid cells
#' @return a list containing the grid as a SpatialPolygonsDataFrame, the
#' resolution, and the grid dimensions.

make_grid <- function(tracts, resolution){
	bbx <- bbox(tracts)
	length_scale <- 1/sqrt(85 * 111) * resolution
	dims <- (c(bbx[1,2] - bbx[1,1], bbx[2,2] - bbx[2,1]) / length_scale) %>%
		as.integer()
	dims <- dims + 1
	grid <- GridTopology(c(bbx[1,1], bbx[2,1]), c(length_scale,length_scale), dims)
	coords <- coordinatevalues(grid)
	grid <- grid %>%
		as.SpatialPolygons.GridTopology() %>%
		as('SpatialPolygonsDataFrame')
	proj4string(grid) <- proj4string(tracts)
	list(grid = grid, resolution = resolution, dims = dims)
}

#' Construct a data frame that relates the grid cells to the original tracts.
#' There is one row for each instance of a tract overlapping a grid cell.
#' This operation can be quite computationally expensive.
#' @param grid a SpatialPolygonsDataFrame
#' @param tracts a SpatialPolygonsDataFrame
#' @return a data frame relating grid to tracts.

make_grid_tract <- function(grid, tracts){

	# get the list of tracts associated with a fixed cell.
	get_cell_window <- function(cell, tracts, grid){
		window <- tracts[grid[cell,],]
		if(nrow(window) != 0){
			poly_i     <- raster::intersect(grid[cell,], window)
			areas      <- map_dbl(poly_i@polygons, ~.@area)
			full_areas <- map_dbl(poly_i@polygons, ~.@area)
			tract_names     <- poly_i@data$GEOID
		}else{
			tract_names = c(NA)
			areas = 0
			full_areas = 1
		}
		return(data_frame(cell = cell,
						  tract = tract_names,
						  area = areas,
						  full_area = full_areas))
	}

	# subset the cells to only those that overlap at least one element of the
	# grid.
	tracts@data$group <- 'group'
	unioned <- maptools::unionSpatialPolygons(tracts, IDs = tracts@data$group)
	cells_to_compute <- grid[unioned,] %>% row.names

	other <- row.names(grid) %>%
		setdiff(cells_to_compute)

	# map get_piece to each cell, producing a data frame
	grid_tract <- cells_to_compute %>%
		purrr::map(get_cell_window, tracts, grid) %>%
		reduce(rbind)

	# attach the other cells with NA data.
	if(length(other) > 0){
		other <- other %>%
			data.frame(cell = ., tract = NA, area = 0, full_area = 1) %>%
			tbl_df()
		grid_tract <- grid_tract %>%
			rbind(other)
	}

	# create a weight column reflecting the proportion of area contained in the
	# cell.
	grid_tract <- grid_tract %>%
		mutate(weight = area / full_area) %>%
		tbl_df() %>%
		mutate(tract = as.character(tract))

	return(grid_tract)
}

#' Construct a data frame of coordinates for a grid SPDF
#' @param grid a square grid SpatialPolygonsDataFrame
#' @param dims the dimensions of grid
#' @param resolution the spatial resolution of the grid, i.e. the units
#' corresponding to a single side of a cell.
#' @param nonspatial_coords a list of nonspatial coordinates. The standard
#' use case is the addition of time variables.
#' @return a data frame consisting of a rectangular coordinates for a grid
#' SPDF.

make_coordinates <- function(grid,  dims, resolution, nonspatial_coords = NULL){
	nonspatial_coords %>%
		append(list(cell = row.names(grid)))%>%
		expand.grid() %>%
		mutate(i = as.integer(stringr::str_sub(cell, 2)),
			   x = ((i-1) %% dims[1])*resolution,
			   y = - (as.integer((i -1) / dims[1])) * resolution,
			   cell = as.character(cell)) %>%
		select(-i)
}


#' Construct a data frame consisting of data originally keyed to a tracts
#'  object, re-keyed to the grid. The grid and tracts must be related by
#'  grid_tract.
#' @param grid_tract a data frame, produced by `make_grid_tract()`, relating
#' a grid SPDF to a tracts SPDF
#' @param data the data frame keyed to the tract object.
#' @param coords the coordinates of the grid, as produced by
#' `make_coordinates()`
#' @return a data frame giving counts from data, re-keyed to grid, as well as
#' spatial (rectangular) coordinates of the grid variables.
make_grid_data <- function(grid_tract, data, coords){
	nonspatial_names <- setdiff(names(coords), c('x','y'))
	nonspatials <- coords[nonspatial_names]

	name_order <- c('x', 'y') %>%
		append(nonspatial_names) %>%
		append(c('group', 'n'))

	grid_tract %>%
		left_join(nonspatials) %>%
		left_join(data) %>%
		mutate(n = n * weight) %>%
		group_by_(.dots = append(nonspatial_names, 'group')) %>%
		summarise(n = sum(n * weight, na.rm = T)) %>%
		ungroup() %>%
		complete(t, cell, group, fill = list(n = 0)) %>%
		left_join(coords[c('cell', 'x', 'y', 't')]) %>%
		select_(.dots = name_order)

}



#' Construct grid_data from tracts and data, given a fixed grid resolution. This
#' is a wrapper for make_grid, make_grid_tract, make_coordinates, and
#' make_grid_data.
#'
#' @param tracts an SPDF
#' @param resolution the resolution of the grid at which to aggregate, in km
#' @param data a data frame consisting at least three columns:
#' - tract, the key relating the data frame to the tracts SPDF
#' - group, the group labels
#' - n a count for each group label in each tract
#' - (Optional) additional nonspatial coordinates, such as a time coordinate
#' @return a list containing grid_tract and grid_Data.
#' @export
grid_aggregate <- function(tracts, resolution, data){
	grid_info <- tracts %>% make_grid(resolution)
	grid_tract <- make_grid_tract(grid_info$grid, tracts)

	nonspatial_names <- names(data) %>%
		setdiff(c('group', 'n', 'tract'))

	if(length(nonspatial_names) != 0){
		nonspatial_coords <- data[nonspatial_names] %>%
			as.list() %>%
			map(unique)
	}else{
		nonspatial_coords = NULL
	}

	coords <- make_coordinates(grid              = grid_info$grid,
							   dims              = grid_info$dims,
							   resolution        = grid_info$resolution,
							   nonspatial_coords = nonspatial_coords)

	grid_data <- make_grid_data(grid_tract, data, coords)


	coord_names <- c('x', 'y') %>% append(names(coords) %>% setdiff(c('cell', 'x', 'y')))


	key_val <- paste0('paste(',
					  paste(coord_names, collapse = ','),
					  ", sep = '_')")

	grid_data <- grid_data %>%
		arrange_(.dots = coord_names) %>%
		mutate_(.dots = setNames(key_val, 'coord_key'))

	return(list(grid_tract = grid_tract, grid_data = grid_data))
}

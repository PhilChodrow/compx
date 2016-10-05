#' info
#' @name info
#' @docType package
#' @import dplyr tidyr

NULL

#' @export

make_grid <- function(tracts, resolution){
	radius = 1/sqrt(85 * 111) * resolution # (roughly resolution km after lat-lon conversion)
	xx = sp::spsample(tracts, type="hexagonal", cellsize=radius)
	xxpl = sp::HexPoints2SpatialPolygons(xx)
	xxpl <- as(xxpl, 'SpatialPolygonsDataFrame')
	print(paste0(nrow(tracts@data), ' tracts || ', length(xx), ' grid cells'))
	cell_area <- rgeos::gArea(xxpl[1])

	d <- data_frame(cell = character(), tract = numeric(), area = numeric())

	for(x in row.names(xxpl@data)){
		window     <- tracts[xxpl[x,],]
		poly_i     <- raster::intersect(xxpl[x,], window)
		areas      <- sapply(poly_i@polygons, function(y) y@area)
		full_areas <- sapply(window@polygons, function(y) y@area)
		d <- rbind(d, data_frame(cell = x, tract = window@data$GEOID, area = areas, full_area = full_areas))
	}

	names(d) <- c('cell', 'tract', 'area', 'full_area')
	d <- d %>%
		mutate(weight = area / full_area)
	xxpl <- as(xxpl, 'SpatialPolygonsDataFrame')
	xxpl@data$id = row.names(xxpl@data)

	list(grid_tract = d, grid_polys = xxpl)
}

#' @export
info_analysis <- function(tracts, columns, resolution = NULL, grid_tract = NULL, grid_polys = NULL){
	if(is.null(grid_tract)){
		grid       <- make_grid(tracts, resolution)
		grid_tract <- grid$grid_tract
		grid_polys <- grid$grid_polys
	}

	data <- tracts@data[,c(columns, 'GEOID')]

	cells <- grid_tract %>%
		mutate(tract = as.character(tract)) %>%
		left_join(data, by = c('tract' = 'GEOID')) %>%
		mutate_at(columns, function(y) y * .$weight) %>%
		select_(.dots = c('cell', columns)) %>%
		data.table() %>%
		setkey('cell')

	df <- cells[, list(info = mutual_info(.SD),
					   pop = sum(.SD)),
				by = .(cell),
				.SDcols = columns]
	cell_data <- cells[,lapply(.SD, sum, na.rm = T), by = .(cell)]
	df <- df %>% left_join(cell_data, by = 'cell')

	join_cols <- setdiff(names(df), names(grid_polys))

	grid_polys@data <- grid_polys@data %>% left_join(df[,join_cols], by = c('id' = 'cell'))

	return(list(H_Y  = tracts@data[,columns] %>% colSums() %>% simplex_normalize() %>% H,
				I_XY = tracts@data[,columns] %>% mutual_info(),
				J_XY = weighted.mean(df$info, df$pop),
				grid_tract = grid_tract,
				grid       = grid_polys))
}

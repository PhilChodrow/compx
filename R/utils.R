#' utils
#'
#' @name utils
#' @docType package
#'
#' @import dplyr acs readr tidyr stringr

NULL

# convert to matrix (main input for the algorithm we are going to be working with)

#' Make a set of joined data based on specified census information at the block level.
#'
#' @param state the state for which to obtain data.
#' @param counties a vector of ints reflecting counties to use.
#' @param table_name character name of the census table t use
#' @param end_year the final year of the data set to use
#' @param span the number of years the data should extend back
#'
#' @return merged_df a data frame including both geographical and data fields.
#' @importFrom stringr str_pad
#' @export

get_census_data <- function(state, table_nums = c("B15003", 'B03002', 'C24010'), lookup_table = NULL, end_year=2014, span=5){
	library(acs, warn.conflicts = FALSE)

	counties <- acs.fetch(geography=geo.make(state=state, county="*"),
						  endyear = 2014,
						  table.number="B01003")
	counties <- as.numeric(geography(counties)[[3]])

	polys <- tigris::block_groups(state = state, county = counties, cb=TRUE)

	for(table_num in table_nums){
		data <- acs::acs.fetch(endyear = end_year,
							   span = span,
							   geography = acs::geo.make(state = state,
							   						  county = counties,
							   						  tract = '*',
							   						  block.group = '*'),
							   table.number = table_num,
							   col.names = "pretty")

		df <- data.frame(data@estimate) %>%
			tbl_df
		colnames(df) <- colnames(data@estimate)

		df <- cbind(df, data@geography) %>%
			tbl_df

		df <- df %>%
			gather(key = column, value = n, -(NAME:blockgroup)) %>%
			left_join(lookup, by = c('column' = 'old')) %>%
			filter(new != 'DELETE') %>%
			group_by(NAME, new, state, county, tract, blockgroup) %>%
			summarise(n = sum(n)) %>%
			spread(new, n, fill = 0) %>%
			mutate(GEOID = paste0(str_pad(state, 2, 'left', pad = '0'),
								  str_pad(county, 3, 'left', pad = '0'),
								  str_pad(tract, 6, 'left', pad = '0'),
								  str_pad(blockgroup, 1, 'left', pad = '0'))) %>%
			ungroup() %>%
			select(-NAME, -state, -county, -tract, -blockgroup)

		polys@data <- left_join(polys@data, df)
	}
	return(polys)
}


#' @export
checkerboard_illustration <- function(){
	expand.grid(x = 1:8, y = 1:8) %>%
		mutate(`(a)` = 1,
			   `(b)` = .5,
			   `(c)` = (x > 4)*1,
			   `(d)` = (x + y) %% 2) %>%
		gather(key = model, value = p, -x, -y) %>%
		ggplot(aes(x = x, y = y)) +
		theme_minimal() +
		theme(axis.ticks = element_blank(),
			  axis.text.x = element_blank(),
			  axis.text.y = element_blank(),
			  panel.background = element_rect(),
			  panel.grid.major = element_line(size = 0),
			  panel.grid.minor = element_line(size = 0),
			  plot.margin=unit(c(0,0,0,0),"mm")) +
		geom_tile(aes(fill = p)) +
		scale_fill_continuous(low = 'white', high = 'black ', limits=c(0,1)) +
		facet_grid(~model) +
		xlab('') +
		ylab('') +
		guides(fill=FALSE)
}


#' @export
method_illustration <- function(city, radius_km = 1){

	races <- c('Black', 'Hispanic', 'Asian', 'White', 'Other')
	tracts <- readOGR(dsn = paste0('cities/',city), layer = 'geo', verbose = FALSE)
	tracts <- tracts[tracts@data$total > 0,]

	radius = 1/sqrt(85 * 111) * radius_km # (roughly 1 km after lat-lon conversion)

	xx = spsample(tracts, type="hexagonal", cellsize=radius)
	xxpl = HexPoints2SpatialPolygons(xx)

	# Define information measures on the grid
	h <- function(i){
		window <- tracts[xxpl[i,],]@data[,c(races, 'total', 'area')]
		window <- window / window$area # total becomes density
		c(mean(window$total), 4 * mutual_info(window[,races])/radius_km^2) # returns estimated density and mutual info
	}

	# Compute measures on the grid and collect as df
	df <- 1:length(xxpl) %>%
		as.matrix() %>%
		apply(MARGIN = 1, FUN = h) %>%
		t %>%
		as.data.frame()
	names(df) <- c('density', 'info')
	df$id <- paste0('ID', row.names(df))

	hexgrid <- fortify(xxpl) %>%
		left_join(df, by = c('id' = 'id'))

	bbox <- c(min(hexgrid$long), min(hexgrid$lat), max(hexgrid$long), max(hexgrid$lat))
	map <- get_map(location = bbox, maptype = 'terrain-background',)

	p <- ggmap(map,darken = .5)
	p <- p +
		geom_polygon(data = hexgrid, aes(x = long, y = lat, group = group, fill = info), alpha = .8) +
		scale_fill_continuous(low = 'white', high = 'steelblue', name= expression(J[Y](x)), limits = c(0, log(length(races)) / (radius_km^2))) +
		theme(axis.ticks = element_blank(),
			  axis.text.x = element_blank(),
			  axis.text.y = element_blank()) +
		xlab('') +
		ylab('') +
		theme(legend.justification=c(.8,0),
			  legend.position=c(1,0),
			  legend.background = element_rect(fill = alpha('blue', 0)),
			  legend.title = element_text(colour="white"),
			  legend.text = element_text(colour="white"),
			  legend.key.size = unit(2, "mm"),
			  legend.text = element_text(size = rel(.1))) +
		ggtitle(paste0(city, ': Mean = ', round(weighted.mean(df$info, df$density), 2)))

}


#' @export
city_list <- function(file){
	cities <- read_csv(file)

	cities <- cities %>%
		separate(name, c('name','state'), sep = ',') %>%
		mutate(state = substr(county_name, start  = nchar(county_name) - 1, stop  = nchar(county_name)),
			   county_num = as.integer(substr(county_num, start = 3, stop = 5))) %>%
		select(-county_name)

	cities <- cities %>% aggregate(county_num ~ name + state + code, c, data = .) %>%
		tbl_df %>%
		left_join(cities[,c('name', 'state')]) %>%
		filter(!duplicated(x = .[,c('code', 'state')])) %>%
		mutate(state = as.character(state),
			   short_name = word(name, sep = '-'))
	cities
}

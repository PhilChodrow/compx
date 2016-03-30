#' utils
#'
#' @name utils
#' @docType package
#'
#' @import dplyr

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

get_census_data <- function(states, counties, table_num, end_year=2012, span=5){
	library(acs, warn.conflicts = FALSE)
		data <- acs::acs.fetch(endyear = end_year,
						   span = span,
						   geography = acs::geo.make(state = states, county = counties, tract = '*', block.group = '*'),
						   table.number = table_num,
						   col.names = "pretty")

	polys <- tigris::block_groups(state = states, county = counties, cb=TRUE)

	df <- cbind(data.frame(data@geography),
				data.frame(data@estimate)) %>%
		tbl_df() %>%
		mutate(GEOID = paste0(str_pad(state, 2, 'left', pad = '0'),
									 str_pad(county, 3, 'left', pad = '0'),
									 str_pad(tract, 6, 'left', pad = '0'),
									 str_pad(blockgroup, 1, 'left', pad = '0')))

	df_merged <- tigris::geo_join(polys, df, "GEOID", "GEOID")

	return(df_merged)
}

#' Get the centroids for a bunch of spatial data.
#' @param data the SpatialPolygon data for which to compute centroids
#' @return a matrix of centroid coordinates
#' @export

get_centroids <- function(data){
	centroids <- rgeos::gCentroid(data,byid=TRUE)
	centroids <- centroids@coords
	colnames(centroids) <- c('lon', 'lat')
	return(centroids)
}

#' Get matching df for a bunch of spatial data
#' @param data the SpatialPolygon data for which to compute centroids
#' @return df the complete set of data corresponding to the SpatialPolygons, keyed by ID to match with them.
#' @export

get_table <- function(data){
	exclude <- c('GEOID', 'AWATER', 'COUNTYFP', 'TRACTCE', 'BLKGRPCE', 'AFFGEOID', 'NAME', 'LSAD', 'ALAND','NAME.1', 'state', 'county', 'tract', 'blockgroup', 'GEOID.1', 'STATEFP')
	tab <- data@data
	tab <- tab[,(!names(tab) %in% exclude)]
	tab <- tab[,!(grepl('Total', names(tab)))]
	return(tab)
}

#' @export
format_data <- function(acs_data){
	P <- acs_data %>%
		get_table() %>%
		as.matrix()
	P <- P + 1 # FUDGE FOR DEALING WITH ZERO ENTRIES
	P <- P %>%
		t %>%
		scale(., center=FALSE, scale = colSums(.)) %>%
		t %>%
		na.omit
	X <- get_centroids(acs_data)
	X <- X[row.names(P),]
	data <- list(X = X, P = P)
}

#' @export
get_dims <- function(data, K){
	n <- ncol(data$X)
	I <- nrow(data$P)
	J <- ncol(data$P)
	list(n = n, I = I, J = J, K = K)
}




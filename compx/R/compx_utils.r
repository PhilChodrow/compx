#' compx
#'
#' @name compx_utils
#' @docType package
NULL


# convert to matrix (main input for the algorithm we are going to be working with)

library(tigris)
library(acs)
library(stringr)
library(dplyr)
library(leaflet)
library(rgeos)

#' Make a set of joined data based on specified census information at the block level.
#'
#' @param state the state for which to obtain data.
#' @param counties a vector of ints reflecting counties to use.
#' @param table_name character name of the census table t use
#' @param end_year the final year of the data set to use
#' @param span the number of years the data should extend back
#'
#' @return merged_df a data frame including both geographical and data fields.
#'
#' @export

get_census_data <- function(state, counties, table_num, end_year=2012, span=5){

	data <-acs.fetch(endyear = 2012,
					   span = 5,
					   geography = geo.make(state = state, county = counties, tract = '*', block.group = '*'),
					   table.number = table_num,
					   col.names = "pretty")

	polys <- block_groups(state = state, county = counties, cb=TRUE)

	df <- cbind(data.frame(data@geography),
					   data.frame(data@estimate)) %>%
		tbl_df %>%
		mutate(GEOID = paste0(str_pad(state, 2, 'left', pad = '0'),
							  str_pad(county, 3, 'left', pad = '0'),
							  str_pad(tract, 6, 'left', pad = '0'),
							  str_pad(blockgroup, 1, 'left', pad = '0')))

	df_merged <- geo_join(polys, df, "GEOID", "GEOID")

	return(df_merged)
}

#' Get the centroids for a bunch of spatial data.
#'
#' @param data the SpatialPolygon data for which to compute centroids
#'
#' @return a matrix of centroid coordinates
#'
#' @export

get_centroids <- function(data){
	require(rgeos)
	centroids <- gCentroid(data,byid=TRUE)
	centroids <- centroids@coords
	return(centroids)
}

#' Get matching df for a bunch of spatial data
#'
#' @param data the SpatialPolygon data for which to compute centroids
#'
#' @return df the complete set of data corresponding to the SpatialPolygons, keyed by ID to match with them.
#'
#' @export

get_table <- function(data){
	data@data <- data@data %>%
		mutate(ID = row.names(data))
	return(data@data)
}


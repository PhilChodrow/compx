#' visualization
#' @name visualization
#' @docType package
#' @import dplyr ggplot2

NULL

#' Construct a function that returns the total spatial influence at a point
#' @param pars  list of model parameters
#' @export

inf_constructor <- function(pars){
	influence <- function(x1, x2){
		x <- c(x1, x2)
		densities <- normal_vec(x = x, pars = pars)
		as.numeric(pars$C %*% densities)
	}
}

#' Make a contour plot of spatial influence, overlaid on points
#' @param data the empirical data
#' @param pars the parameters
#' @export

spatial_plot <- function(data, pars, grid_density = 50){
	X <- data$X %>%
		data.frame

	lims <- X %>%
		summarise(x0 =  min(lon),
				  x1 =  max(lon),
				  y0 =  min(lat),
				  y1 =  max(lat))

	x_grid <- seq(lims$x0, lims$x1, (lims$x1 - lims$x0)/grid_density)
	y_grid <- seq(lims$y0, lims$y1, (lims$y1 - lims$y0)/grid_density)
	gg <- expand.grid(x = x_grid, y = y_grid)

	influence <- inf_constructor(pars)
	gg$z <- mapply(influence, gg$x, gg$y)

	p <- ggplot(gg, aes(x = x, y = y, z = z)) + stat_contour(aes(color = ..level..))
	p + geom_point(data = X, aes(x = lon, y = lat, z = 0))
}

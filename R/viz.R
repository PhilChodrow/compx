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
		densities <- normal_vec(x, pars$Mu, pars$Sigma)
		as.numeric(pars$C %*% densities)
	}
}

#' Make a contour plot of spatial influence, overlaid on points
#' @param data the empirical data
#' @param pars the parameters
#' @export

spatial_plot <- function(data, pars){
	X <- data$X %>% do.call(rbind, .) %>%
		data.frame

	lims <- X %>%
		summarise(x0 = 10 * min(X1),
				  x1 = 10 * max(X1),
				  y0 = 10 * min(X2),
				  y1 = 10 * max(X2))
	m <- 50 # grid density

	x_grid <- seq(lims$x0, lims$x1, (lims$x1 - lims$x0)/50)
	y_grid <- seq(lims$y0, lims$y1, (lims$y1 - lims$y0)/50)
	gg <- expand.grid(x = x_grid, y = y_grid)

	influence <- inf_constructor(pars$Mu, pars$Sigma, pars$C)
	gg$z <- mapply(influence, gg$x, gg$y)

	p <- ggplot(gg, aes(x = x, y = y, z = z)) + stat_contour(aes(color = ..level..))
	p + geom_point(data = X, aes(x = X1, y = X2, z = 0))
}

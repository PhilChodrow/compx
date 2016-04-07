library(compx)
library(alabama)
library(dplyr)

context('Optimization Functionality')

states   <- 'MA'
counties <- c(025)

data <- get_census_data(state = states, counties = counties, table_num = "B19001") %>%
	format_data()

i <- 10
data$P <- data$P[1:i,]
data$X <- data$X[1:i,]

dims <- get_dims(data, K = 5)
pars <- random_params(dims)
v <- to_vector(pars)

problem <- make_problem(data, dims, obj_fun = sum_of_squares)

sol <- auglag(par = v,
	   fn = problem$objective,
	   hin = problem$hin,
	   heq = problem$heq)

pars <- from_vector(sol$par,dims)
ests <- est(data, pars)
spatial_plot(data, pars)



# Some fixes achieved, but still seeing pathological convergence to a value regardless of K -- 3.904 in this case. Not quite sure what this number means, but I guess it's the objective value for optimally predicitng using just K = 1 representative distribution? So the algorithm currently isn't figuring out how to use its extra parameters.
# One possibility: are the spatial influence functions not matching up correctly? So the alg has to `normalize` by killing one of the two distributions?

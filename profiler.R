library(alabama)
library(tidyr)
library(compx)
library(dplyr)
library(profvis)

# Data dimensions
i <- 10
j <- 5
n <- 1
k <- 5
data_dims <- list(I = i, J = j, K = k, n = n)
rand_data <- random_1d_data(data_dims)


K = 2
model_dims <- get_dims(rand_data, K = K)
problem <- make_problem(rand_data, model_dims)

pars <- warm_start(rand_data, model_dims)
vec <- to_vec(pars)


p <- profvis({
opt_pars <- auglag(par = vec,
				   fn = problem$objective,
				   heq = problem$heq,
				   hin = problem$hin,
				   gr = problem$gradient,
				   control.optim = list(reltol = .0001),
				   control.outer=list(kkt2.check=FALSE,
				   				   eps = .001,
				   				   trace = TRUE))
})

p

# Todo: analytically compute the derivative of the squaring function
# Optimize grad_term, which should be able to run much more efficiently -- no reason for it to take this much longer to evaluate the gradient than the objective.



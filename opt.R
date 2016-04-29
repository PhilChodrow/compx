library(dplyr)
library(compx)
library(alabama)

i <- 15
j <- 2
n <- 1
k <- 10

q <- runif(k*j,0,1) %>%
	matrix(k,j) %>%
	apply(MARGIN = 1, FUN = simplex_normalize) %>%
	t

mu <- 1:k %>% matrix
sig <- .5

b <- runif(k, 0, 1)

v <- cbind(mu, sig) %>% t %>% c

x <- seq(0,k+1,by = .1) %>% matrix

dims <- list(I = i, J = j, K = k, n = n)
pars <- list(Q = q, b = b, V = v)

gen_vec <- to_vec(pars)

p <- x %>% apply(MARGIN = 1, FUN = psi, vec = gen_vec, dims = dims) %>% t

heatmap(p, Rowv = NA, Colv = NA, scale = 'column')

data <- list(X = x, P = p)


K = 2
model_dims <- get_dims(data, K = K)
problem <- make_problem(data, model_dims)

pars <- warm_start(data, model_dims)
vec <- to_vec(pars)



opt_pars <- auglag(par = vec,
	   fn = problem$objective,
	   heq = problem$heq,
	   hin = problem$hin
	   # ,
	   # gr = problem$gradient
	   # ,
	   # control.outer = list(sig0 = 100,  i.scale = 1/100) # not sure if this is right, but we REALY don't want to violate the inequality constraints!
	   # ,
	   # hin.jac = problem$grad_hin
	   # ,
	   # heq.jac = problem$grad_heq
	   )


par(mfrow = c(2,2))

Psi(X = data$X, vec = opt_pars$par, dims = dims) %>%
	t %>%
	heatmap(Rowv = NA, Colv = NA, scale = 'none')

 heatmap(p, Rowv = NA, Colv = NA, scale = 'none')

 dat <- p %>% data.frame %>% tbl_df
 mod <- Psi(X = data$X, vec = opt_pars$par, dims = dims) %>%
 	t %>% data.frame %>% tbl_df



library(ggplot2)


problem$objective(opt_pars$par) # is this just being calculated wrong?
# sure looks like it!
# WORKING HYPOTHESIS: the objective function is wrong, possibly on the DKL side,
# in that it's not capturing approximations to the data, but rather doing something else.
# we can see this by the fact that it has very low values, essentially zero, for a
# collection of distributions that are not even close to the same. So, we need
# to change the way that the objective is calculated and ensure that it is zero on
# identical arrays, but only on identical arrays.

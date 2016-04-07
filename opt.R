library(compx)
library(alabama)


states   <- 'MA'
counties <- c(025)

data <- get_census_data(state = states, counties = counties, table_num = "B19001") %>%
	format_data()

i <- 10
data$P <- data$P[1:i,]
data$X <- data$X[1:i,]

dims <- get_dims(data, K = 2)
n <- dims$n
K <- dims$K
J <- dims$J
I <- dims$I

center <- data$X %>% colMeans()

x_range <- max(data$X[,'lon']) - min(data$X[,'lon'])
y_range <- max(data$X[,'lat']) - min(data$X[,'lat'])

# initialize spatial responsibility functions
f <- function(i){
	mu_x <- runif(1, center['lon'] - x_range, center['lon'] + x_range)
	mu_y <- runif(1, center['lat'] - x_range, center['lat'] + x_range)
	sigma <- c(sqrt(x_range), 0, sqrt(y_range))
	v <- c(mu_x, mu_y, sigma)
}

V <- 1:K %>%
	matrix %>%
	apply(1, f) %>%
	c

# V <- rnorm(n*(n+3) / 2 * K)


b <- runif(K, 0, 1)
Q <- random_representatives(dims)
par <- list(b = b, V = V, Q = Q)
vec <- to_vec(par)

par <- from_vec(vec, dims)


n_params <- K * (n*(n+3)/2 + 1 + J)

problem <- make_problem(data, dims)

obj_term(data$X[1,], data$P[1,], vec)
psi(x = data$X[1,], vec = vec, dims = dims)

problem$objective(vec)

problem$gradient(vec)


opt_pars <- auglag(par = vec,
	   fn = problem$objective,
	   heq = problem$heq,
	   hin = problem$hin
	   # ,
	   # gr = problem$gradient
	   # ,
	   control.outer = list(sig0 = 100,  i.scale = 1/100) # not sure if this is right, but we REALY don't want to violate the inequality constraints!
	   # ,
	   # hin.jac = problem$grad_hin
	   # ,
	   # heq.jac = problem$grad_heq
	   )



opt_pars$par %>% from_vec(dims)

opt_pars$gradient

problem$hin(opt_pars$par) > 0


problem$gradient(opt_pars$par) - opt_pars$gradient

jacobian(problem$objective, opt_pars$par) %>% c

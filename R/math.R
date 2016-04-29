#' math
#' @name math
#' @docType package
#' @import dplyr mvtnorm matrixcalc numDeriv magic mvnfast gdata

NULL

# -----------------------------------------------------------------------------
# INFORMATION GEOMETRY
# -----------------------------------------------------------------------------

#' Check whether a vector is an element of the probability simplex
#' @param p vector the vector to check
#' @param allow_zero boolean whether to allow elements with zero entries (boundary of simplex)
#' @return boolean whether p is a valid probability distribution
#' @export

simplex_check <- function(p, allow_zero = TRUE){
	if(min(is.na(p)) == 1) return(FALSE)
	nonneg <- ifelse(allow_zero, min(p >= 0), min(p > 0))
	normed <- abs(sum(p) - 1) < .1 # numerical tolerance
	nonneg & normed
}

#' Normalize a nonnegative vector so that it lies in the probability simplex
#' @param p vector the vector
#' @export

simplex_normalize <- function(p){
	if (min(p >= 0) == 0){
		stop("Vector has negative entries")
	}
	normed <- p / sum(p)
	if (min(p > 0) == 0){
		warning('p lies on simplex boundary')
	}
	normed
}

#' Find the Kullback-Leibler divergence of two empirical distributions.
#' @param p the 'true' distribution
#' @param q the estimated distribution
#' @return numeric the KL divergence between p and q
#' @export

DKL <- function(p,q){
	if(!is.numeric(p) | !is.numeric(q)){
		return(NaN)
	}

	if(!simplex_check(p) | !simplex_check(q)){
		message <- paste0('p or q are not on the simplex: sum(p) = ', sum(p), ' and sum(q) = ', sum(q))
		warning(message)
	}
	if(length(p) != length(q)){
		stop('Distribution alphabets are different size')
	}

	drop <- p < 10^(-10)
	p <- p[!drop]
	q <- q[!drop]

	return(as.numeric(p %*% log(p/q)))
}

#' Find the entropy of a distribution
#' @export

H <- function(p){
	if(!simplex_check(p)){
		message <- paste0('p is not on the simplex: sum(p) = ', sum(p))
		warning(message)
	}
	drop <- p < 10^(-10)
	p <- p[!drop]
	as.numeric(- p %*% log(p))
}


# -----------------------------------------------------------------------------
# SPATIAL RESPONSIBILITY FUNCTIONS
# -----------------------------------------------------------------------------


#' Compute the gradient of the squaring function for matrices.
#' @param sig the input vector
#' @return the jacobian
#' @export
square_grad <- function(sig){
	S <- sig %>% UT_ravel()
	# for each element in the lower triangle
	f <- function(k){
		M <- ((1:length(sig) == k)*1) %>% UT_ravel()
		(M %*% t(S) + S %*% t(M)) %>% UT_unravel()
	}
	1:length(sig) %>%
		matrix %>%
		apply(MARGIN = 1, FUN = f)
}



#' Compute the value of a normal distribution at point x with parameters v
#' v must have length n*(n+3)/2, where n is the length of x
#' @param x the location at which to evaluate the normal density
#' @param v the vector of parameters, in format c(mu, sigma)
#' @return numeric, the value of the normal density at x
#' @export

lambda <- function(x, v){
	n <- length(x)
	l <- n*(n+3)/2
	# assert_that(length(v) == l)
	pars <- v2p(v)
	dmvn(X = x, mu = pars$mu, sigma = pars$sigma %*% pars$sigma)
}

#' Compute the value of multiple normal distributions at point x with parameters V
#' V must have length divisible by n*(n+3)/2, where n is the length of x
#' The number of normal densities to evaluate is inferred from the length of x and
#' the length of V.
#' @param x the location at which to evaluate the normal density
#' @param V the vector of parameters, in format c(mu1, sigma1,...,muK,sigmaK)
#' @return numeric vector, the value of the K normal densities at x
#' @export

Lambda <- function(x,V){
	n <- length(x)
	l <- n*(n+3)/2
	K <- length(V)/l

	# assert_that(K %% 1 == 0)

	M <- t(matrix(V, l, K))
	apply(M, MARGIN = 1, lambda, x = x)

}

#' Compute the gradient of the normal distribution at point x with respect to its parameters
#' v. v must have length n*(n+3)/2, where n is the length of x.
#' @param x the point at which to evaluate the gradient.
#' @param v the parameters at which to evaluate the gradient
#' @return numeric vector, the gradient at x and v.
#' @export

d_lambda <- function(x, v){
	n <- length(x)
	l <- n*(n+3)/2
	assert_that(length(v) == l)

	par <- v2p(v)

	grads <- mvnorm.grad(x = x,
						 mu = par$mu,
						 sigma = par$sigma %*% par$sigma)

	sig_base <- grads$sigma.grad %>%
		UT_unravel()

	sig_grad <- sig_base %*% square_grad(UT_unravel(par$sigma))
	c(grads$mu.grad, sig_grad) %>% matrix() %>% t
}

#' Compute the gradient of a vector of normal distributions at point x with respect
#' to parameters V.
#' V must have length divisible by n*(n+3)/2, where n is the length of x.
#' The number of normal densities to evaluate is inferred from the length of x and
#' the length of V.
#' @param x the point at which to evaluate the gradient.
#' @param V the parameters at which to evaluate the gradient
#' @return numeric matrix, the gradient at x and V.
#' @export

d_Lambda <- function(x, V){

	n <- length(x)
	l <- n*(n+3)/2
	K <- length(V)/l

	assert_that(K %% 1 == 0)

	split(V, ceiling(seq_along(V)/l)) %>%
		lapply(d_lambda, x = x) %>%
		do.call(adiag, .)
}

#' Compute the gradient of a normal distribution with respect to its parameters.
#' @param x a vector denoting location at which gradient is computed
#' @param mua vector denoting mean of multivariate normal
#' @param sigma covariance matrix
#' @param log logical variable indicating whether density or log-density should be used
#' @return list of vector and matrix, the gradient with respect to the mean and covariance matrix.
#' @author Ravi Varadhan, Johns Hopkins University, June 24, 2009

mvnorm.grad <- function(x, mu, sigma, log=FALSE) {

	siginvmu <- c(solve(sigma, x-mu))
	mu.grad <- siginvmu
	sigma.grad <- tcrossprod(siginvmu) - solve(sigma)  # is there a way to avoid inverting `sigma'?
	diag(sigma.grad) <- 0.5 * diag(sigma.grad)  # this is required to account for symmetry of covariance matrix

	if (log) list(mu.grad=mu.grad, sigma.grad=sigma.grad) else
	{
		f <- dmvn(X = x, mu=mu, sigma=sigma)
		list(mu.grad=mu.grad*f, sigma.grad=sigma.grad*f)
	}
}

#' Compute the normalized Hadamard product of two vectors.
#' @param x first vector
#' @param y second vector
#' @return numeric
#' @export

phi <- function(x,y) x * y / x %*% y

#' Compute the jacobian of phi
#' @param x first vector
#' @param y second vector
#' @return numeric matrix
#' @export

d_phi <- function(x,y){
	dot <- as.numeric(x %*% y)
	dx <- t(dot * diag(y) - y %*% t(x * y)) / dot^2
	dy <- t(dot * diag(x) - x %*% t(x * y)) / dot^2
	cbind(dx, dy)
}

# -----------------------------------------------------------------------------
# OBJECTIVE FUNCTIONS
# -----------------------------------------------------------------------------

#' Compute psi, the estimated distributon at a point x with parameters vec
#' @param x
#' @param vec
#' @param dims
#' @return matrix, the estimate
#' @export
psi <- function(x, vec, dims){
	par <- from_vec(vec, dims)
	phi(par$b, Lambda(x, par$V)) %*% par$Q
}

#' @export
Psi <- function(X, vec, dims){
		t(apply(X, MARGIN = 1, FUN = psi, vec = vec, dims = dims))
}

#' Compute d_psi, the gradient of psi at a point x with respect to its parameters vec
#' @param x
#' @param vec
#' @param dims
#' @return matrix, the gradient
#' @export

d_psi <- function(x, vec, dims){
	par <- from_vec(vec, dims)
	b <- par$b
	V <- par$V
	Q <- par$Q
	K <- dims$K
	J <- dims$J
	# derivative with respect to b
	db <- t(d_phi(b, Lambda(x,V))[,1:K]) %*% Q %>% t

	# derivative with respect to V
	dV <- t(d_phi(b, Lambda(x,V))[,(K+1):(2*K)] %*% d_Lambda(x,V)) %*% Q %>% t

	# derivative with respect to Q
	dQ <- phi(par$b, Lambda(x,par$V)) %>%
		rep(J) %>%
		split(., ceiling(seq_along(.)/K)) %>%
		lapply(function(v) t(matrix(v))) %>%
		do.call(adiag, .)

	cbind(db, dV, dQ)
}


#' Construct the objective function for a problem given data and problem dimensions
#' @param data
#' @param dims
#' @export

obj_constructor <- function(data, dims){
	# refactor this so that we compute the matrix of estimates first, and then loop
	# over that and P.

	obj <- function(vec){
		estimates <- Psi(data$X, vec, dims)

		1:dim(data$P)[1] %>%
			matrix %>%
			apply(1, function(i) DKL(data$P[i,], estimates[i,])) %>%
			sum
	}

	obj
}
#' Construct the gradient of the objective function for a problem given data and problem dimensions
#' @param data
#' @param dims
#' @export

grad_constructor <- function(data, dims){
	grad_term <- function(x, p, vec, dims){
		first <- p / psi(x, vec, dims)
		second <- d_psi(x, vec, dims)
		- first %*% second
	}
	obj_grad <- function(vec){
		X <- data$X
		P <- data$P
		1:dim(P)[1] %>%
			matrix %>%
			apply(1, function(i) grad_term(X[i,], P[i,], vec, dims)) %>%
			rowSums()
	}
}





#' gradients
#' @name gradients
#' @docType package
#' @import dplyr mvtnorm matrixcalc numDeriv

NULL

#' @export
grad_obj <- function(P, Q, obj_fun = 'DKL'){
	DKL_grad <- function(p,q){
		p/q
	}
	sum_squares_grad <- function(p,q){
		2*(p - q)
	}

	grads <- list(DKL = DKL_grad, sum_of_squares = sum_squares_grad)
	chosen_grad <- grads[[obj_fun]]
	1:dim(P)[2] %>%
		matrix %>%
		apply(1, function(i) chosen_grad(P[,i], Q[,i])) %>%
		colSums
}

#' cheating to use the jacobian function, may need to work on that
#' @export
square_grad <- function(sig){
	f <- function(sig){
		X <- sig %>% UT_ravel()
		X %*% X %>% UT_unravel() # X is itself symmetric
	}
	jacobian(func = f, x = sig)
}

#' @export
d_lambda <- function(par, n, x){
	mu <- par[1:n]
	sig <- par[(n+1):length(par)]

	sig_mat <- sig %>% UT_ravel()

	Sigma <- t(sig_mat) %*% sig_mat

	grads <- mvnorm.grad(x = x,
					 mu = mu,
					 sigma = Sigma)

	sig_base <- grads$sigma.grad %>%
		UT_unravel()

	sig_grad <- sig_base %*% square_grad(sig)
	c(grads$mu.grad, sig_grad) %>% matrix()
}


####################################################
mvnorm.grad <- function(x, mu, sigma, log=FALSE) {
	# x = a vector denoting location at which gradient is computed
	# mu = a vector denoting mean of multivariate normal
	# sigma = covariance matrix
	# log = logical variable indicating whether density or log-density should be used
	# Note:  gradient is computed with respect to location parameters, and with respect to variances and covariances
	#
	# Author:  Ravi Varadhan, Johns Hopkins University, June 24, 2009
	#
	siginvmu <- c(solve(sigma, x-mu))
	mu.grad <- siginvmu
	sigma.grad <- tcrossprod(siginvmu) - solve(sigma)  # is there a way to avoid inverting `sigma'?
	diag(sigma.grad) <- 0.5 * diag(sigma.grad)  # this is required to account for symmetry of covariance matrix

	if (log) list(mu.grad=mu.grad, sigma.grad=sigma.grad) else
	{
		f <- dmvnorm(x, mean=mu, sigma=sigma, log=FALSE)
		list(mu.grad=mu.grad*f, sigma.grad=sigma.grad*f)
	}
}


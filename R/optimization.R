#' optimization
#' @name optimization
#' @docType package
#' @import dplyr matrixcalc assertthat

NULL

#' Create the objective and constraints for the problem
#' @param
#' @export

make_problem <- function(data, dims){

	K <- dims$K
	J <- dims$J

	objective <- obj_constructor(data = data, dims = dims)
	gradient  <- grad_constructor(data = data, dims = dims)

	heq <- function(vec){
		pars <- from_vec(vec, dims)
		rowSums(pars$Q) - 1
	}

	# needs checking, but the command
	# jacobian(function(Q) rowSums(Q) - 1, Q) gives K identity matrices of size J
	grad_heq <- function(vec){
		pars <- from_vec(vec, dims)
		q_component <- rep(K,J) %>%
			matrix %>%
			lapply(function(i) diag(i)) %>%
			do.call(cbind,.)

		padding <- matrix(0, K, length(pars$b) + length(pars$V))
		cbind(padding, q_component)
	}

	#
	hin <- function(vec){
		pars <- from_vec(vec, dims)
		cbind(pars$Q, pars$b) %>% as.numeric()
	}

	grad_hin <- function(vec){
		pars <- from_vec(vec, dims)
		b_component <- pars$b %>% length %>% diag
		b_pad <- matrix(0, length(pars$b), length(pars$V) +length(pars$Q))
		Q_component <- pars$Q %>% length %>% diag
		Q_pad <- matrix(0,length(pars$Q), length(pars$b) + length(pars$V))
		rbind(cbind(b_component, b_pad), cbind(Q_pad, Q_component))
	}

	list(objective = objective,
		 gradient = gradient,
		 heq = heq,
		 grad_heq = grad_heq,
		 hin = hin,
		 grad_hin = grad_hin)
}

#' Find a warm-start, feasible solution at which to initialize the problem.
#' @param data the data
#' @param dims the dimension of the problem, but only K (the number of component distributions) is used
#' @return V a vector of candidate parameters
#' @export

warm_start <- function(data, dims){

	pars <- list()

	pars$b <- rep(1.0/dims$K, dims$K) # b is uniform to start
	mu <- kmeans(data$X, dims$K)$centers # plausible centroids

	# all covariance matrices equal at first
	lims <- apply(data$X, 2, range)
	ranges <- lims[2,] - lims[1,]
	area <- prod(ranges)
	d <- area / (2*dims$K) # find the diagonal element

	sig <- 1:dims$K %>% # construct K identical, diagonal matrices
		matrix %>%
		apply(1, function(i) UT_unravel(diag(d, dims$n))) %>%
		matrix

	# compute plausible starting values for component distributions Q.

	M <- cbind(mu, sig)
	make_component <- function(k){
		v <- M[k,]
		weights <- apply(data$X, MARGIN = 1, lambda, v = v)
		(weights %*% data$P) / sum(weights)
	}
	pars$Q <- 1:dims$K %>%
		matrix %>%
		apply(MARGIN = 1, make_component) %>%
		t
	assert_that(all.equal(rowSums(pars$Q), rep(1, dims$K)))

	pars$V <- M %>% t %>% c

	pars

}

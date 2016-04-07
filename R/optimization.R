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

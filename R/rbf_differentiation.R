make_matrices <- function(data, tracts){
	coord_vars <- c('x', 'y')
	spatial <- coords_df(tracts, km = T)
	all <- spatial %>%
		left_join(data, by = c('geoid' = 'tract')) %>%
		spread(key = group, value = n, fill = 0)
	if('t' %in% names(all)){
		coord_vars <- coord_vars %>% append(c('t'))
		all <- all %>%
			mutate(key = paste(geoid, t, sep = '_'))
	}else{
		all <- all %>%
			mutate(key = geoid)
	}
	group_vars <- unique(data$group)

	list(X = coord_vars, N = group_vars) %>%
		map(~all[,.] %>% as.matrix()) %>%
		map(`row.names<-`, all$key)
}

distance_matrix <- function(X, sigma = c(1, 1, 1)){
	1:dim(X)[2] %>%
		map(~ 1/(2*sigma[.]) * dist(X[,.])) %>%
		map(as.matrix) %>%
		reduce(`+`)
}

kernel_matrix <- function(X, sigma = c(1,1,10)){
	distances <- distance_matrix(X, sigma)
	K <- exp(-distances^2)
	return(K)
}

rbf_smoother <- function(X, N, sigma = c(1,1,10)){
	K <- kernel_matrix(X, sigma)
	unnorm <- t(t(N) %*% as.matrix(K))
	normed <- unnorm / rowSums(unnorm) * rowSums(N)
	return(list(N = normed, K = K))
}

rbf_derivatives <- function(X,N, sigma = rep(1, dim(X)[2])){

	Sigma <- diag(1/(2*sigma))
	Beta <- solve(Sigma)
	row.names(Beta) <- colnames(X)
	colnames(Beta) <- colnames(X)
	smoothed <- rbf_smoother(X, N, sigma)
	N_bar <- smoothed$N
	K <- smoothed$K

	dN_bar <- function(i){
		t(N) %*% ((sweep(X, 2, X[i,]) %*% Beta) * K[,i])
	}

	dp <- function(i){
		(dN_bar(i) * sum(N_bar[i,]) - t(N_bar[i,,drop = F]) %*% colSums(dN_bar(i))) / (sum(N_bar[i,]))^2
	}

	d_alpha <- 1:nrow(X) %>%
		map(dp)
	list(D_alpha = d_alpha, N_bar = N_bar)
}

#' Needs to handle time variables gracefully, and figure out problem
#' with DKL_

compute_metric_rbf <- function(tracts, data, sigma = c(1,1,10), hessian = euc_){

	has_temporal <- 't' %in% names(data)

	mats <- make_matrices(data, tracts)
	X <- mats$X
	N <- mats$N
	derivs <- rbf_derivatives(X, N, sigma = sigma)
	D_alpha <- derivs$D_alpha
	N_bar <- derivs$N_bar

	H <- (N_bar / rowSums(N_bar)) %>%
		split(., 1:nrow(.)) %>%
		map(hessian)

	g <- map2(D_alpha, H, NA_multiply)
	out <- data_frame(key = row.names(N), total = rowSums(N), g = g)
	if(has_temporal){
		out <- out %>%
			separate(key, into = c('geoid', 't'), sep = '_') %>%
			mutate(t = as.numeric(t))
	}else{
		out <- out %>% rename(geoid = key)
	}
	return(out)
}

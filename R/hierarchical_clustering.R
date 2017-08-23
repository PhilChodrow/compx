#' hierarchical_clustering
#' @name hierarchical_clustering
#' @docType package
#' @import tidyverse igraph sf
NULL


#' Conduct hierarchical clustering of either an sf geography with demographic data or an information graph
#' @param input either an sf data.frame of polygons or an information graph as constructed by construct_information_graph
#' @param data if input is an sf data.frame, the demographic data for each geographic unit
#' @return a, an hclust object
#' @export
h_clust <- function(input, data){

	if('sf' %in% class(input)){
		tract_lookup <- id_lookup(input)
		adj <- st_relate(input, pattern = '****1****', sparse = TRUE)

		data_lookup <- data %>%
			nest(-tract) %>%
			mutate(n = map(data, ~.$n)) %>%
			select(-data) %>%
			left_join(tract_lookup, by = c('tract' = 'geoid'))

		adj <- 1:length(adj) %>%
			map(~data_frame(from = as.character(.),
							to = as.character(adj[[.]]))) %>%
			reduce(rbind) %>%
			filter(from != to) %>%
			left_join(data_lookup, by = c('from' = 'row')) %>%
			left_join(data_lookup, by = c('to' = 'row'), suffix = c('_from', '_to')) %>%
			rename(from_n = n_from, to_n = n_to)
		n <- nrow(input)
		names <- input$cluster

		M <- data %>%
			group_by(group) %>%
			summarise(n = sum(n)) %>%
			select(n) %>%
			unlist()

	}else{
		adj <- input %>%
			as.directed(mode = 'mutual') %>%
			as_long_data_frame() %>%
			dplyr::as_data_frame()

		M <- V(input)$n %>% reduce(`+`)
		n <- length(V(input))
		names <- V(input)$name
	}

	divergence <- function(n,m){
		p <- n / sum(n)
		q <- m / sum(m)
		p_bar <- sum(n) / sum(m + n) * p + sum(m) / sum(m + n) * q
		r <- M / sum(M)
		sum(n) / sum(M) * DKL(p, r) +
			sum(m) / sum(M)*DKL(q, r) -
			sum(m + n) / sum(M) * DKL(p_bar, r)
	}

	adj <- adj %>%
		mutate(dist = map2_dbl(from_n, to_n, divergence)) %>%
		mutate(from = as.integer(from),
			   to = as.integer(to))

	merge <- matrix(nrow = n-1, ncol = 2)
	height <- numeric(length = n-1)
	cluster_stage <- 1

	for(i in 1:(n-1)){
		new_cluster_name <- n + i

		ix <- which.min(adj$dist)

		height[cluster_stage] <- min(adj$dist)

		i <- adj$from[ix]
		j <- adj$to[ix]

		if(i <= n){
			eye <- -i
		}else{
			eye <- i - n
		}

		if(j <= n){
			jay <- -j
		}else{
			jay <- j - n
		}

		merge[cluster_stage, c(1,2)] <- c(eye, jay)

		op <- adj %>%
			filter(from %in% c(i,j) | to %in% c(i,j))

		ij_n <- adj %>%
			filter(from == i, to == j) %>%
			select(from_n, to_n)

		ij_n <- ij_n$from_n[[1]] + ij_n$to_n[[1]]

		adj <- adj %>% anti_join(op, by = c('from' = 'from', 'to' = 'to'))

		op <- op %>%
			filter(!(from == i & to == j) & !(from == j & to == i)) %>%
			select(-dist)

		test <- op %>%
			mutate(new = rep(list(ij_n), nrow(.)))

		test$from_n <- ifelse(test$from %in% c(i,j), test$new, test$from_n)
		test$to_n   <- ifelse(test$to %in% c(i,j), test$new, test$to_n)

		test <- test %>%
			mutate(from = ifelse(from %in% c(i,j), new_cluster_name, from),
				   to = ifelse(to %in% c(i,j), new_cluster_name, to))

		test <- test %>%
			select(-new) %>%
			mutate(dist = map2_dbl(from_n, to_n, divergence))

		adj <- adj %>%
			rbind(test) %>%
			distinct(from, to, .keep_all = TRUE)

		cluster_stage <- cluster_stage + 1
	}

	a <- list()
	a$merges <- merge
	a$height <- height
	a$labels <- names
	a$order  <- names
	class(a) <- "hclust"
	a
}

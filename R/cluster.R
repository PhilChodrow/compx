#' cluster
#' @name cluster
#' @docType package
#' @import rgdal dplyr rgeos Matrix reshape2

NULL








#' @export
update_adj <- function(adj, i, j){
	eye = min(i,j); jay = max(i,j)
	adj[eye,] <- adj[eye,] + adj[jay,]
	adj[,eye] <- adj[,eye] + adj[,jay]
	adj <- adj[-jay,-jay]
	adj <- (adj > 0) * 1
	if(length(adj) == 1){
		adj <- as.matrix(adj)
	}
	adj
}


#' @export
info_loss <- function(df, p_Y, N, i, j){

	q_i <- as.numeric(df[i,])
	q_j <- as.numeric(df[j,])
	q_ij <- q_i + q_j

	tot_i <- sum(q_i)
	tot_j <- sum(q_j)

	p_i <- tot_i / N
	p_j <- tot_j / N
	p_ij <- p_i + p_j

	q_i <- q_i / tot_i
	q_j <- q_j / tot_j
	q_ij <- q_ij / (tot_i + tot_j)

	p_i * DKL(p = q_i, q = p_Y) +
	p_j * DKL(p = q_j, q = p_Y) -
	p_ij * DKL(p = q_ij, q = p_Y)

}

#' @export
make_dist <- function(df, columns, adj){
	p_Y <- colSums(df) %>% simplex_normalize()
	N <- sum(df)
	indices <- which(adj > 0, arr.ind = T)

	dists <- matrix(nrow = nrow(df), ncol = nrow(df))
	for(k in 1:nrow(indices)){
		if(indices[k,1] != indices[k,2]){
			dists[indices[k,1],indices[k,2]] <- info_loss(df,
														  p_Y,
														  N,
														  indices[k,1],
														  indices[k,2])
		}
	}
	dists
}

#' @export
update_df <- function(df, i, j){
	eye = min(i,j); jay = max(i,j)
	new_df <- df
	new_df[eye,] <- new_df[eye,] + new_df[jay,]
	new_df[-jay,]
}

#' @export
new_dists <- function(df, i, adj){
	p_Y <- colSums(df) %>% simplex_normalize()
	N <- sum(df)
	which(adj[i,] > 0) %>%
		mapply(FUN = function(j) info_loss(df, p_Y, N, i, j), .)
}


#' This is getting there, but the resulting tree cannot be cut or converted with as.dendrogram, though it can be plotted
#' @export
info_clust <- function(spdf, columns, spatial = TRUE){

	# initialize control structures
	df <- spdf@data[,columns] %>% tbl_df
	control <- (gRelate(spdf, byid = TRUE, pattern = '****1****')) %>%
		melt() %>%
		tbl_df()
	names(control) <- c('i', 'j', 'adjacent')
	control <- control %>% mutate(i = i+1, j = j+1)

	p_Y <- colSums(df) %>% simplex_normalize()
	N <- sum(df)

	loss <- function(i, j, adjacent){
		if(adjacent & i != j){
			info_loss(df, p_Y, N, i, j)
		}else{
			NA
		}
	}
	control$dist <- mapply(FUN = loss, control$i, control$j, control$adjacent)

	# update adjacency
	# update dists

	a <- list()
	a$merge <- matrix(0, 1, 2)
	a$merge <- a$merge[-1,]
	a$height <- numeric()
	a$order <- 1:(nrow(df))
	a$labels <- 1:(nrow(df))
	merge_list <- list()
	# main loop
	for(i in 1:nrow(df)){
		merge_list[[i]] <- -i
	}

	x <- numeric()
	N = nrow(df) - 1
	for(k in 1:N){ # alright for a first draft, but we need to make sure that this eventually gets down to 1.
		print(pryr::object_size(control))

		update <- control %>%
			filter(dist == min(dist, na.rm = T), i != j) %>%
			head(1)

		if(nrow(update) == 0){
			eye <- 1
			jay <- 2

			height <- info_loss(df = df,
					  colSums(df) / sum(colSums(df)),
					  N = sum(df),
					  i = eye,
					  j = jay)

			a$height <- c(a$height, height)
		}else{
			eye <- min(update$i, update$j)
			jay <- max(update$i, update$j)
			a$height <- c(a$height, update$dist)
		}

		a$merge <- rbind(a$merge, c(merge_list[[eye]], merge_list[[jay]]))
		merge_list[[eye]] <- k
		# merge_list <- merge_list[-jay]

		# df update
		df <- update_df(df, i, j)

		# adjacency update
		control$adjacent[control$i == eye] <- control$adjacent[control$i == eye] |
			control$adjacent[control$i == jay]

		control <- control %>% filter(i != jay,
									  j != jay,
									  i != j)

		# dists update
		to_update <- control %>% filter(control$i == eye)
		control$dist[control$i == eye] <- mapply(FUN = loss,
												 to_update$i,
												 to_update$j,
												 to_update$adjacent)
	}
	a$height <- cumsum(a$height)
	class(a) <- 'hclust'
	a

}

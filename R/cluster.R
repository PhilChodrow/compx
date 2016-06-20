#' cluster
#' @name cluster
#' @docType package
#' @import rgdal dplyr rgeos Matrix

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

	df <- spdf@data[,columns] %>% tbl_df
	if(spatial){
		adj <- (gRelate(spdf, byid = TRUE, pattern = '****1****') * 1) %>% Matrix(sparse = TRUE)
	} else{
		adj <- 1 %>% matrix(nrow(df), nrow(df))
	}

	dists <- make_dist(df, races, adj)

	a <- list()
	a$merge <- matrix(0, 1, 2)
	a$merge <- a$merge[-1,]
	a$height <- numeric()
	a$order <- 1:(nrow(df))
	a$labels <- 1:(nrow(df))
	# main loop
	merge_list <- list()
	for(i in 1:nrow(df)){
		merge_list[[i]] <- -i
	}

	x <- numeric()
	N = nrow(df) - 1
	for(k in 1:N){ # alright for a first draft, but we need to make sure that this eventually gets down to 1.
		print(pryr::object_size(adj))
		min_val <- min(dists, na.rm = T)
		if(is.infinite(min_val)){
			i <- 1
			j <- 2

			height <- info_loss(df = df,
					  colSums(df) / sum(colSums(df)),
					  N = sum(df),
					  i = i,
					  j = j)

			a$height <- c(a$height, height)
		}else{
			to_combine <- which(dists == min_val, arr.ind = T)[1,]
			i <- min(to_combine)
			j <- max(to_combine)
			a$height <- c(a$height, min_val)
		}
		a$merge <- rbind(a$merge, c(merge_list[[i]], merge_list[[j]]))
		merge_list[[i]] <- k
		merge_list <- merge_list[-j]


		adj <- update_adj(adj, i, j)
		df <- update_df(df, i, j)
		dists <- dists[-j, -j]
		if(k == N){
			dists <- as.matrix(dists)
		}
		new <- new_dists(df, i, adj)
		dists[i, adj[i,] > 0] <- new
		dists[adj[,i] > 0, i] <- new
		diag(dists) <- NA
		# x[k] <- mutual_info(df)
	}

	a$height <- cumsum(a$height)
	class(a) <- 'hclust'
	a

}

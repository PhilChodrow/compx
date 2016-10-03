#' cluster
#' @name cluster
#' @docType package
#' @import rgdal dplyr rgeos reshape2 data.table utils

NULL

#' Compute constrained, greedy information-maximizing agglomerative clustering on a contingency table.
#' @param df a data frame in which each row is a separate observation and each column a category.
#' @param constraint a matrix of constraints. Let constraint be a matrix of TRUE to do unconstrained clustering.
#' @return a an object of class hclust
#' @export

info_clust <- function(df, constraint){
	k <- 1
	N <- nrow(df)
	pb <- txtProgressBar(style = 3)

	merge_lookup <- empty_merge_lookup(nrow(df))
	all_merges <- data.frame(i = integer(),
							 j = integer(),
							 loss = numeric())

	while(nrow(df) > 1){
		merge_list <- make_merge_list(df, constraint)
		constraint <- update_constraint(constraint, merge_list)
		if(nrow(merge_list) != 0){
			df <- update_df(df, merge_list)
			setTxtProgressBar(pb, (N - nrow(df) + 1) / N)
			all_merges <- rbind(all_merges,
								data.frame(i = unlist(merge_lookup[merge_list$i]),
										   j = unlist(merge_lookup[merge_list$j]),
										   loss = merge_list$loss))
			for(m in 1:nrow(merge_list)){
				merge_lookup[[merge_list[m,1]]] <- k # this seems right?
				k <- k + 1
			}
			merge_lookup[merge_list[,2]] <- NULL
		}
	}
	# this needs some work, TBD
	all_merges <- all_merges %>% sort_merges()

	a <- empty_hclust(N)
	a$height <- all_merges$loss
	a$height <- cumsum(a$height)
	a$merge <- as.matrix(all_merges[,c('i','j')])
	class(a) <- 'hclust'
	a
}



#' Compute the information loss associated with merging two rows of a data frame.
#' @param df a data frame in which each column represents a category and each row the observed counts in that category.
#' @param i the first row to merge
#' @param i the second row to merge

info_loss <- function(df, i, j){

	p_Y <- colSums(df) %>% simplex_normalize()
	N <- sum(df)

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

	p_i  * DKL(q_i,  p_Y, drop_threshold = 10^(-10)) +
	p_j  * DKL(q_j,  p_Y, drop_threshold = 10^(-10)) -
	p_ij * DKL(q_ij, p_Y, drop_threshold = 10^(-10))

}

#' Make the list of merges for each outer iteration
#' @param df the dataframe at this stage of the algorithm.
#' @param constraint the constraint at this stage of the algorithm
#' @return merge_list a data frame with columns i and j (to be merged) and loss, the information loss associated with each.

make_merge_list <- function(df, constraint){
	# each node needs to be merged no more than once in each outer iteration.

	df <- dplyr::mutate(df, cluster = row_number())

	merge_list <- data_frame(i = integer(),
							 j = integer(),
							 loss = numeric())

	find_smallest <- function(i){
		candidates <- which(constraint[i,] == T)
		candidates <- candidates[candidates != i]
		losses <- candidates %>%
			as.matrix() %>%
			apply(MARGIN = 1, FUN = function(j) info_loss(df, i, j))

		return(c(candidates[which.min(losses)], min(losses)))
	}

	# main loop
	for(i in df$cluster){
		if(!(i %in% c(merge_list$i,merge_list$j))){
			min_jay <- find_smallest(i)[1]
			if(!(min_jay %in% c(merge_list$i,merge_list$j))){
				the_min <- find_smallest(min_jay)
				min_eye <- the_min[1]
				if(i == min_eye){
					merge_list <- rbind(merge_list,
										data.frame(i = min_eye, j = min_jay, loss = the_min[2]))
				}
			}
		}
	}
	merge_list
}

#' Update the constraint matrix using the current merge list. Note that, if no merges are possible, then the constraint matrix inserts a TRUE at positions [1,2] and [2,1]. In the geographical context, this scenario corresponds to disconnected components on the map, and the constraint matrix declares two of them to be connected arbitrarily.
#' @param constraint the current state of the constraint matrix.
#' @param merge_list the list of observations to be merged.
#' @return matrix the updated constraint matrix

update_constraint <- function(constraint, merge_list){

	if(nrow(merge_list) == 0){
		constraint[1,2] <- T
		constraint[2,1] <- T
		return(constraint)
	}else{
		constraint[merge_list[,1],] <- constraint[merge_list[,1],] + constraint[merge_list[,2],]
		constraint[,merge_list[,1]] <- constraint[,merge_list[,1]] + constraint[,merge_list[,2]]
		as.matrix((constraint[-merge_list[,2], -merge_list[,2]] > 0) * 1)
	}
}

#' Update the data frame of observations using the current merge list.
#' @param df the current state of the df
#' @param merge_list the list of observations to be merged.
#' @return df the updated data frame
update_df <- function(df, merge_list){

	df[merge_list[,1],] <- df[merge_list[,1],] + df[merge_list[,2],]
	df[-merge_list[,2],]

}

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------

#' Initialize an empty hclust object to be populated by the loop in info_clust
#' @param n the number of leaves
#' @return a a list with entries appropriate for conversion to an hclust object.
empty_hclust <- function(n){
	a <- list()
	a$merge <- matrix(0, 1, 2)
	a$merge <- a$merge[-1,]
	a$height <- numeric()
	a$order <- 1:n
	a$labels <- 1:n
	a
}

#' Create an empty merge lookup for use populating a$merge
#' @param n the number of leaves
#' @return merge_lookup a list
empty_merge_lookup <- function(n){

	merge_lookup <- list()
	for(i in 1:n){
		merge_lookup[[i]] <- -i
	}
	merge_lookup
}



#' Sort the full list of merges produced in the main loop to ensure the greediness of the algorithm. While this version of agglomerative clustering is fast, the merge list generated is not on its own greedy. We can greedify the algorithm using the sorting procedure below, ensuring that each entry of a$height is minimal given the previous merges. *Sorting the list does not change the cluster topology*, just the order of the merges in the tree.
#' @param merges the data frame of merges
#' @return sorted the sorted data frame of merges

sort_merges <- function(merges){

	# merges <- merges %>% tbl_df

	sorted <- data.frame(i = integer(),
						 j = integer(),
						 loss = numeric(),
						 original_order = integer())

	merges <- dplyr::mutate(merges,
							original_order = row_number(),
							has_i = (i < 0 | i %in% sorted$original_order),
							has_j = (j < 0 | j %in% sorted$original_order),
							predecessors = (has_i & has_j))

	# updates
	while(nrow(merges) > 0){
		merges <- merges %>%
			filter(!(i %in% sorted$i)) %>%
			mutate(has_i = i < 0 | i %in% sorted$original_order,
				   has_j = j < 0 | j %in% sorted$original_order,
				   predecessors = has_i & has_j)

		update <- merges %>%
			filter(predecessors) %>%
			filter(loss == min(loss)) %>%
			select(i, j, loss, original_order)


		sorted <- rbind(sorted, update)
	}

	order_lookup <- sorted %>%
		select(original_order) %>%
		mutate(new_order = row_number()) %>%
		arrange(original_order) %>%
		select(new_order)

	order_lookup <- order_lookup$new_order

	sorted$i[sorted$i > 0] <- order_lookup[sorted$i[sorted$i > 0]]
	sorted$j[sorted$j > 0] <- order_lookup[sorted$j[sorted$j > 0]]
	rownames(sorted) <- 1:nrow(sorted)
	sorted

}

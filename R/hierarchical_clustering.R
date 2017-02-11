#' hierarchical_clustering
#' @name hierarchical_clustering
#' @docType package
#' @import tidyverse igraph
NULL

#' Compute the merge list for a single pass of agglomerative hierarchical clustering. 
#' The result is a df in which each row corresponds to a pair of nodes
#' and a pair is present if each node is the nearest neighbor
#' of the other. Also includes stage information (labels) 
#' and the loss associated with agglomerating the two nodes. 
#' @param h the graph to cluster. Must have "stage" information for the nodes
#' @param current_stage the number of agglomerations performed so far
#' @param divergence the divergence function used to compute distances between nodes
#' @return a data frame giving the merges and loss associated with each one
make_merge_list <- function(h, divergence, current_stage = 0){

	h %>%
		as_long_data_frame() %>%
		tbl_df() %>%
		select(from, to, from_n, to_n, from_stage, to_stage) %>%
		mutate(div = map2_dbl(from_n, to_n, divergence)) %>%
		group_by(from) %>%
		mutate(smallest_from = div == min(div)) %>%
		group_by(to) %>%
		mutate(smallest_to = div == min(div),
			   to_merge = smallest_from & smallest_to) %>%
		ungroup() %>%
		filter(to_merge) %>%
		select(from, to, from_stage, to_stage, div) %>%
		mutate(pair = row_number()) %>%
		gather(key = dir, value = name, -pair, -div, -from_stage, -to_stage) %>%
		distinct(name, .keep_all = T) %>%
		group_by(pair) %>%
		filter(n() == 2) %>%
		ungroup() %>%
		spread(key = dir, value = name) %>%
		select(-pair) %>%
		arrange(div) %>%
		mutate(new_id = row_number(),
			   new_stage = new_id + current_stage)
}

#' Merge nodes in a graph based on an input merge list
#' @param h the graph to merge
#' @param merge_list a data frame giving the merges and
#' losses associated with each merge
#' @return a graph h with specified nodes merged
merge_graph <- function(h, merge_list){
	offset <- max(merge_list$new_id) + 1

	map_vec <- rep(0, length(V(h)))
	for(i in 1:nrow(merge_list)){
		this_merge <- merge_list[i,]
		map_vec[this_merge$from] <- this_merge$new_id
		map_vec[this_merge$to] <- this_merge$new_id
	}
	map_vec[map_vec == 0] <- offset:(offset + sum(map_vec == 0)-1)


	h <- h %>% contract.vertices(map_vec, vertex.attr.comb = list(name = 'concat',
															 n = function(n) reduce(n, .f = `+`),
															 stage = function(stage) reduce(stage, .f = `*`),
															 'ignore')) %>%
		simplify(remove.multiple = T)

	h <- h %>% set_vertex_attr('stage',
						  index = V(h)[1:nrow(merge_list)],
						  value = merge_list$new_stage)

	h
}

#' Construct an empty object of class 'hclust'
#' @param n the number of objects to be clustered
empty_hclust <- function(n){
	a <- list()
	a$merge <- matrix(0, 1, 2)
	a$merge <- a$merge[-1,]
	a$height <- numeric()
	a$order <- 1:n
	a$labels <- 1:n
	a
}

#' Sort a merge list so that the corresponding hierarchical clustering is indeed greedy. 
#' @param merges the merge list to sort
#' @return the merge list sorted, with appropriately relabeled indices. 
sort_merges <- function(merges){

	merges <- merges %>%
		rename(i = from_stage, j = to_stage, loss = div)

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
	sorted %>%
		rename(from_stage = i, to_stage = j, div = loss)

}

#' Cluster a graph using information-maximizing agglomerative clustering
#' and a custom, user-specified divergence function.
#' @param g the graph to cluster
#' @param divergence the divergence function to use for clustering. 
#' Must take as an argument a pair of node attributes n. 
#' @return an object of class `hclust` summarising the clustering, including
#' cluster labels and the cluster heights. Suitable for use with `cutree` and 
#' all other methods used with `hclust` objects, although `plot` is a bit janky
#' because we haven't implemented ordering. 
#' @export
info_cluster <- function(g, divergence){
	a <- empty_hclust(length(V(g)))
	h <- g

	h <- h %>%
		set_vertex_attr(name = 'stage', value = -(1:length(V(h))))

	label_lookup <- data_frame(name = V(g)$name,
							   id = 1:length(V(g)))

	all_merges <- data_frame(from_stage = integer(), to_stage = integer(), div = numeric())
	# construct a data frame with columns for merge names
	# and height, should be basically simple...?
	current_stage <- 0
	while(length(V(h)) > 1){
		merge_list <- make_merge_list(h, divergence, current_stage)
		h <- merge_graph(h, merge_list)
		current_stage <- current_stage + nrow(merge_list)
		to_add <- merge_list %>%
			select(from_stage, to_stage, div)
		all_merges <- all_merges %>% rbind(to_add)
	}

	all_merges <- sort_merges(all_merges)
	a$merge <- all_merges[,c('from_stage', 'to_stage')] %>% as.matrix()
	a$height <- all_merges$div %>% unlist() %>% cumsum()
	a$labels <- V(g)$name
	class(a) <- 'hclust'
	a
}

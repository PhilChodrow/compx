#' smoother
#' @name smoother
#' @docType package
#' @import rgdal tidyverse rgeos

NULL

#
# #' Spatially smooth a set of vector attributes using a Gaussian radial basis
# #' function spatial smoothing kernel.
# #' @export
#
# rbf_smoother <- function(data, tracts, sigma, local = F){
#
# 	tracts <- tracts[tracts@data$GEOID %in% data$tract,]
#     # compute geographic distances
#     ids <- compx::id_lookup(tracts, key_col = 'GEOID')
#     coords <- coords_df(tracts, km = T)
#     centroids <- coords %>%
#         select(-geoid) %>%
#         as.matrix()
#     rownames(centroids) <- coords$geoid
#     dists <- dist(centroids) %>%
#         as.matrix()
#
#     # construct the kernel matrix
#     K <- exp(- dists^2 / (2 * sigma))
#     if(local){
#     	adj <- gRelate(tracts, byid = TRUE, pattern = '****1****')
#     	K[adj == 0] <- 0
#     }
#
#     # Convert data to matrix form
#     P_df <- data %>%
#         spread(key = group, value = n, fill = 0)
#     P <- P_df %>% select(-tract) %>%
#         as.matrix()
#     row.names(P) <- P_df$tract
#     P <- P[row.names(K),] # to ensure the orders match
#
#     # compute the smoother and normalize
#     smoothed_P <- K %*% P
#     smoothed_P <- (smoothed_P / rowSums(smoothed_P)) * rowSums(P)
#
#     # convert back into data frame form
#     smoothed_df <- smoothed_P %>%
#         as.data.frame() %>%
#         rownames_to_column('tract') %>%
#         gather(key = group, value = n, -tract) %>%
#         tbl_df()
# }

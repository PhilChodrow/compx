library(compx, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(rgdal, quietly = TRUE, warn.conflicts = FALSE)
library(rgeos, quietly = TRUE, warn.conflicts = FALSE)

context('Clustering')

races <- c('Asian', 'Black', 'Hispanic', 'Other', 'White')

# load('data/boston.RData')
# load('data/constraint.RData')

h <- info_clust(tracts@data[,races], constraint = constraint)

testthat::test_that("Clustering correctly measures information loss", {
	expect_equal(tracts@data[,races] %>% mutual_info(), h$height[length(h$height)])
})

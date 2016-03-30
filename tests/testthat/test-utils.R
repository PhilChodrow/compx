library(compx)

context('Test data preparation utils')
states   <- 'MA'
counties <- c(025)
acs_data <- get_census_data(state = states, counties = counties, table_num = "B19001")
data <- format_data(acs_data)

test_that("Data has appropriate length",{
	expect_equal(length(data), 2)
})

test_that("Spatial and distribution data have same number of rows",{
	expect_equal(nrow(data$P), nrow(data$X))
})

test_that("No NAs",{
	expect_equal(sum(is.na(data$P)), 0)
})

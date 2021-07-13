#Error statements
test_that("missing tree", {
  expect_error(
    mapsce(example_cn, example_ccf, example_mutational_ccf, tree = NULL),
    "missing tree"
    )
})


test_that("missing CN", {
  expect_error(
    mapsce(copy_number = NULL, example_ccf, example_mutational_ccf, example_tree),
    "missing copy number"
  )
})

test_that("CN with NA values", {
  expect_error(
    mapsce(copy_number = c(example_cn[-1], NA), example_ccf, example_mutational_ccf, example_tree),
    "copy number NA"
  )
})


test_that("CN with NA values", {
  expect_error(
    mapsce(example_cn, as.matrix(example_ccf[,1]), example_mutational_ccf, example_tree),
    "requires multi region data"
  )
})

test_that("mismatch in CN values vs cluster ccf", {
  expect_error(
    mapsce(example_cn, example_ccf[,-1], example_mutational_ccf, example_tree),
    "mismatch in number of observed copy numbers vs number of regions"
  )
})

# Mapsce example
test_that("correct example mapsce output", {
  expect_s3_class(
    mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree),
    "tbl_df"
  )
})

# Mapsce switching to mapsce for 2 regions
test_that("correct example mapsce output", {
  expect_s3_class(
    mapsce(example_cn_2r, example_ccf_2r, example_mutational_ccf, example_tree_2r),
    "tbl_df"
  )
})

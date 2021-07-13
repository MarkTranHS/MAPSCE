#Error statements
test_that("missing tree", {
  expect_error(
    mapsce2r(example_cn_2r, example_ccf_2r, tree = NULL),
    "missing tree"
  )
})


test_that("missing CN", {
  expect_error(
    mapsce2r(copy_number = NULL, example_ccf_2r, example_tree_2r),
    "missing copy number"
  )
})

test_that("CN with NA values", {
  expect_error(
    mapsce2r(copy_number = c(example_cn_2r[-1], NA), example_ccf_2r, example_tree_2r),
    "copy number NA"
  )
})


test_that("CN with NA values", {
  expect_error(
    mapsce2r(example_cn_2r, as.matrix(example_ccf_2r[,1]), example_tree_2r),
    "requires multi region data"
  )
})

test_that("mismatch in CN values vs cluster ccf", {
  expect_error(
    mapsce2r(example_cn_2r[-1], example_ccf_2r, example_tree_2r),
    "mismatch in number of observed copy numbers vs number of regions"
  )
})

test_that("more than 2 regions", {
  expect_error(
    mapsce2r(example_cn, example_ccf, example_tree_2r),
    "this mapsce mode needs 2 regions only"
  )
})

# mapsce2r example
test_that("correct example mapsce2r output", {
  expect_s3_class(
    mapsce2r(example_cn_2r, example_ccf_2r, example_tree_2r),
    "tbl_df"
  )
})


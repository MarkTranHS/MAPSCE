#Error statements
test_that("missing tree", {
  expect_error(
    get_consensus(mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree), NULL),
    "missing tree"
  )
})


test_that("data frame as input", {
  expect_error(
    get_consensus(NULL, example_tree),
    "needs data frame/tibble as input data"
  )
})

test_that("mapsce output as input", {
  expect_error(
    get_consensus(mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree)[,-1], example_tree),
    "number of columns in data doesn't match mapsce output"
  )
})

test_that("mismatch between data input and tree", {
  expect_error(
      get_consensus(mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree), example_tree_2r),
    "mismatch between data input and number of tree nodes"
  )
})


# get_consensus example
test_that("correct example mapsce output", {
  expect_type(
    get_consensus(mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree), example_tree),
    "double"
  )
})

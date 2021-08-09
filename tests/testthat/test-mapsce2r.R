#Testing error statements
test_that("missing tree", {
  expect_error(
    mapsce2r(example_cn_2r, example_ccf_2r, tree = NULL),
    "missing tree"
  )
  expect_error(
    mapsce2r(copy_number = NULL, example_ccf_2r, example_tree_2r),
    "missing copy number"
  )
  expect_error(
    mapsce2r(copy_number = c(example_cn_2r[-1], NA), example_ccf_2r, example_tree_2r),
    "copy number NA"
  )
  expect_error(
    mapsce2r(example_cn_2r, as.matrix(example_ccf_2r[,1]), example_tree_2r),
    "requires multi region data"
  )
  expect_error(
    mapsce2r(example_cn_2r[-1], example_ccf_2r, example_tree_2r),
    "mismatch in number of observed copy numbers vs number of regions"
  )
  expect_error(
    mapsce2r(example_cn, example_ccf, example_tree_2r),
    "this mapsce mode needs 2 regions only"
  )
})

# Testing mapsce2r example
test_that("correct example mapsce2r output", {
  mapsce_result_2r <- mapsce2r(example_cn_2r, example_ccf_2r, example_tree_2r)
  expect_s3_class(
    mapsce_result_2r,
    "tbl_df"
  ) #output class for mapsce2r example

  expect_true(
    all(sort(colnames(mapsce_result_2r)) == c("after", "before", "bf",
                                     "bic", "branch", "btn",
                                     "evid", "index", "nclones",
                                     "nregions", "null", "rss")
    )
  ) # checking all colnames of mapsce2r output

  expect_equal(
    mapsce_result_2r %>% dplyr::pull(branch) %>% sort(),
    sort(unique(as.vector(example_tree_2r)))
  ) # checking  clones in the tree vs output

  expect_equal(
    all(is.na(mapsce_result_2r %>% dplyr::pull(btn))),
    TRUE
  ) # checking whether all btns are NA for mapsce2r

  expect_equal(
    mapsce_result_2r %>% dplyr::filter(before == after) %>% dplyr::pull(null),
    "yes"
  ) # checking null conditions
  expect_true(
    mapsce_result_2r %>% dplyr::filter(null == "yes") %>% dplyr::pull(evid) %in% c(0,2)
  ) # checking null conditions

  expect_true(
    any(mapsce_result_2r %>% dplyr::filter(evid == 1 | evid == 2) %>% dplyr::pull(index) == 1)
  ) #checking whether best result is always included
})

test_that("toy_examples", {
  this_cn <- c(0.9, 0.02)
  this_ccf <- matrix(c(100, 0, 100, 100), nrow = 2)
  dimnames(this_ccf) <- list(c("1", "2"), c("R1", "R2"))
  this_tree <- matrix(c("1", "2"), nrow = 1)
  res <- mapsce2r(this_cn, this_ccf, this_tree)
  expect_equal(
    res$branch[1], "2"
  )

  this_cn <- c(0.98, 0.5)
  this_ccf <- matrix(c(100, 80, 0, 100, 60, 50), nrow = 3)
  dimnames(this_ccf) <- list(c("1", "2", "3"), c("R1", "R2"))
  this_tree <- matrix(c("1", "2", "2", "3"), nrow = 2)
  res <- mapsce2r(this_cn, this_ccf, this_tree)
  expect_equal(
    res$branch[1], "3"
  )
})

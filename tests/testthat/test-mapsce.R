#Testing error statements
test_that("stop errors", {
  expect_error(
    mapsce(example_cn, example_ccf, example_mutational_ccf, tree = NULL),
    "missing tree"
    )
  expect_error(
    mapsce(copy_number = NULL, example_ccf, example_mutational_ccf, example_tree),
    "missing copy number"
  )
  expect_error(
    mapsce(copy_number = c(example_cn[-1], NA), example_ccf, example_mutational_ccf, example_tree),
    "copy number NA"
  )
  expect_error(
    mapsce(example_cn, as.matrix(example_ccf[,1]), example_mutational_ccf, example_tree),
    "requires multi region data"
  )
  expect_error(
    mapsce(example_cn, example_ccf[,-1], example_mutational_ccf, example_tree),
    "mismatch in number of observed copy numbers vs number of regions"
  )
})

# Testing mapsce example
test_that("correct example mapsce output", {
  mapsce_result <- mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree)
  mapsce_result_2r <- mapsce(example_cn_2r, example_ccf_2r, example_mutational_ccf, example_tree_2r)
  expect_s3_class(
    mapsce_result,
    "tbl_df"
  ) # output class for mapsce example
  expect_s3_class(
    mapsce_result_2r,
    "tbl_df"
  ) # output class for mapsce switching to mapsce for 2 regions

  expect_true(
    all(sort(colnames(mapsce_result)) == c("after", "before", "bf",
                                        "bic", "branch", "btn",
                                        "evid", "index", "nclones",
                                        "nregions", "null", "rss")
    )
  ) # checking all colnames of mapsce output

  expect_equal(
    mapsce_result %>% dplyr::pull(branch) %>% sort(),
    sort(unique(as.vector(example_tree)))
  ) # checking  clones in the tree vs output

  expect_equal(
    mapsce_result %>% dplyr::filter(null == "yes") %>% dplyr::pull(btn),
    101
  ) # checking null conditions
  expect_equal(
    mapsce_result %>% dplyr::filter(before == after) %>% dplyr::pull(null),
    "yes"
  ) # checking null conditions
  expect_true(
    mapsce_result %>% dplyr::filter(null == "yes") %>% dplyr::pull(evid) %in% c(0,2)
  ) # checking null conditions

  expect_true(
    any(mapsce_result %>% dplyr::filter(evid == 1 | evid == 2) %>% dplyr::pull(index) == 1)
  ) #checking whether best result is always included
})


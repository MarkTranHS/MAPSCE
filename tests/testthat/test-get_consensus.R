#Testing error statements
test_that("error stops for input data", {
  mapsce_result <- mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree, print_duration = F)

  expect_error(
    get_consensus(mapsce_result, NULL),
    "missing tree"
  )
  expect_error(
    get_consensus(NULL, example_tree),
    "needs data frame/tibble as input data"
  )
  expect_error(
    get_consensus(mapsce_result[,-1], example_tree),
    "number of columns in data doesn't match mapsce output"
  )
  expect_error(
    get_consensus(mapsce_result, example_tree_2r),
    "mismatch between data input and number of tree nodes"
  )
})

#Testing the example and correct output
test_that("correct example mapsce output", {
  mapsce_result <- mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree, print_duration = F)
  consensus_result <- get_consensus(mapsce_result, example_tree)
  consensus_result_only <- get_consensus(mapsce_result, example_tree, consensus.only = T)

  expect_type(
    consensus_result,
    "double"
  ) #checking the type of consensus_result
  expect_equal(
    class(consensus_result),
    c("matrix", "array")
  )

  expect_equal(
    nrow(mapsce_result),
    ncol(consensus_result)
  ) #checking number of clones match up between mapsce and consensus results
  expect_equal(
    ncol(consensus_result),
    length(unique(as.vector(example_tree)))
  ) #checking number of clones match up between tree and consensus result

  expect_equal(
    rownames(consensus_result)[nrow(consensus_result)],
    "consensus"
  ) #checking last row of consensus result is the consensus

  expect_equal(
    mapsce_result %>% dplyr::filter(evid == 1) %>% nrow(),
    nrow(consensus_result) - 1
  ) #checking good results number vs number of results used in get_consensus

  expect_equal(
    is.vector(consensus_result_only),
    TRUE
  ) #checking consensus.only output

  consensus_result <- get_consensus(mapsce_result, example_tree, consensus = F)

  expect_equal(
    length(consensus_result),
    length(unique(as.vector(example_tree)))
  ) #checking results when consensus is FALSE
  expect_equal(
    mapsce_result %>% dplyr::filter(evid == 1) %>% nrow(),
    nrow(consensus_result)
  ) #checking results when consensus is FALSE
})




#' mapsce2r
#'
#' `mapsce2r` returns a matrix with a branch test for every single branch of a tree
#'
#' @param copy_number a numeric vector with copy number values ordered by sample
#' @param cluster_ccf a matrix with mean CCF values for clones in rows and samples (same order as in copy_number) in columns
#' @param tree a matrix listing all the branches in the tree, where the first column is the ancestral node and the second column is the descendant clone. By definition, the root node will only be present in the first column. The clone IDs must correspond to the cluster IDs in cluster_ccf and mutation_ccf.
#' @param print_raw_matrix logical, printing of raw results
#' @param print_duration logical, printing of the time taken to run
#' @param clone_ccf logical, forcing mapsce to use clone CCF instead of just cluster CCF
#' @param consensus logical, forcing mapsce to run the get_consensus function
#' @param before a numeric vector, forcing a CN before
#' @param after a numeric vector, forcing a CN after
#' @return a tibble with column names:
#' \itemize{
#'   \item branch. branch ID
#'   \item before. copy number state before
#'   \item after. copy number state after
#'   \item rss. residual square of sums
#'   \item btn. how many bootstrapped BICs are better than the null's BI
#'   \item nregions. number of regions
#'   \item nclones. number of clones or branches.
#'   \item null. whether the branch is a trunk or not.
#'   \item bic.  bic calculated from mean RSS
#'   \item bf. the strength of evidence of bayes factor comparison
#'   \item good_result. 1 for good results based on the Bayes Factor comparison
#'   \item index. ordering of results.
#' }
#'
#'@examples
#'data(example_data)
#'mapsce2r(example_cn_2r, example_ccf_2r, example_tree_2r)
#'#returns a tibble with branch 2 as the best result
#'
#' @export
#'

## MAPSCE algorithm
mapsce2r <- function(copy_number,
                     cluster_ccf,
                     tree,
                     print_raw_matrix = F,
                     print_duration = T,
                     clone_ccf = F,
                     consensus = F,
                     before = 0,
                     after = 0){
  start.time <- Sys.time() #timing

  #Stop conditions
  if(length(tree) == 0) {
    stop("missing tree")
  }
  if(length(copy_number) == 0 ){
    stop("missing copy number")
  }
  if(any(is.na(copy_number))){
    stop("copy number NA")
  }
  if(ncol(cluster_ccf)==1){
    stop("requires multi region data")
  }
  if(ncol(cluster_ccf) != length(copy_number)){
    stop("mismatch in number of observed copy numbers vs number of regions")
  }
  if(length(copy_number) != 2){
    stop("this mapsce mode needs 2 regions only")
  }

  #List of all clone IDs from tree - continue results
  clones <- unique(c(tree[, 1], tree[, 2]))
  number_of_clones <- length(clones)
  nregions <- (ncol(cluster_ccf))
  results_matrix <- matrix(ncol=7, nrow=number_of_clones) #preparing matrix of results
  colnames(results_matrix) <- c("branch", "before", "after", "rss", "nregions", "nclones", "null")
  results_matrix <- as.data.frame(results_matrix)


  row_number <- 1 #indexing the row number

  pyclone_ccf <- cluster_ccf

  #Get recuded ccf and clone fraction
  reduced_pyclone_ccf <- reduce_pyclone_ccf(pyclone_ccf, tree)
  clone_fraction <- get_clone_fraction_from_ccf(reduced_pyclone_ccf, tree)

  #Set non-negative restraint
  clone_fraction[clone_fraction < 0] <- 0
  #print(clone_fraction[1,])
  nested_set <- get_nested_set_from_tree(tree)

  #Make copy_number non-negative
  copy_number[copy_number<0] <- 0

  #Branch testing with QP
  for (this_clone in clones) {
    if(all(reduced_pyclone_ccf[this_clone,] == 0)){
      print(paste("It's clone number", this_clone))
      print("This clone has no mutations in boostrapped clusterCCF")
      results_matrix <- results_matrix[-row_number,]
      next()
    }

    # From the tree object, we read the clone IDs as strings. We want to use
    # the names of the columns in the nested_set matrix
    subtree <- c(this_clone, names(which(nested_set[this_clone, ] == 1)))

    # Null hypothesis testing
    if (length(subtree) == number_of_clones) {
      new_matrix <- matrix(rep(1, ncol(clone_fraction)), 1)
      Dmat <- new_matrix %*% t(new_matrix)
      dvec <- copy_number %*% t(new_matrix)
      Amat <- diag(1)
      bvec <- c(after)
      k <- 1 #nr of parameters
      null <- "yes"
    } else {
      rest.of.the.tree <- setdiff(clones, subtree)

      if(clone_ccf == T){ #run clone CCF instead of using only cluster CCF
        new_matrix <- rbind(colSums(matrix(clone_fraction[subtree, ],
                                           length(subtree), ncol(clone_fraction))),
                            colSums(matrix(clone_fraction[rest.of.the.tree, ],
                                           number_of_clones - length(subtree), ncol(clone_fraction))))
      } else {
        new_matrix <- rbind(reduced_pyclone_ccf[this_clone,],
                            100 - reduced_pyclone_ccf[this_clone,])
      }

      new_matrix <- new_matrix / 100 # Divide by 100 as the CCF are in 100%

      # Run QP with restriction that all values must be positive (or equal to 0)
      Dmat <- new_matrix %*% t(new_matrix)
      dvec <- copy_number %*% t(new_matrix)
      Amat <- diag(2) # All positive
      bvec <- c(after,before)
      k <- 2 # Number of parameters
      null <- "no"
    }

    sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec))

    if (class(sol) != "try-error") {

      predicted <- t(new_matrix) %*% sol$solution
      diff <- (copy_number - t(new_matrix) %*% sol$solution)
      rsum <- sum(diff ^ 2)

      if(length(sol$solution) == 1) {sol$solution <- c(sol$solution, sol$solution)}

      results_matrix[row_number,1] <- this_clone #clone
      results_matrix[row_number,2] <- sol$solution[2] #before
      results_matrix[row_number,3] <- sol$solution[1] #after
      results_matrix[row_number,4] <- rsum #rss
      # results_matrix[row_number,5] <- as.numeric(length(copy_number) * log(rsum/length(copy_number)) +
      #                                              k*log(length(copy_number))) #bic
      results_matrix[row_number,5] <- nregions #number of regions
      results_matrix[row_number,6] <- number_of_clones #number of clones
      results_matrix[row_number,7] <- null #null
      row_number <- row_number+1
    }
  }

  results_matrix <- results_matrix %>% dplyr::as_tibble()

  # # Summarising results
  # null_bic <- results_matrix %>%
  #   dplyr::filter(null == "yes") %>%
  #   dplyr::pull(bic)
  # null_bic <- rep(null_bic, nrow(results_matrix))
  # null_bic <- tibble::tibble(null_bic=null_bic)
  #
  # summarised_results <- results_matrix %>%
  #   cbind(null_bic) %>%
  #   dplyr::group_by(branch) %>%
  #   dplyr::mutate(bic = as.numeric(bic)) %>%
  #   dplyr::mutate(bic = ifelse(bic==-Inf, -1*10^100, bic)) %>%
  #   dplyr::arrange(bic)
  #
  # ## Bayes factors evidence
  # all_branches <- summarised_results %>%
  #   dplyr::pull(branch) %>%
  #   unique %>%
  #   as.character()
  # top_bic <- tibble::tibble(top_bic = summarised_results %>% .[1,] %>%
  #                             dplyr::pull(bic) %>%
  #                             rep(., length(all_branches)))
  #
  # summarised_results <- summarised_results %>%
  #   dplyr::ungroup() %>%
  #   cbind(top_bic) %>%
  #   dplyr::as_tibble() %>%
  #   dplyr::arrange(bic) %>%
  #   dplyr::mutate(bic_diff = exp((bic - top_bic)/2)) %>%
  #   dplyr::mutate(bf = ifelse(bic_diff < 6, 1, 0))
  #
  #
  # ## Adding strength of evidence
  # summarised_results <- summarised_results %>%
  #   dplyr::mutate(good_result = ifelse(bf == 1, 1, 0))

  # ## Not rejected null
  # all_evidence <- summarised_results %>% dplyr::pull(bf)
  # if(all(all_evidence == 1)){
  #   summarised_results <- summarised_results %>% dplyr::mutate(good_result = ifelse(null == "yes", 2, 0))
  # } else {
  #   if(!rlang::is_empty(intersect(which(all_evidence == 1), which(summarised_results %>% dplyr::pull(null) == "yes")))){
  #     summarised_results <- summarised_results %>% dplyr::mutate(good_result = ifelse(null == "yes", 2, 0))
  #   }
  # }



  if(print_raw_matrix==T){print(summarised_results)}
  if(print_duration==T){print(Sys.time()-start.time)}

  summarised_results <- results_matrix %>%
    dplyr::mutate(before = before,
                  after = after,
                  rss = rss,
                  nregions = nregions,
                  nclones = nclones) %>%
    dplyr::arrange(rss) %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    dplyr::mutate(good_result = ifelse(index == 1, 1, 0))

  summarised_results <- summarised_results %>%
    dplyr::arrange(desc(good_result), rss)

  summarised_results <- summarised_results %>%
    dplyr::mutate(before = before,
                    after = after,
                    rss = rss,
                    btn = NA,
                    bic = NA,
                    nregions = nregions,
                    nclones = nclones,
                    bf = NA,
                    good_result = good_result,
                    index = index
    )
  if(consensus == T){
    consensus_mapping <- get_consensus(summarised_results, tree)
    v_viability <- sum(!is.na(consensus_mapping["consensus",]))
    no_results <-  nrow(consensus_mapping) - 1
    nclones <- nrow(summarised_results)
    top <- F
    if(v_viability < nclones - no_results){
      if(all(summarised_results$good_result == 2)){
      summarised_results[2:nrow(summarised_results), "good_result"] <- 0
      top <- T
      }
    }
    summarised_results <- summarised_results %>%
      dplyr::mutate(viability = v_viability, nresults = no_results, top_only = top)
  }

  return(summarised_results)
}

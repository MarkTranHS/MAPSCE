#' mapsce
#'
#' `mapsce` returns a matrix with a branch test for every single branch of a tree
#'
#' @param copy_number a numeric vector with copy number values ordered by sample
#' @param cluster_ccf a matrix with mean CCF values for clones in rows and samples (same order as in copy_number) in columns
#' @param mutation_ccf a tibble or data frame where the first columns contain the mutational CCF values for the  samples (same order as in copy_number) in columns. Tw of the columns in this object need to be "PycloneCluster" and "CleanCluster", the former being the clone assignment and the latter being the PyClone filter (1 for valid clusters, 0 for clusters to be ignored)
#' @param tree a matrix listing all the branches in the tree, where the first column is the ancestral node and the second column is the descendant clone. By definition, the root node will only be present in the first column. The clone IDs must correspond to the cluster IDs in cluster_ccf and mutation_ccf.
#' @param bootstraps number of bootstraps for reclustering of the CCF
#' @param print_raw_matrix printing of raw results
#' @param print_duration printing of the time taken to run
#' @return a tibble with column names:
#' \itemize{
#'   \item branch. branch ID
#'   \item before. copy number state before
#'   \item after. copy number state after
#'   \item rss. average residual square of sums of all bootstraps
#'   \item btn. how many bootstrapped BICs were better than the null's BIC
#'   \item nregions. number of regions
#'   \item nclones. number of clones or branches.
#'   \item null. whether the branch is a trunk or not.
#'   \item bic. average bic for all bootstraps calculated from mean RSS
#'   \item bf. the strength of evidence of bayes factor comparison
#'   \item evid. 1 for good results based on the Bayes Factor comparison and btn
#'   \item index. ordering of results.
#' }
#'
#'@examples
#'data(example_data)
#'mapsce(example_cn, example_ccf, example_mutational_ccf, example_tree)
#'#returns a tibble with branch 7 as the best result
#'
#' @export
#'

## MAPSCE algorithm
mapsce <- function(copy_number, cluster_ccf, mutation_ccf, tree, bootstraps=100, print_raw_matrix="no", print_duration="yes"){
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
  if(length(copy_number) == 2){
    print("running mapsce2r - mapsce for 2 regions")
    summarised_results <- mapsce2r(copy_number, cluster_ccf, tree)
    return(summarised_results)
  }

  #List of all clone IDs from tree - continue results
  clones <- unique(c(tree[, 1], tree[, 2]))
  number_of_clones <- length(clones)
  nregions <- (ncol(cluster_ccf))
  results_matrix <- matrix(ncol=8, nrow=bootstraps*number_of_clones) #preparing matrix of results
  colnames(results_matrix) <- c("branch", "before", "after", "rss", "bic", "nregions", "nclones", "null")


  row_number <- 1 #indexing the row number

  for(this_bootstrap in 1:bootstraps){

    #Resampling mutations into cluster CCF
    all_clusters <- mutation_ccf %>%
      dplyr::filter(CleanCluster == 1) %>%
      dplyr::pull(PycloneCluster) %>%
      unique() %>%
      sort()
    pyclone_ccf <- matrix(ncol=nregions, nrow=length(all_clusters)) #preparing matrix for resampling
    for (this_cluster in all_clusters) {
      all_mutations <- mutation_ccf %>%
        dplyr::filter(PycloneCluster == this_cluster)
      for (this_region in 1:nregions) {
        pyclone_ccf[this_cluster,this_region] <- c(mean(sample(dplyr::pull(all_mutations[,this_region]),
                                                               nrow(all_mutations), replace=T)))
      }
    }
    pyclone_ccf <- pyclone_ccf*100
    rownames(pyclone_ccf) <- all_clusters
    colnames(pyclone_ccf) <- paste0("R", 1:nregions)

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

      # From the tree object, we read the clone IDs as strings. We want to use
      # the names of the columns in the nested_set matrix
      subtree <- c(this_clone, names(which(nested_set[this_clone, ] == 1)))

      # Null hypothesis testing
      if (length(subtree) == number_of_clones) {
        new_matrix <- matrix(rep(1, ncol(clone_fraction)), 1)
        Dmat <- new_matrix %*% t(new_matrix)
        dvec <- copy_number %*% t(new_matrix)
        Amat <- diag(1)
        bvec <- c(0)
        k <- 1 #nr of parameters
        null <- "yes"
      } else {
        rest.of.the.tree <- setdiff(clones, subtree)

        # New_matrix is a 2-rows matrix with the sum of the clone fractions for the subtree and
        # for the rest of the tree
        new_matrix <- rbind(colSums(matrix(clone_fraction[subtree, ],
                                           length(subtree), ncol(clone_fraction))),
                            colSums(matrix(clone_fraction[rest.of.the.tree, ],
                                           number_of_clones - length(subtree), ncol(clone_fraction))))
        new_matrix <- new_matrix / 100 # Divide by 100 as the CCF are in 100%

        # Run QP with restriction that all values must be positive (or equal to 0)
        Dmat <- new_matrix %*% t(new_matrix)
        dvec <- copy_number %*% t(new_matrix)
        Amat <- diag(2) # All positive
        bvec <- c(0,0)
        k <- 2 # Number of parameters
        null <- "no"
      }

      sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec))

      if (class(sol) != "try-error") {

        predicted <- t(new_matrix) %*% sol$solution
        diff <- (copy_number - t(new_matrix) %*% sol$solution)
        rsum <- sum(diff ^ 2)

        if(length(sol$solution) == 1) {sol$solution <- c(sol$solution, sol$solution)}

        results_matrix[row_number,1] <- as.numeric(this_clone) #clone
        results_matrix[row_number,2] <- as.numeric(sol$solution[2]) #before
        results_matrix[row_number,3] <- as.numeric(sol$solution[1]) #after
        results_matrix[row_number,4] <- as.numeric(rsum) #rss
        results_matrix[row_number,5] <- as.numeric(length(copy_number) * log(rsum/length(copy_number)) +
                                                     k*log(length(copy_number))) #bic
        results_matrix[row_number,6] <- as.numeric(nregions) #number of regions
        results_matrix[row_number,7] <- as.numeric(number_of_clones) #number of clones
        results_matrix[row_number,8] <- null #null
        row_number <- row_number+1
      }
    }
  }
  results_matrix <- results_matrix %>% dplyr::as_tibble()

  # Summarising bootstrapped results
  null_bic <- results_matrix %>%
    dplyr::filter(null == "yes") %>%
    dplyr::pull(bic) %>%
    unique()
  null_bic <- rep(null_bic, nrow(results_matrix))
  null_bic <- tibble::tibble(null_bic=null_bic)

  summarised_results <- results_matrix %>%
    dplyr::as_tibble() %>%
    cbind(null_bic) %>%
    dplyr::group_by(branch) %>%
    dplyr::summarize(before = mean(as.numeric(before), na.rm=T),
              after = mean(as.numeric(after), na.rm=T),
              rss = mean(as.numeric(rss), na.rm=T),
              btn = ifelse(unique(null) == "yes", 101, sum(as.numeric(bic) <= as.numeric(null_bic)) / length(bic) * 100),
              bic = mean(as.numeric(bic), na.rm=T),
              nregions = mean(as.numeric(nregions), na.rm=T),
              nclones = mean(as.numeric(nclones), na.rm=T),
              null = unique(null)) %>%
    dplyr::arrange(bic)

  summarised_results <- summarised_results %>%
    dplyr::select(-bic) %>%
    dplyr::mutate(bic = nregions * log(rss/nregions) + ifelse(null == "yes", 1, 2) * log(nregions))

  summarised_results <- summarised_results %>%
    dplyr::mutate(bic = ifelse(bic==-Inf, -1*10^100, bic))

  ## Bayes factors evidence
  all_branches <- summarised_results %>%
    dplyr::pull(branch) %>%
    unique %>%
    as.character()
  top_bic <- tibble::tibble(top_bic = summarised_results %>% .[1,] %>%
                              dplyr::pull(bic) %>%
                              rep(., length(all_branches)))

  summarised_results <- summarised_results %>%
    cbind(top_bic) %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(bic) %>%
    dplyr::mutate(bic_diff = exp((bic - top_bic)/2)) %>%
    dplyr::mutate(bf = ifelse(bic_diff < 6, 1, 0))

  ## Adding null BIC
  summarised_results <- summarised_results %>%
    cbind(null_bic[1:length(all_branches),]) %>%
    dplyr::as_tibble()

  ## Adding strength of evidence
  summarised_results <- summarised_results %>%
    dplyr::mutate(evid = ifelse(bf == 1, 1, 0)) %>%
    dplyr::mutate(evid = ifelse(btn < 95, 0, evid))

  ## Not rejected null
  all_evidence <- summarised_results %>% dplyr::pull(bf)
  all_btn <- summarised_results %>% dplyr::pull(btn)
  if(all(all_evidence == 1)){
    summarised_results <- summarised_results %>% dplyr::mutate(evid = ifelse(null == "yes", 2, 0))
  } else {
    if(all(all_btn < 95)){
      summarised_results <- summarised_results %>% dplyr::mutate(evid = ifelse(null == "yes", 2, 0))
    } else {
        if(rlang::is_empty(intersect(which(all_evidence == 1), which(all_btn >= 95)))){
          summarised_results <- summarised_results %>% dplyr::mutate(evid = ifelse(null == "yes", 2, 0))
        }
      }
  }


  if(print_raw_matrix=="yes"){print(summarised_results)}
  if(print_duration=="yes"){print(Sys.time()-start.time)}

  summarised_results <- summarised_results %>%
    dplyr::arrange(desc(evid), bic) %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    dplyr::select(-top_bic, -bic_diff, -null_bic)

  return(summarised_results)
}

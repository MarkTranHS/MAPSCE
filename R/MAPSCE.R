#' mapsce
#'
#' `mapsce` returns branch test results.
#'
#' @param patient patient ID
#' @param copy_number observed copy numbers
#' @param cluster_ccf CCF matrix
#' @param mutation_ccf mutational CCF matrix
#' @param tree tree matrix
#' @param bootstraps number of bootstraps to run
#' @param print_raw_matrix logical
#' @param print_duration printing time running
#' @return Returns tibble with summarised bootstrapped mapping results
#'
#' @export
#'

## MAPSCE algorithm
mapsce <- function(patient, copy_number, cluster_ccf, mutation_ccf, tree, bootstraps=100, print_raw_matrix="no", print_duration="yes"){
  start.time <- Sys.time() #timing

  #Stop conditions
  if(length(tree) == 0) {
    stop("missing tree")
  }
  if(length(copy_number) == 0 ){
    stop("missing copy number")
  }
  if(any(is.na(copy_number))){
    stop("missing observed copy number")
  }
  if(ncol(cluster_ccf)==1){
    stop("requires multi region data")
  }
  if(ncol(cluster_ccf) != length(copy_number)){
    stop("mismatch in number of observed copy numbers vs number of regions")
  }

  #Adjusting CCF to match the tree branches if there are less branches than trees
  adjust_clone_fraction_from_ccf <- function(clone_fraction, tree, node) {
    for (descendent in tree[tree[, 1] == node, 2]) {
      clone_fraction[node, ] <- clone_fraction[node, ] - clone_fraction[descendent, ]
      clone_fraction <- adjust_clone_fraction_from_ccf(clone_fraction, tree, descendent)
    }
    return(clone_fraction)
  }
  #Setting trunk from tree to feed into recursive function
  get_clone_fraction_from_ccf <- function(ccf, tree) {
    trunk <- setdiff(tree[, 1], tree[, 2])
    stopifnot(length(trunk) == 1)
    clone_fraction <- adjust_clone_fraction_from_ccf(ccf, tree, trunk)
    return(clone_fraction)
  }
  #Reducing pyclone trees
  reduce_pyclone_ccf <- function(ccf, tree) {
    nodes <- unique(c(tree[, 1], tree[, 2]))
    pyclone_ccf <- ccf[nodes, , drop = FALSE]
    return(pyclone_ccf)
  }
  #Making a nested set matrix - 1/0
  get_nested_set_from_tree <- function(tree) {
    clones <- unique(c(tree[, 1], tree[, 2]))
    number_of_clones <- length(clones)

    nested_set <- matrix(0, number_of_clones, number_of_clones)
    rownames(nested_set) <- clones
    colnames(nested_set) <- clones
    nested_set[tree] <- 1
    n <- 0
    while (n < sum(nested_set)) {
      n <- sum(nested_set)
      ancestors <- list()
      for (i in 1:number_of_clones) {
        ancestors[[i]] <- which(nested_set[, i] == 1)
      }
      for (i in 1:number_of_clones) {
        for (j in ancestors[[i]]) {
          nested_set[ancestors[[j]], i] <- 1
        }
      }
    }
    return(nested_set)
  }

  #List of all clone IDs from tree - continue results
  clones <- unique(c(tree[, 1], tree[, 2]))
  number_of_clones <- length(clones)
  number_of_regions <- (ncol(cluster_ccf))
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
    pyclone_ccf <- matrix(ncol=number_of_regions, nrow=length(all_clusters)) #preparing matrix for resampling
    for (this_cluster in all_clusters) {
      all_mutations <- mutation_ccf %>%
        dplyr::filter(PycloneCluster == this_cluster)
      for (this_region in 1:number_of_regions) {
        pyclone_ccf[this_cluster,this_region] <- c(mean(sample(dplyr::pull(all_mutations[,this_region]),
                                                               nrow(all_mutations), replace=T)))
      }
    }
    pyclone_ccf <- pyclone_ccf*100
    rownames(pyclone_ccf) <- all_clusters
    colnames(pyclone_ccf) <- paste(patient, paste("R", 1:number_of_regions, sep=""))

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
        results_matrix[row_number,6] <- as.numeric(number_of_regions) #number of regions
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
              btn = ifelse(unique(null) == "yes", 101, sum(bic <= null_bic) / length(bic) * 100),
              bic = mean(as.numeric(bic), na.rm=T),
              number_of_regions = mean(as.numeric(number_of_regions), na.rm=T),
              nclones = mean(as.numeric(nclones), na.rm=T),
              null = unique(null)) %>%
    dplyr::arrange(bic)

  summarised_results <- summarised_results %>%
    dplyr::select(-bic) %>%
    dplyr::mutate(bic = number_of_regions * log(rss/number_of_regions) + ifelse(null == "yes", 1, 2) * log(number_of_regions))

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
    dplyr::mutate(index = dplyr::row_number())

  return(summarised_results)
}


#https://debruine.github.io/tutorials/your-first-r-package-with-unit-tests.html
#cannot do 1.4.3 Imports
#usethis::use_pipe() leads to error
#✓ Setting active project to '/Users/marktran/Desktop/R/MAPSCE'
#Error: Project 'MAPSCE' does not use roxygen2.
#`use_pipe()` can not work without it.
#You might just need to run `devtools::document()` once, then try again.

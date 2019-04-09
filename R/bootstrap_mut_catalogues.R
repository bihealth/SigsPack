#' Bootstraps a given mutational catalogue
#'
#' Bootstraps a given mutational catalogue  by replicating samples from the
#' original catalogue's distribution of mutational features. The output can be
#' input to signature_exposure.
#'
#' @param n Amount of bootstrapped replicates that are created (they will by
#' default have the same amount of mutations as the original catalogue)
#' @param original Mutational catalogue (matrix) of a sample that is taken as
#' the distribution from which the replicates are sampled
#' @param m Amount of mutations the replicates are supposed to have
#' (e.g. if this differs from the original or if the original is provided as
#' probabilities instead of total counts)
#'
#' @return matrix containing the mutational catalogues of the replicates
#'
#' @examples
#' data(cosmicSigs)
#' reps <- bootstrap_mut_catalogues(n = 150, original = create_mut_catalogues(
#'                                   10, 500)[["catalogues"]][,1])
#'
#' @export
bootstrap_mut_catalogues <- function(n, original, m = NULL){
  if (!is.numeric(n) || n < 1){
    stop('n has to be a positive number', call. = TRUE)
  }

  if (! is.numeric(original)){
    stop('The original mutational catalogue is expected to contain only
         numbers.'
         ,call. = TRUE)
  }

  if(is.null(m)){
    m <- sum(original)
  }
  else if (!is.numeric(n) || n < 1){
    stop('m has to be a positive number', call. = TRUE)
  }
  #=============================================================================
  
  nfeatures <- length(original)
  boot_replicates <- matrix(0, ncol = n, nrow = nfeatures)
  original <- original/m

  for(i in seq_len(n)){
    # partly from SignatureEstimation
    mutations_sampled <- sample(seq(nfeatures), m, replace = TRUE,
                                prob = original)
    rep <- as.numeric(table(factor(mutations_sampled,
                                        levels = seq(nfeatures))))
    boot_replicates[,i] <- rep
  }
  
  rep_names <- paste0('replicate_', seq_len(n))
  colnames(boot_replicates) <- rep_names
  rownames(boot_replicates) <- rownames(original)
  return(boot_replicates)
}

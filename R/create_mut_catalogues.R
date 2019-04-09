#' Creates simulated mutational catalogues
#'
#' Creates mutational catalogues from chosen mutational signature profiles to
#' gain simualted data sampled from a distirbution with with known signature
#' contribution
#'
#' @param n Amount of mutational catalogues that will be created
#' @param m Amount of mutations that each catalogues will have
#' @param P Matrix (f x k) containing the signature profiles of k signatures
#' which will be used to create the catalogues
#' k is the amount of signature profiles that P contains
#' f is the amount of features that the profiles are defined on
#' (default: COSMIC signature matrix 96 x 30) (optional)
#' @param sig_set Numeric vector containing the index of the columns from P
#' (which are the signature profiles) that will be used to create the catalogues
#' (optional). Defaults to all columns of P.
#' @param c_exposure Numeric vector specifying the contribution of each 
#' signature to the distribution each sample is drawn from. Samples will have an
#' appr. exposure of c_exposure[idx] to signature sig_set[idx] 
#' (if sum(c_exposure)==1). (optional)
#' Default: random exposure.
#'
#' @return List containing a matrix (f x n) with the simulated catalogues and a
#' matrix (k x n) detailing the signature exposure of each catalogue
#'
#' @examples
#' data(cosmicSigs)
#' sim_data <- create_mut_catalogues(10, 300, cosmicSigs, c(2,6,15,27))
#' sim_data <- create_mut_catalogues(1000, 500)
#' sim_data <- create_mut_catalogues(1, 500, sig_set = c(1,4,29), c_exposure = c(0.25, 0.65, 0.1))
#'
#' @export
create_mut_catalogues <- function(n, m, P = cosmicSigs, sig_set = NULL,
                                  c_exposure = NULL){
  if (! is.matrix(P)){
    stop('P should be a matrix', call. = TRUE)
  }

  if (!is.numeric(n) || n < 1){
    stop('n has to be a positive number', call. = TRUE)
  }

  if (!is.numeric(m) || m < 1){
    stop('m has to be a positive number', call. = TRUE)
  }

  nsigs <- ncol(P)
  nfeatures <- nrow(P)

  if (is.null(sig_set)){
    sig_set  <- c(seq_len(nsigs))
  }
  else if (!is.numeric(sig_set) || length(sig_set) < 1){
    stop('sig_set has to be a vector containing numeric values', call. = TRUE)
  }
  else if (max(sig_set) > nsigs){
    stop('The signature matrix P has less columns then sig_set is accessing',
         call. = TRUE)
  }

  if (!is.null(c_exposure)) {
    if(!is.numeric(c_exposure)){
      stop('c_exposure has to be a vector containing numeric values',
           call. = TRUE)
    }
    else if(length(c_exposure) != length(sig_set)){
      stop('c_exposure has to be a vector of the same length as sig_set',
           call. = TRUE)
    }
    if(sum(c_exposure) != 1){
      warning('c_exposure does not add up to one, the weights will be scaled 
              accordingly')
      c_exposure = c_exposure/sum(c_exposure)
    }
    
  }

  #============================================================================
  #create catalogues
  catalogues = matrix(0, nrow = nfeatures, ncol = n)
  exposure = matrix(0, nrow = nsigs, ncol = n)

  for (i in seq_len(n)) {
    v = matrix(0,nrow = nsigs, ncol = 1)
    for (j in seq_len(length(sig_set))){
      if (is.null(c_exposure)){
        v[sig_set[j]] <- stats::runif(1, 0, 1)
      }
      else{
        v[sig_set[j]] <- c_exposure[j]
      }
    }
    v <- v/sum(v)
    p <- P %*% v
    p <- p/sum(p)
    exposure[,i] <- v

    #sample
    mutations_sampled <- sample(seq(nfeatures), m, replace = TRUE, prob = p)
    sample <- as.numeric(table(factor(mutations_sampled,
                                         levels = seq(nfeatures))))
    catalogues[,i] <- sample
  }

  sample_names <- paste0('sample_', seq_len(n))
  rownames(catalogues) <- rownames(P)
  rownames(exposure) <- colnames(P)
  colnames(catalogues) <- sample_names
  colnames(exposure) <- sample_names

  sim_data <- list(catalogues, exposure)
  names(sim_data) <- c('catalogues', 'exposure')

  return(sim_data)
}


#' Estimates the signature exposure of a mutational catalogue
#'
#' Estimates the signature exposure of a mutational catalogue by reconstructing
#' it from a chosen set of signatures, by default, the used method is quadratic
#' programming
#'
#' @param mut_cat matrix (f by n) conatining one or several mutational 
#' catalogues of samples whose signature exposures are to be estimated
#' f is the amount of features that the profiles are defined on
#' n is the number of catalogues
#' @param P Matrix (f by k) containing the signature profiles of k
#' signatures whose exposure is to be found
#' k is the amount of signature profiles that P contains
#' f is the amount of features that the profiles are defined on
#' (default: COSMIC signature matrix 96 x 30)
#' @param sig_set Numeric vector containing the index of the columns from the
#' signature matrix (which are the signature profiles) that will be used
#' reconstruct the mutational catalogue.
#' @param FUN Function to estimate the signature exposure. Default: Quadratic 
#' programming (P has to be of full rank to use this method)
#' @param ... control parameters that will be passed to FUN
#'
#' @return List of the estimated signature exposure, the reconstructed profile
#' of the sample, the cosine similarity between the two and the error
#'
#' @examples
#' data(cosmicSigs)
#' signature_exposure(create_mut_catalogues(10,500)[['catalogues']])
#' signature_exposure(create_mut_catalogues(10,500)[['catalogues']], sig_set = c(2,7,16,28,30))
#' signature_exposure(as.matrix(create_mut_catalogues(10,500)[['catalogues']][,1]))
#' 
#' 
#'
#' @export
signature_exposure <- function(mut_cat, P = cosmicSigs, sig_set = NULL, FUN=decomposeQP, ... ){
    
    if (! is.matrix(P)){
        stop('P should be a matrix', call. = TRUE)
    }
    if (ncol(P) == 1){
        stop("Matrix 'P' must have at least 2 columns.", call. = TRUE)
    }
    if (! is.numeric(mut_cat)){
        stop('The mutational catalogue is expected to contain only numbers.',
             call. = TRUE)
    }
    if (nrow(mut_cat) != nrow(P)){
        stop("The mutational catalogues (mut_cat) and the signatures (P) must
         have the same number of rows.", call. = TRUE)
    }
    nsigs <- ncol(P)
    nfeatures <- nrow(P)
    
    # in case the catalogue is a vector
    mut_cat <- as.matrix(mut_cat)
    sample_names <- paste0('sample_', seq_len(ncol(mut_cat)))
    colnames(mut_cat) <- sample_names
    
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
    
    P_adj <- P[,sig_set]
    
    
    # some parts are from 'findSigExposure' in SigEs
    mut_cat <- apply(mut_cat, 2, function(x) {
        x/sum(x)
    })
    estimate <- apply(mut_cat, 2, function(x) {
        FUN(x, P_adj, ...)
    })
    reconstructed <- t(t(estimate)%*%t(P_adj))
    estimate_adj = matrix(0L,  nrow = dim(P)[2], ncol = dim(mut_cat)[2])
    estimate_adj[sig_set,] <- estimate
    rownames(estimate_adj) <- colnames(P)
    colnames(estimate_adj) <- colnames(mut_cat)
    
    
    # error
    errors <- (mut_cat - reconstructed)**2
    errors <- as.matrix(colSums(errors))
    rownames(errors) <- colnames(mut_cat)
    colnames(errors) <- c('error')
    cossim <- numeric(length = ncol(mut_cat))
    for(i in seq_len(ncol(mut_cat))){
        cossim[i] <- as.numeric(mut_cat[,i] %*% reconstructed[,i]/(sqrt(mut_cat[,i] %*% mut_cat[,i]) * sqrt(reconstructed[,i] %*% reconstructed[,i])))
    }
    res <- list(estimate_adj, reconstructed, cossim, errors)
    names(res) <- c('exposures', 'reconstructed profiles', 'cosine similarity',
                    'errors')
    
    return(res)
}

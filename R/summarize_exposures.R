#' Signature exposure analysis of a mutational catalogue
#'
#' Bootstrapps a mutational catalogue and details the results about the 
#' signature estimation in a table and a boxplot
#'
#' @param mut_cat mutational catalogue (96 by 1)
#' @param P Matrix (f x k) containing the signature profiles of k
#' signatures whose exposure is to be found
#' k is the amount of signature profiles that P contains
#' f is the amount of features that the profiles are defined on
#' (default: COSMIC signature matrix 96 by 30)
#' @param m numeric, amount of mutations in the profile (needed for 
#' bootstrapping in case mut_cat doesn't contain absolute counts.)
#' @param plotting boolean, TRUE by default, if True,  a boxplot showing the 
#' distribution of estimated signature exposure for all the re-samples, 
#' highlighting the one of the original mutational catalogue, thus providing 
#' insights on the reliability of the estimates
#'
#' @return table (k by 6) containing statistics for each signature about the 
#' estimated exposure over 1000 bootstrapped resamples
#' - original exposure estimated exposure of the original mutational profile
#' - min minimum estimated exposure to this signature
#' - 1. quartile
#' - median 
#' - 3. quartile
#' - max maximum estimated exposure to this signature
#'
#' if the plotting option is chosen, a boxplot will also be displayed
#' 
#' @examples
#' summarize_exposures(create_mut_catalogues(10,500)[['catalogues']][,1], plotting=FALSE)
#' 
#' @export
summarize_exposures <- function(mut_cat, P = cosmicSigs,  plotting=TRUE, m=NULL){

  if (! is.numeric(mut_cat)){
    stop('The mutational catalogue is expected to contain only numbers.',
         call. = TRUE)
  }
  mut_cat <- as.matrix(mut_cat)

  QPoriginal <- signature_exposure(mut_cat, P)
  
  if(is.null(m)){
    m <- colSums(mut_cat)
  }
  bootstrapped <- bootstrap_mut_catalogues(1000, mut_cat[,1], m)
  
  QPbootstrapped <- signature_exposure(bootstrapped, P)
  
  exposure_stats <- matrix(nrow = ncol(P), ncol = 6)
  rownames(exposure_stats) = colnames(P)
  colnames(exposure_stats) = c('original exposure', 'min','1. quartile',
                               'median', '3. quartile', 'max')
  original_exposures <- QPoriginal[['exposures']][,1]
  boot_exposures <- QPbootstrapped[['exposures']]
  
  exposure_stats[,1] <- original_exposures
  exposure_stats[,2] <- Biobase::rowMin(boot_exposures)
  exposure_stats[,3] <- apply(boot_exposures, 1, function(x){
    stats::quantile(x, 0.25)
  })
  exposure_stats[,4] <- Biobase::rowMedians(boot_exposures)
  exposure_stats[,5] <- apply(boot_exposures, 1, function(x){
    stats::quantile(x, 0.75)
  })
  exposure_stats[,6] <- Biobase::rowMax(boot_exposures)
  
  #==========================================================================
  #Plot
  if(plotting){
    plot_bootstrapped_exposure(boot_exposures, as.matrix(QPoriginal[['exposures']]))
  }
  #==========================================================================
  
  return(exposure_stats)
}

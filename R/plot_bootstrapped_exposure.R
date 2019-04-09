#' Plot signature exposure estimation for several samples
#'
#' Creates a boxplot illustrating the results of the signature estimation for 
#' several mutational catalogues (e.g. bootstrapped re-samples or a cohort). 
#' The plot shows the distribution of estimated signature exposure for all the 
#' catalogues, highlighting the one of the original mutational catalogue if one
#' is provided.
#'
#' @param bootstrapped_exposure matrix (n by k) containing the signature 
#' exposure of several mutational catalogues (bootstrapped re-samples)
#' n is the amount of re-samples
#' k is the amount of signature profiles that P contains
#' @param original_estimation matrix (n by k) containing the signature exposure
#' of one or several mutational catalogues (typically, the original sample) that
#' will be diplayed on top of the boxplots
#' @param title character, title of the plot
#' @param box_col color option of the boxplot
#' @param point_col color option of the points that indicate the exposure of the 
#' original profile(s). Should be a vector, if several original profiles are 
#' plotted
#' @param sig_names character vector, names of the signatures to be displayed as
#' x-axis labels
#' @param sample_names character vector, names of the original profile(s)
#'
#' @return
#' Diplays a boxplot
#' 
#' @note The function can of course also be used to plot the distribution 
#' of estimated signature exposures in a cohort instead of one bootstrapped 
#' sample.
#' 
#' @examples
#' # prepare input
#' data(cosmicSigs)
#' mut_cat <- create_mut_catalogues(4,400)
#' exposures <- signature_exposure(bootstrap_mut_catalogues(
#' 1000, mut_cat[['catalogues']][,1]))[['exposures']]
#' original_exposure <- signature_exposure(mut_cat[['catalogues']])[['exposures']]
#' 
#' plot_bootstrapped_exposure(exposures, as.matrix(original_exposure[,1]))
#' 
#' plot_bootstrapped_exposure(exposures, as.matrix(original_exposure), 
#' title='Example', box_col='grey')
#' 
#' @export
plot_bootstrapped_exposure <- function(bootstrapped_exposure,
                                      original_estimation=NULL, 
                                      title=NULL, 
                                      box_col=NULL, 
                                      point_col=NULL, 
                                      sig_names=NULL,
                                      sample_names=NULL){

    if(!is.numeric(bootstrapped_exposure)){
      stop('The mutational catalogues should contain only numbers', call. = TRUE)
    }
    if(!is.null(original_estimation)){
        if(nrow(bootstrapped_exposure) != nrow(original_estimation)){
            stop('The bootstrapped exposures and the original estimation have to have the same amount of rows (signatures)', call. = TRUE)
        }
        else if(!is.numeric(original_estimation)){
          stop('The mutational catalogues should contain only numbers', call. = TRUE)
        }
      
        if(is.null(sample_names)){
          if (ncol(original_estimation) == 1){
            sample_names = 'original profile'
        }
        else{
          sample_names = as.character(seq_len(ncol(original_estimation)))
        }
      }
    }
    else{
        original_estimation <- c()
    }
    if(is.null(title)){
        title <- paste0('Bootstrap distribution ', sample_names[1])
    }
    else{
        if(!is.character(title)){
            stop('Title has to be a character vector', call. = TRUE)
        }
    }
    if(is.null(box_col)){
        box_col <- "aquamarine3"
    }
    if(is.null(point_col)){
        point_col <- c("darkorange3", "darkmagenta", "red", "blue", "black")
        
    }
    if(is.null(sig_names)){
        sig_names <- seq_len(nrow(bootstrapped_exposure))
    }
    else{
        if(length(sig_names)!=nrow(bootstrapped_exposure)){
            stop('Amount of sig_names has to match the amount of signatures', call.=TRUE)
        }
    }

  #============================================================================
  # Plot
  
    boxplot(t(bootstrapped_exposure), range = 0,
            ylim = c(0, (max(bootstrapped_exposure)+0.10)), frame=F, 
            ylab = "Exposure", main = title, cex.axis = 1, xaxt = "n", 
            xlab = "Signatures", col = box_col)
    axis(1, seq_len(nrow(bootstrapped_exposure)), 
         labels = as.character(sig_names), las = 1, tick=F, cex.axis = 0.7, 
         mgp = c(2,0.01,0))
    if(!is.null(original_estimation)){
        for(i in seq_len(ncol(original_estimation))){
            points(original_estimation[,i], pch = '*', col = point_col[i], 
                   cex = 2, lwd = 2)
        }
        legend('topleft', legend=sample_names,
               col=point_col, pch = '*', cex=1.2, box.lty=0)
    }
}
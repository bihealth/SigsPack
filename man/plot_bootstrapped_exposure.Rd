% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_bootstrapped_exposure.R
\name{plot_bootstrapped_exposure}
\alias{plot_bootstrapped_exposure}
\title{Plot signature exposure estimation for several samples}
\usage{
plot_bootstrapped_exposure(bootstrapped_exposure,
  original_estimation = NULL, title = NULL, box_col = NULL,
  point_col = NULL, sig_names = NULL, sample_names = NULL)
}
\arguments{
\item{bootstrapped_exposure}{matrix (n by k) containing the signature 
exposure of several mutational catalogues (bootstrapped re-samples)
n is the amount of re-samples
k is the amount of signature profiles that P contains}

\item{original_estimation}{matrix (n by k) containing the signature exposure
of one or several mutational catalogues (typically, the original sample) that
will be diplayed on top of the boxplots}

\item{title}{character, title of the plot}

\item{box_col}{color option of the boxplot}

\item{point_col}{color option of the points that indicate the exposure of the 
original profile(s). Should be a vector, if several original profiles are 
plotted}

\item{sig_names}{character vector, names of the signatures to be displayed as
x-axis labels}

\item{sample_names}{character vector, names of the original profile(s)}
}
\value{
Diplays a boxplot
}
\description{
Creates a boxplot illustrating the results of the signature estimation for 
several mutational catalogues (e.g. bootstrapped re-samples or a cohort). 
The plot shows the distribution of estimated signature exposure for all the 
catalogues, highlighting the one of the original mutational catalogue if one
is provided.
}
\note{
The function can of course also be used to plot the distribution 
of estimated signature exposures in a cohort instead of one bootstrapped 
sample.
}
\examples{
# prepare input
data(cosmicSigs)
mut_cat <- create_mut_catalogues(4,400)
exposures <- signature_exposure(bootstrap_mut_catalogues(
1000, mut_cat[['catalogues']][,1]))[['exposures']]
original_exposure <- signature_exposure(mut_cat[['catalogues']])[['exposures']]

plot_bootstrapped_exposure(exposures, as.matrix(original_exposure[,1]))

plot_bootstrapped_exposure(exposures, as.matrix(original_exposure), 
title='Example', box_col='grey')

}

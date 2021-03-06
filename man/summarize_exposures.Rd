% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_exposures.R
\name{summarize_exposures}
\alias{summarize_exposures}
\title{Signature exposure analysis of a mutational catalogue}
\usage{
summarize_exposures(mut_cat, P = get(utils::data("cosmicSigs", package =
  "SigsPack")), plotting = TRUE, m = NULL)
}
\arguments{
\item{mut_cat}{mutational catalogue (96 by 1)}

\item{P}{Matrix (f x k) containing the signature profiles of k
signatures whose exposure is to be found
k is the amount of signature profiles that P contains
f is the amount of features that the profiles are defined on
(default: COSMIC signature matrix 96 by 30)}

\item{plotting}{boolean, TRUE by default, if True,  a boxplot showing the 
distribution of estimated signature exposure for all the re-samples, 
highlighting the one of the original mutational catalogue, thus providing 
insights on the reliability of the estimates}

\item{m}{numeric, amount of mutations in the profile (needed for 
bootstrapping in case mut_cat doesn't contain absolute counts.)}
}
\value{
table (k by 6) containing statistics for each signature about the 
estimated exposure over 1000 bootstrapped resamples
- original exposure estimated exposure of the original mutational profile
- min minimum estimated exposure to this signature
- 1. quartile
- median 
- 3. quartile
- max maximum estimated exposure to this signature

if the plotting option is chosen, a boxplot will also be displayed
}
\description{
Bootstrapps a mutational catalogue and details the results about the 
signature estimation in a table and a boxplot
}
\examples{
summarize_exposures(create_mut_catalogues(10,500)[['catalogues']][,1], plotting=FALSE)

}

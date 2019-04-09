#' decomposeQP function
#'
#' This function is taken from the package 'SignatureEstimation' by Xiaoqing 
#' Huang and Damian Wojtowicz 
#' (source: https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#signatureestimation)
#' The function allows to get the optimal solution by using dual method to solve
#' the quadratic programming problem.
#'
#'
#' @param m observed tumor profile vector for a single patient/sample, 96 by 1. 
#' m is normalized.
#' @param P signature profile matrix, 96 by N (N = # signatures, COSMIC: N=30)
#' @param ... control parameter that can be passed into the solve.QP 
#' function
#'
#' @return matrix containing estimated signature exposures
#'
#' @examples
#' data(cosmicSigs)
#' mut_cat<- (create_mut_catalogues(1,500)[['catalogues']])/500
#' decomposeQP(mut_cat, cosmicSigs)
#' @export
decomposeQP <- function (m, P, ...)
{
    N = ncol(P)
    G = t(P) %*% P
    C <- cbind(rep(1, N), diag(N))
    b <- c(1, rep(0, N))
    d <- t(m) %*% P
    out = quadprog::solve.QP(Dmat = G, dvec = d, Amat = C, bvec = b,
                             meq = 1)
    exposures = out$solution
    exposures[exposures < 0] = 0
    exposures = exposures/sum(exposures)
    return(exposures)
}

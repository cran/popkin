#' Extract FST from a population-level kinship matrix or vector of inbreeding coefficients
#'
#' This function simply returns the weighted mean inbreeding coefficient \eqn{f_j^T}.
#' If no weights are provided, the regular mean \eqn{f_j^T} is returned.
#' If a kinship matrix \eqn{\Phi^T} is provided, then \eqn{f_j^T} are extracted from its diagonal using \code{\link{inbr}} (assumes the diagonal of \eqn{\Phi^T} is \eqn{\phi_{jj}^T = \frac{1}{2}(1+f_j^T)}{\phi_jj^T = (1+f_j^T)/2} as \code{\link{popkin}} returns, and not \eqn{f_j^T} as \code{\link{inbrDiag}} returns).
#'
#' The returned weighted mean inbreeding coefficient equals the generalized \eqn{F_{ST}}{FST} if all individuals are "locally outbred" (i.e. if the self relatedness of every individual stems entirely from the population structure rather than due partly to having unusually closely related parents, such as first or second cousins).
#' Note most individuals in population-scale human data are locally outbred.
#' If there are locally inbred individuals, the returned value will overestimate \eqn{F_{ST}}{FST}.
#' 
#' @param x The vector of inbreeding coefficients \eqn{(f_j^T)}, or the kinship matrix \eqn{\Phi^T} if \code{x} is a matrix.
#' @param w Weights for individuals (optional, defaults to uniform weights)
#'
#' @return \eqn{F_{ST}}{FST}
#'
#' @examples
#' ## Get FST from a genotype matrix
#' 
#' ## Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' 
#' ## NOTE: for BED-formatted input, use BEDMatrix!
#' ## "file" is path to BED file (excluding .bed extension)
#' # library(BEDMatrix)
#' # X <- BEDMatrix(file) # load genotype matrix object
#'
#' ## estimate the kinship matrix "Phi" from the genotypes "X"!
#' Phi <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#' w <- weightsSubpops(subpops) # can weigh individuals so subpopulations are balanced
#' Fst <- fst(Phi, w) # use kinship matrix and weights to calculate fst
#' 
#' Fst <- fst(Phi) # no weights implies uniform weights
#'
#' inbr <- inbr(Phi) # if you extracted inbr for some other analysis...
#' Fst <- fst(inbr, w) # ...use this inbreeding vector as input too!
#'
#' @export
fst <- function(x, w) {
    ## if input is a matrix, let's assume it is the kinship matrix, so extract the inbreeding coefficients first
    if (class(x) == 'matrix') {
        x <- inbr(x)
    }
    ## now x is a vector of inbreeding coefficients
    if (missing(w)) {
        return( mean(x) ) # no weights implies uniform weights
    } else {
        return( drop( x %*% w ) ) # weighted mean
    }
}

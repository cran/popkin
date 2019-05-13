#' Rescale kinship matrix to set a given kinship value to zero.
#'
#' Rescales the input kinship matrix \eqn{\Phi^T} so that the value \eqn{\phi_{\mbox{min}}^T}{\phi_min^T} in the original kinship matrix becomes zero, using the formula
#' \deqn{\Phi^{T'} = \frac{\Phi^T - \phi_{\mbox{min}}^T}{1 - \phi_{\mbox{min}}^T}.}{\Phi^T' = (\Phi^T - \phi_min^T)/(1 - \phi_min^T).}
#' This is equivalent to changing the ancestral population \eqn{T} into \eqn{T'} such that \eqn{\phi_{\mbox{min}}^{T'} = 0}{\phi_min^T' = 0}.
#' If subpopulation labels \code{subpops} are provided, they are used to estimate \eqn{\phi_{\mbox{min}}^T}{\phi_min^T}.
#' If both \code{subpops} and \code{min_kinship} are provided, only \code{min_kinship} is used.
#' If both \code{subpops} and \code{min_kinship} are omitted, the adjustment is equivalent to \code{min_kinship=min(kinship)}.
#'
#' @param kinship An \eqn{n \times n}{n-by-n} kinship matrix.
#' @param subpops The length-\eqn{n} vector of subpopulation assignments for each individual.
#' @param min_kinship A scalar kinship value to define the new zero kinship.
#'
#' @return The rescaled \eqn{n \times n}{n-by-n} kinship matrix, with the desired level of relatedness set to zero.
#'
#' @examples
#' # Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' subpops2 <- 1:3 # alternate labels treat every individual as a different subpop
#' 
#' # NOTE: for BED-formatted input, use BEDMatrix!
#' # "file" is path to BED file (excluding .bed extension)
#' ## library(BEDMatrix)
#' ## X <- BEDMatrix(file) # load genotype matrix object
#'
#' # suppose we first estimate kinship without subpopulations, which will be more biased
#' kinship <- popkin(X) # calculate kinship from genotypes, WITHOUT subpops
#' # then we visualize this matrix, figure out a reasonable subpopulation partition
#'
#' # now we can adjust the kinship matrix!
#' kinship2 <- rescale_popkin(kinship, subpops)
#' # prev is faster but otherwise equivalent to re-estimating kinship from scratch with subpops:
#' # kinship2 <- popkin(X, subpops) 
#'
#' # can also manually set the level of relatedness min_kinship we want to be zero:
#' min_kinship <- min(kinship) # a naive choice for example
#' kinship2 <- rescale_popkin(kinship, min_kinship = min_kinship)
#'
#' # lastly, omiting both subpops and min_kinship sets the minimum value in kinship to zero
#' kinship3 <- rescale_popkin(kinship2)
#' # equivalent to both of:
#' # kinship3 <- popkin(X)
#' # kinship3 <- rescale_popkin(kinship2, min_kinship = min(kinship))
#'
#' @export
rescale_popkin <- function(kinship, subpops = NULL, min_kinship = NA) {
    # die if this is missing
    if (missing(kinship))
        stop('`kinship` matrix is required!')
    
    # run additional validations
    validate_kinship(kinship)

    if (is.na(min_kinship))
        min_kinship <- min_mean_subpops(kinship, subpops)
    
    # finally, perform a simple IBD rescaling
    kinship <- (kinship - min_kinship)/(1 - min_kinship) # return this matrix!
}

# stick deprecated function name here

#' @title Rescale kinship matrix to set a given kinship value to zero.
#' @description Rescale kinship matrix to set a given kinship value to zero.
#' @param kinship A kinship matrix
#' @param subpops Vector of subpopulation assignments for each individual.
#' @param min_kinship A scalar kinship value to define the new zero kinship.
#' @return The rescaled kinship matrix, with the desired level of relatedness set to zero.
#'
#' @name rescalePopkin-deprecated
#' @usage rescalePopkin(kinship, subpops = NULL, min_kinship = NA)
#' @seealso \code{\link{popkin-deprecated}}
#' @keywords internal
NULL

#' @rdname popkin-deprecated
#' @section \code{rescalePopkin}:
#' For \code{rescalePopkin}, use \code{\link{rescale_popkin}}.
#'
#' @export
rescalePopkin <- function(kinship, subpops = NULL, min_kinship = NA) {
    # mark as deprecated
    .Deprecated('rescale_popkin')
    # return as usual, to not break things just yet
    rescale_popkin(kinship, subpops, min_kinship)
}
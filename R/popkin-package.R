#' A package for estimating kinship and FST under arbitrary population structure
#'
#' The heart of this package is the \code{\link{popkin}} function, which estimates the kinship matrix of all individual pairs from their genotype matrix.
#' Inbreeding coefficients, the generalized \eqn{F_{ST}}{FST}, and the individual-level pairwise \eqn{F_{ST}}{FST} matrix are extracted from the kinship matrix using \code{\link{inbr}}, \code{\link{fst}}, and \code{\link{pwfst}}, respectively.
#' \code{\link{fst}} accepts weights for individuals to balance subpopulations obtained with \code{\link{weightsSubpops}}.
#' Kinship matrices can be renormalized (to change the most recent common ancestor population or MRCA) using \code{\link{rescalePopkin}}.
#' Lastly, kinship and pairwise \eqn{F_{ST}}{FST} matrices can be visualized using \code{\link{plotPopkin}} (with the help of \code{\link{inbrDiag}} for kinship matrices only).
#' 
#' @examples
#' ## estimate and visualize kinship and FST from a genotype matrix
#'
#' ## Construct toy data
#' X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow=3, byrow=TRUE) # genotype matrix
#' subpops <- c(1,1,2) # subpopulation assignments for individuals
#' subpops2 <- 1:3 # alternate labels treat every individual as a different subpop
#' 
#' ## NOTE: for BED-formatted input, use BEDMatrix!
#' ## "file" is path to BED file (excluding .bed extension)
#' # library(BEDMatrix)
#' # X <- BEDMatrix(file) # load genotype matrix object
#'
#' ## estimate the kinship matrix "Phi" from the genotypes "X"!
#' ## all downstream analysis require "Phi", none use "X" after this
#' Phi <- popkin(X, subpops) # calculate kinship from X and optional subpop labels
#'
#' ## plot the kinship matrix, marking the subpopulations
#' ## note inbrDiag replaces the diagonal of Phi with inbreeding coefficients
#' plotPopkin( inbrDiag(Phi), labs=subpops )
#'
#' ## extract inbreeding coefficients from Phi
#' inbr <- inbr(Phi)
#' 
#' ## estimate FST
#' w <- weightsSubpops(subpops) # weigh individuals so subpopulations are balanced
#' Fst <- fst(Phi, w) # use kinship matrix and weights to calculate fst
#' Fst <- fst(inbr, w) # estimate more directly from inbreeding vector (same result)
#'
#' ## estimate and visualize the pairwise FST matrix
#' pwF <- pwfst(Phi) # estimated matrix
#' legTitle <- expression(paste('Pairwise ', F[ST])) # fancy legend label
#' plotPopkin(pwF, labs=subpops, legTitle=legTitle) # NOTE no need for inbrDiag() here!
#'
#' ## rescale the kinship matrix using different subpopulations (implicitly changes the MRCA)
#' Phi2 <- rescalePopkin(Phi, subpops2)
#'
#' @docType package
#' @name popkin-package
#' @aliases popkin-package
"_PACKAGE"

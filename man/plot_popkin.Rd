% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_popkin.R
\name{plot_popkin}
\alias{plot_popkin}
\title{Visualize one or more kinship matrices}
\usage{
plot_popkin(kinship, titles = NULL, col = NULL, col_n = 100,
  mar = NULL, mar_pad = 0.2, diag_line = FALSE,
  panel_letters = toupper(letters), panel_letters_cex = 1.5,
  ylab = "Individuals", ylab_adj = NA, ylab_line = 0,
  layout_add = TRUE, layout_rows = 1, leg_title = "Kinship",
  leg_cex = 1, leg_mar = NULL, leg_n = 5, names = FALSE,
  names_cex = 1, names_line = NA, labs = NULL, labs_cex = 1,
  labs_las = 0, labs_line = 0, labs_sep = TRUE, labs_lwd = 1,
  labs_col = "black", labs_ticks = FALSE, labs_text = TRUE,
  labs_even = FALSE, ...)
}
\arguments{
\item{kinship}{A numeric kinship matrix or a list of matrices.
Note \code{kinship} may contain \code{NULL} elements (makes blank plots with titles; good for placeholders or non-existent data)}

\item{titles}{Titles to add to each matrix panel (default is no titles)}

\item{col}{Colors for heatmap (default is a red-white-blue palette symmetric about zero constructed using RColorBrewer).}

\item{col_n}{The number of colors to use in the heatmap (applies if \code{col = NULL}).}

\item{mar}{Margins for each panel (if a list) or for all panels (if a vector).
Margins are in \code{c(bottom,left,top,right)} format that \code{\link[graphics]{par}('mar')} expects.
Note the padding \code{mar_pad} below is also added to every margin if set.
If \code{NULL}, the original margin values are used without change, and are reset for every panel that has a \code{NULL} value.
The original margins are also reset after plotting is complete.}

\item{mar_pad}{Margin padding added to all panels (\code{mar} above and \code{leg_mar} below).
Default 0.2.
Must be a scalar or a vector of length 4 to match \code{\link[graphics]{par}('mar')}.}

\item{diag_line}{If \code{TRUE} adds a line along the diagonal (default no line).
May also be a vector of booleans to set per panel (lengths must agree).}

\item{panel_letters}{Vector of strings for labeling panels (default A-Z).
No labels are added if \code{panel_letters = NULL}, or when there is only one panel except if \code{panel_letters} is set to a single letter in that case (this behavior is useful if goal is to have multiple external panels but popkin only creates one of these panels).}

\item{panel_letters_cex}{Scaling factor of panel letters (default 1.5).

AXIS LABEL OPTIONS}

\item{ylab}{The y-axis label (default "Individuals").
If \code{length(ylab) == 1}, the label is placed in the outer margin (shared across panels);
otherwise \code{length(ylab)} must equal the number of panels and each label is placed in the inner margin of the respective panel.}

\item{ylab_adj}{The value of "adj" passed to \code{\link[graphics]{mtext}}.
If \code{length(ylab) == 1}, only the first value is used, otherwise \code{length(ylab_adj)} must equal the number of panels.}

\item{ylab_line}{The value of "line" passed to \code{\link[graphics]{mtext}}.
If \code{length(ylab) == 1}, only the first value is used, otherwise \code{length(ylab_line)} must equal the number of panels.

LAYOUT OPTIONS}

\item{layout_add}{If \code{TRUE} (default) then \code{\link[graphics]{layout}} is called internally with appropriate values for the required number of panels for each matrix, the desired number of rows (see \code{layout_rows} below) plus the color key legend.
The original layout is reset when plotting is complete and if \code{layout_add = TRUE}.
If a non-standard layout or additional panels (beyond those provided by \code{plot_popkin}) are desired, set to FALSE and call \code{\link[graphics]{layout}} yourself beforehand.}

\item{layout_rows}{Number of rows in layout, used only if \code{layout_add = TRUE}.

LEGEND (COLOR KEY) OPTIONS}

\item{leg_title}{The name of the variable that the heatmap colors measure (default "Kinship").}

\item{leg_cex}{Scaling factor for \code{leg_title} (default 1).}

\item{leg_mar}{Margin vector (in \code{c(bottom,left,top,right)} format that \code{\link[graphics]{par}('mar')} expects) for the legend panel only.
If not provided, the margins used in the last panel are preserved with the exception that the left margin is set to zero (plus the value of \code{mar_pad}, see above).}

\item{leg_n}{The desired number of ticks in the legend y-axis (input to \code{\link{pretty}}, see that for more details).

INDIVIDUAL LABEL OPTIONS}

\item{names}{If \code{TRUE}, the column and row names are plotted in the heatmap.}

\item{names_cex}{Scaling factor for the column and row names.}

\item{names_line}{Line where column and row names are placed.

SUBPOPULATION LABEL OPTIONS}

\item{labs}{Subpopulation labels for individuals.
Use a matrix of labels to show groupings at more than one level (for a hierarchy or otherwise).
If input is a vector or a matrix, the same subpopulation labels are shown for every heatmap panel; the input must be a list of such vectors or matrices if the labels vary per panel.}

\item{labs_cex}{A vector of label scaling factors for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_las}{A vector of label orientations (in format that \code{\link[graphics]{mtext}} expects) for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_line}{A vector of lines where labels are placed (in format that \code{\link[graphics]{mtext}} expects) for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_sep}{A vector of booleans that specify whether lines separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_lwd}{A vector of line widths for the lines that divide subpopulations (if \code{labs_sep = TRUE}) for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_col}{A vector of colors for the lines that divide subpopulations (if \code{labs_sep = TRUE}) for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_ticks}{A vector of booleans that specify whether ticks separating the subpopulations are drawn for each level of labs, or a list of such vectors if labels vary per panel.}

\item{labs_text}{A vector of booleans that specify whether the subpopulation labels are shown for each level of labs, or a list of such vectors if labels vary per panel.
Useful for including separating lines or ticks without text.}

\item{labs_even}{A vector of booleans that specify whether the subpopulations labels are drawn with equal spacing for each level of labs, or a list of such vectors if labels vary per panel.
When \code{TRUE}, lines mapping the equally-spaced labels to the unequally-spaced subsections of the heatmap are also drawn.}

\item{...}{Additional options passed to \code{\link[graphics]{image}}.
These are shared across panels}
}
\description{
This function plots one or more kinship matrices and a shared legend for the color key.
Many options allow for fine control of individual or subpopulation labeling.
This code assumes input matrices are symmetric.
}
\details{
\code{plot_popkin} plots the input kinship matrices as-is.
For best results, a standard kinship matrix (such as the output of \code{\link{popkin}}) should have its diagonal rescaled to contain inbreeding coefficients (\code{\link{inbr_diag}} does this) before \code{plot_popkin} is used.

This function permits the labeling of individuals (from row and column names when \code{names = TRUE}) and of subpopulations (passed through \code{labs}).
The difference is that the labels passed through \code{labs} are assumed to be shared by many individuals, and lines (or other optional visual aids) are added to demarcate these subgroups.
}
\examples{
# Construct toy data
X <- matrix(c(0,1,2,1,0,1,1,0,2), nrow = 3, byrow = TRUE) # genotype matrix
subpops <- c(1,1,2) # subpopulation assignments for individuals

# NOTE: for BED-formatted input, use BEDMatrix!
# "file" is path to BED file (excluding .bed extension)
## library(BEDMatrix)
## X <- BEDMatrix(file) # load genotype matrix object

# estimate the kinship matrix from the genotypes "X"!
kinship <- popkin(X, subpops) # calculate kinship from X and optional subpop labels

# simple plot of the kinship matrix, marking the subpopulations only
# note inbr_diag replaces the diagonal of kinship with inbreeding coefficients
# (see vignette for more elaborate examples)
plot_popkin( inbr_diag(kinship), labs = subpops )

}
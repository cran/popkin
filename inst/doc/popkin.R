## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- cache = FALSE, include = FALSE-------------------------------------
## copied from examples from the "simmer" R package
## after: https://www.enchufa2.es/archives/suggests-and-vignettes.html
## by Iñaki Úcar
required <- c("lfa") # not a CRAN package, only suggested since popkin doesn't need it to run...

if (!all(sapply(required, requireNamespace, quietly = TRUE)))
  knitr::opts_chunk$set(eval = FALSE)

## ------------------------------------------------------------------------
library(popkin)
library(lfa) # for hgdp_subset sample data only
X <- hgdp_subset # rename for simplicity
dim(X)

## ------------------------------------------------------------------------
# shorten subpopulation labels
colnames(X)[colnames(X) == 'AFRICA'] <- 'AFR'
colnames(X)[colnames(X) == 'MIDDLE_EAST'] <- 'MDE'
colnames(X)[colnames(X) == 'EUROPE'] <- 'EUR'
colnames(X)[colnames(X) == 'CENTRAL_SOUTH_ASIA'] <- 'SAS'
colnames(X)[colnames(X) == 'EAST_ASIA'] <- 'EAS'
colnames(X)[colnames(X) == 'OCEANIA'] <- 'OCE'
colnames(X)[colnames(X) == 'AMERICA'] <- 'AMR'
# order roughly by distance from Africa
# (without crossing the Atlantic Ocean)
pop_order <- c('AFR', 'MDE', 'EUR', 'SAS', 'EAS', 'OCE', 'AMR')
# applies reordering
X <- X[, order(match(colnames(X), pop_order))]
subpops <- colnames(X) # extract subpopulations vector

## ------------------------------------------------------------------------
X[1:10,1:15]

## ------------------------------------------------------------------------
kinship <- popkin(X, subpops)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'--------------
# set outer margin for axis labels (left and right are non-zero)
par(oma = c(0,1.5,0,3))
# set inner margin for subpopulation labels (bottom and left are non-zero), add padding
par(mar = c(1,1,0,0)+0.2)
# now plot!
plot_popkin(kinship, labs = subpops)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'--------------
par(oma = c(0,1.5,0,3))
par(mar = c(1,1,0,0)+0.2)
plot_popkin(inbr_diag(kinship), labs = subpops)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'--------------
par(oma = c(0,1.5,0,3))
# increase margins because labels go farther out
par(mar = c(2,2,0,0)+0.2)
plot_popkin(inbr_diag(kinship), labs = subpops, labs_even = TRUE, labs_line = 1, labs_cex = 0.7)

## ------------------------------------------------------------------------
# get weights
w <- weights_subpops(subpops)
# compute FST!
# Note: don't use the output to inbr_diag(kinship) or FST will be wrong!
fst(kinship, w)

## ---- fig.width = 4, fig.height = 2, fig.align = 'center'----------------
inbrs <- inbr(kinship) # vector of inbreeding coefficients
# quick plot
par(mar = c(4, 4, 0, 0.2) + 0.2) # adjust margins
plot(density(inbrs), xlab = 'inbreeding coefficient', main = '') # see their distribution

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'--------------
pairwise_fst <- pwfst(kinship) # compute pairwise FST matrix from kinship matrix
leg_title <- expression(paste('Pairwise ', F[ST])) # fancy legend label
par(oma = c(0,1.5,0,3))
par(mar = c(2,2,0.2,0)+0.2)
# NOTE no need for inbr_diag() here!
plot_popkin(pairwise_fst, labs = subpops, labs_even = TRUE, labs_line = 1, labs_cex = 0.7, leg_title = leg_title)

## ---- fig.width = 3, fig.height = 2, fig.align = 'center'----------------
indexes_AFR <- subpops == 'AFR'
# AFR subset of the kinship matrix
kinship_AFR <- kinship[indexes_AFR, indexes_AFR]

# kinship matrix plot
par(oma = c(0,1.5,0,3))
par(mar = c(0,0,0,0)+0.2) # zero margins for no labels
plot_popkin(inbr_diag(kinship_AFR))

# estimate FST before rescaling (this value will be wrong, too high!)
fst(kinship_AFR)

## ---- fig.width = 3, fig.height = 2, fig.align = 'center'----------------
# rescale kinship_AFR
# since subpops is missing, minimum kinship value is set to zero
# (no averaging between subpopulations)
kinship_AFR <- rescale_popkin(kinship_AFR)

# kinship matrix plot
par(oma = c(0,1.5,0,3))
par(mar = c(0,0,0,0)+0.2) # zero margins for no labels
plot_popkin(inbr_diag(kinship_AFR))

# FST is now correct, relative to the MRCA of AFR individuals
fst(kinship_AFR)

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'--------------
par(oma = c(0,1.5,0,3))
# increase top margin for titles
par(mar = c(2,2,2,0)+0.2)
# dummy labels to have lines in second panel
subpops_AFR <- subpops[indexes_AFR]
plot_popkin(
    inbr_diag(list(kinship, kinship_AFR)), # inbr_diag also works on a list of matrices
    titles = c('All', 'AFR only, rescaled'), # title of each panel
    labs = list(subpops, subpops_AFR), # pass per-panel labels using a list
    labs_even = TRUE, # scalar options are shared across panels
    labs_line = 1,
    labs_cex = 0.5
)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'--------------
# create second level of labels
# first copy first-level labels
blocks <- subpops
# first block is AFR
blocks[blocks == 'AFR'] <- 'B1'
# second block is West Eurasians, broadly defined
blocks[blocks == 'MDE'] <- 'B2'
blocks[blocks == 'EUR'] <- 'B2'
blocks[blocks == 'SAS'] <- 'B2'
# third block is East Eurasians, broadly defined
blocks[blocks == 'EAS'] <- 'B3'
blocks[blocks == 'OCE'] <- 'B3'
blocks[blocks == 'AMR'] <- 'B3'

par(oma = c(0,1.5,0,3))
# increase margins again
par(mar = c(3,3,0,0)+0.2)
# plotting with different options per level is more complicated...
plot_popkin(
	inbr_diag(kinship),
	labs = cbind(subpops,blocks),   # ... labs is now a matrix with levels on columns
	labs_even = c(TRUE, FALSE),     # ... even spacing for first level only
	labs_line = c(1,2),             # ... put second level further out
	labs_cex = c(0.7, 1),           # ... don't shrink second level
	labs_sep = c(FALSE, TRUE),      # ... draw lines inside heatmap for second level only
	ylab_adj = 0.65                 # push up outer margin ylab "Individuals"
	)

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'--------------
par(oma = c(0,1.5,0,3))
par(mar = c(3,3,2,0)+0.2)
plot_popkin(
    inbr_diag(list(kinship, kinship_AFR)),
	titles = c('All', 'AFR only, rescaled'),
	labs = list(cbind(subpops, blocks), subpops_AFR), # list of matrices
	labs_even = c(TRUE, FALSE), # non-list: values are reused for both panels
	labs_line = c(1,2),
	# make label bigger in second panel (custom per-panel values)
	labs_cex = list(c(0.5, 0.7), 1), # list of vectors
	# add lines for first level of second panel (custom per-panel values)
	labs_sep = list(c(FALSE, TRUE), TRUE) # list of vectors
	)


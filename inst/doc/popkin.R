## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- cache = FALSE, include = FALSE------------------------------------------
## copied from examples from the "simmer" R package
## after: https://www.enchufa2.es/archives/suggests-and-vignettes.html
## by Iñaki Úcar
required <- c("lfa") # not a CRAN package, only suggested since popkin doesn't need it to run...

if (!all(sapply(required, requireNamespace, quietly = TRUE)))
  knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
library(popkin)
# for hgdp_subset sample data only
library(lfa)
# rename for simplicity
X <- hgdp_subset
dim(X)

## -----------------------------------------------------------------------------
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
# extract subpopulations vector
subpops <- colnames(X)

## -----------------------------------------------------------------------------
X[1:10,1:15]

## -----------------------------------------------------------------------------
kinship <- popkin(X, subpops)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
plot_popkin(
    kinship,
    labs = subpops,
    # shared bottom and left margin value, to make space for labels
    mar = 1
)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
plot_popkin(
    inbr_diag( kinship ),
    labs = subpops,
    mar = 1
)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
plot_popkin(
    inbr_diag(kinship),
    labs = subpops,
    labs_even = TRUE,
    labs_line = 1,
    labs_cex = 0.7,
    mar = 2
)

## -----------------------------------------------------------------------------
# get weights
w <- weights_subpops(subpops)
# compute FST!
# Note: don't use the output to inbr_diag(kinship) or FST will be wrong!
fst(kinship, w)

## ---- fig.width = 4, fig.height = 2, fig.align = 'center'---------------------
# vector of inbreeding coefficients
inbrs <- inbr(kinship)
# quick plot
# adjust margins
par(mar = c(4, 4, 0, 0.2) + 0.2)
# see their distribution
plot(density(inbrs), xlab = 'inbreeding coefficient', main = '')

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
# compute pairwise FST matrix from kinship matrix
pairwise_fst <- pwfst(kinship)
# fancy legend label
leg_title <- expression(paste('Pairwise ', F[ST]))
# NOTE no need for inbr_diag() here!
plot_popkin(
    pairwise_fst,
    labs = subpops,
    labs_even = TRUE,
    labs_line = 1,
    labs_cex = 0.7,
    leg_title = leg_title,
    mar = c(2, 0.2)
)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
indexes_AFR <- subpops == 'AFR'
# AFR subset of the kinship matrix
kinship_AFR <- kinship[indexes_AFR, indexes_AFR]

# kinship matrix plot
plot_popkin(
    inbr_diag( kinship_AFR ),
    mar = 0
)

# estimate FST before rescaling (this value will be wrong, too high!)
fst( kinship_AFR )

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
# rescale kinship_AFR
# since subpops is missing, minimum kinship value is set to zero
# (no averaging between subpopulations)
kinship_AFR <- rescale_popkin( kinship_AFR )

# kinship matrix plot
plot_popkin(
    inbr_diag( kinship_AFR ),
    mar = c(0, 0.3)
)

# FST is now correct, relative to the MRCA of AFR individuals
fst( kinship_AFR )

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'-------------------
# dummy labels to have lines in second panel
subpops_AFR <- subpops[ indexes_AFR ]
plot_popkin(
    # inbr_diag also works on a list of matrices
    inbr_diag( list( kinship, kinship_AFR ) ),
    # title of each panel
    titles = c('All', 'AFR only, rescaled'),
    # pass per-panel labels using a list
    labs = list( subpops, subpops_AFR ),
    # scalar options are shared across panels
    labs_even = TRUE,
    labs_line = 1,
    labs_cex = 0.5,
    # second value is top margin (to make space for titles)
    mar = c(2, 2)
)

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
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

# plotting with different options per level is more complicated...
plot_popkin(
	inbr_diag( kinship ),
	labs = cbind( subpops, blocks ), # ... labs is now a matrix with levels on columns
	labs_even = c(TRUE, FALSE),      # ... even spacing for first level only
	labs_line = c(1, 2),             # ... put second level further out
	labs_cex = c(0.7, 1),            # ... don't shrink second level
	labs_sep = c(FALSE, TRUE),       # ... draw lines inside heatmap for second level only
	ylab_adj = 0.65,                 # push up outer margin ylab "Individuals"
    mar = 3                          # increase margins again
	)

## ---- fig.width = 6, fig.height = 2.8, fig.align = 'center'-------------------
plot_popkin(
    inbr_diag(list(kinship, kinship_AFR)),
	titles = c('All', 'AFR only, rescaled'),
    # list of matrices
	labs = list( cbind(subpops, blocks), subpops_AFR ),
    # non-list: values are reused for both panels
	labs_even = c(TRUE, FALSE),
	labs_line = c(1, 2),
	# make label bigger in second panel (custom per-panel values)
    # list of vectors
	labs_cex = list(c(0.5, 0.7), 1),
	# add lines for first level of second panel (custom per-panel values)
    # list of vectors
	labs_sep = list(c(FALSE, TRUE), TRUE),
    mar = c(3, 2)
	)

## ---- fig.width = 6.5, fig.height = 2.8, fig.align = 'center'-----------------
plot_popkin(
    inbr_diag(list(kinship, kinship_AFR)),
	titles = c('All', 'AFR only, rescaled'),
	labs = list( cbind(subpops, blocks), subpops_AFR ),
	labs_even = c(TRUE, FALSE),
	labs_line = c(1, 2),
	labs_cex = list(c(0.5, 0.7), 1),
	labs_sep = list(c(FALSE, TRUE), TRUE),
    mar = c(3, 2),
    # this option adds a legend per panel
    leg_per_panel = TRUE
	)

## -----------------------------------------------------------------------------
# number of loci (rows)
m_loci <- nrow( X )
# number of subpopulations (columns)
k_subpops <- length( pop_order )
# initialize matrix
P <- matrix( NA, nrow = m_loci, ncol = k_subpops )
# copy names from data
colnames( P ) <- pop_order
rownames( P ) <- rownames( X )
# navigate subpopulations
for ( u in 1 : k_subpops ) {
    # subpopulation label name
    pop <- pop_order[ u ]
    # columns of interest
    js <- subpops == pop
    # now average genotypes into allele frequency estimates and store
    P[ , u ] <- rowMeans( X[ , js ], na.rm = TRUE ) / 2
}

## -----------------------------------------------------------------------------
coancestry <- popkin_af( P )

## ---- fig.width = 4.2, fig.height = 3, fig.align = 'center'-------------------
# coancestry matrix plot
# NOTE: inbr_diag() is not needed for coancestry!
plot_popkin(
    coancestry,
    mar = 3,
    names = TRUE,
    ylab = 'Subpopulations'
)


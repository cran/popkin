# loads Rdata matrices to test
load('Xs.RData')

# construct artificial function that reads m loci at the time from matrix
# input is so tiny only one chunk will be read probably, but let's have a nice general solution just in case...
ml <- 1 # keeps track of last locus not yet read from (will be edited as global)
Xf <- function(m) {
    # set up range of loci to access
    if (ml > nrow(X)) {
        ml <<- 1 # reset after we're done, so next test just works
        return(NULL) # return NULL if we've read everything (tells popkin to stop)
    }
    is <- ml:min( (ml+m-1), nrow(X)) # this gives correct max length of "m"
    # if (max(is) > nrow(X)) { # reduce range if it's out of bounds of matrix
    #     is <- is[is <= nrow(X)] # reduce here
    #     if (length(is) == 0) return(NULL) # if nothing was left, 
    # }
    ml <<- max(is)+1 # update global with next value to read...
    X[is,] # return this subset!
}
n <- ncol(X) # need to pass it separately for function inputs

# only repeat tests where genotypes X is input!

test_that("function returns precomputed values: popkin_A", {
    expect_silent( obj <- popkin_A( Xf, n = n ) )
    expect_equal( obj$A, A )
    expect_equal( obj$M, M )
})

# higher-level tests now!

test_that("function returns precomputed values: popkin", {
    expect_equal(popkin(Xf, n=n), Phi0)
    expect_equal(popkin(Xf, subpops0), Phi0) # NOTE: n==length(subpops0) is inferred
    expect_equal(popkin(Xf, subpops), Phi)   # NOTE: n==length(subpops) is inferred
})

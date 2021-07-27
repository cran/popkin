plot_admix <- function(
                       Q,
                       col = RColorBrewer::brewer.pal( ncol( Q ), "Paired" ),
                       cex_leg = 0.6, # shrink legend labels
                       bty_leg = 'n', # no box for legend
                       # these can be NA to omit
                       xlab = "Individuals",
                       ylab = "Ancestry",
                       ylab_side = 4 # default put on right
                       ) {
    graphics::barplot(
                  t( Q ),
                  col = col,
                  xlab = '',
                  yaxt = 'n', # to be able to put axis on other side
                  border = NA,
                  space = 0,
                  main = '',
                  legend.text = TRUE, # names from Q, though no legend gets added if names are missing
                  args.legend = list(
                      cex = cex_leg,
                      bty = bty_leg
                  )
              )
    # add y-axis
    graphics::axis( ylab_side )
    if ( !is.na( ylab ) )
        graphics::mtext( ylab, side = ylab_side, line = 2 )
    # add x-axis
    if ( !is.na( xlab ) )
        graphics::mtext( xlab, side = 1, line = 2 )
}

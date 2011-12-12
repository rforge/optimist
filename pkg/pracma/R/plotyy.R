plotyy <- function(x, y1, y2, grid = TRUE,
                   type = "l", lwd = 1, lty = 1,
                   col.y1 = "navy", col.y2 = "maroon", ...) {
    stopifnot(is.numeric(x), is.numeric(y1), is.numeric(y2))

    y1range <- pretty(y1)
    y1l <- min(y1range); y1u <- max(y1range)
    y2range <- pretty(y2)
    y2l <- min(y2range); y2u <- max(y2range)

    y2p <- y1l + (y2 - y2l)/(y2u - y2l) * (y1u - y1l)
    y2pretty <- y1l + (y2range - y2l)/(y2u - y2l) * (y1u - y1l)

    opar <- par(mar = c(5.1, 4.1, 4.1, 3.1))
    plot(range(x), range(c(y1, y2p)), type = "n", yaxt = "n",
         xlab = "x", ylab = "y", ...)
    box(col = "grey")
    axis(side = 2, at=NULL, col = col.y1)
    axis(side = 4, at = y2pretty, labels = y2range, col = col.y2)

    if (grid) grid()
    points(x, y1, type = type, col = col.y1, lwd = lwd, lty = lty)
    lines(x, y2p, type = type, col = col.y2, lwd = lwd, lty = lty)
    par(opar)

    invisible()
}

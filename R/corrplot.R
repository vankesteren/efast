#' Create a clean correlation plot
#'
#' This function is a wrapper around corrplot
#'
#' @param mat <matrix> or <list<matrix>> the matrix or list of matrices to be
#'            displayed.
#' @param ... other arguments passed to corrplot::corrplot()
#' @seealso \code{\link[corrplot]{corrplot}}
#'
#' @importFrom corrplot corrplot
#'
#' @export
cplot <- function(mat, ...) {
  if (is.list(mat)) {
    for (i in seq_along(mat)) cplot(mat[[i]], ...)
    return(invisible())
  }
  corrplot(mat, tl.pos = "n", cl.pos = "n", method = "square",
           addgrid.col = FALSE, ...)
}

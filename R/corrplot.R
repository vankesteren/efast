#' Create a clean correlation plot
#'
#' This function is a thin wrapper around corrplot
#'
#' @param mat <matrix> or <list<matrix>> the matrix or list of matrices to be
#'            displayed.
#' @param tl.pos <string> text label position, see corrplot::corrplot()
#' @param cl.pos <string> color label position, see corrplot::corrplot()
#' @param method <string> drawing method, see corrplot::corrplot()
#' @param addgrid.col <string> grid colour see corrplot::corrplot()
#' @param ... other arguments passed to corrplot::corrplot()
#' @seealso \code{\link[corrplot]{corrplot}}
#'
#' @importFrom corrplot corrplot
#'
#' @export
cplot <- function(mat, tl.pos = "n", cl.pos = "n", method = "square",
                  addgrid.col = NA, ...) {
  if (is.list(mat)) {
    for (i in seq_along(mat)) cplot(mat[[i]], ...)
    return(invisible())
  }
  corrplot(mat, tl.pos = tl.pos, cl.pos = cl.pos, method = method,
           addgrid.col = addgrid.col, ...)
}

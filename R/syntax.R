#' model_syntax
#' @rdname model_syntax
#' @keywords internal
efa_model <- function(x, M = 3) {
  # input: dataset or covariance matrix
  # output: efa model with M factors
  nm <- colnames(x)
  mod <- paste(
    "# Exploratory factor model\n# ------------------------\n",
    paste0("efa('efa1')*F", 1:M, collapse = " +\n"),
    "=~",
    paste(nm, collapse = " + ")
  )

  if (is.data.frame(x)) {
    # put relaxed constraints on the residual variances
    mod <- paste0(mod,
                  "\n\n# Uniqueness constraints\n# ----------------------\n")
    for (p in seq_along(nm)) {
      var_p <- stats::var(x[,p])
      nm_p  <- nm[p]
      mod <- paste0(mod, nm_p,
                    " ~~ lower(", -var_p * 0.05, ")*", nm_p,
                    " + upper(", var_p * 1.05, ")*", nm_p,
                    " + v_", p, "*", nm_p, "\n")
    }
  } else {
    # put relaxed constraints on the residual variances
    mod <- paste0(mod,
                  "\n\n# Uniqueness constraints\n# ----------------------\n")
    for (p in seq_along(nm)) {
      var_p <- x[p,p]
      nm_p  <- nm[p]
      mod <- paste0(mod, nm_p,
                    " ~~ lower(", -var_p * 0.05, ")*", nm_p,
                    " + upper(", var_p * 1.05, ")*", nm_p,
                    " + v_", p, "*", nm_p, "\n")
    }
  }
  return(mod)
}

#' str_model
#' @rdname model_syntax
#' @keywords internal
str_model <- function(x, rstruct) {
  # input: dataset or vector of names, pairs of residual structure
  # output: structure model with residual structure
  # first, check if any of rstruct are not in names(dat)
  if (is.data.frame(x) || is.matrix(x))
    datanames <- colnames(x)
  else
    datanames <- x
  rstrnames <- unique(unlist(rstruct))
  # get names in rstruct which are not in the dataset
  errnames <- setdiff(rstrnames, intersect(datanames, rstrnames))
  if (length(errnames) > 0) stop(
    "There are variables in the specified residual structure (rstruct) ",
    "which do not occur in the dataset. \nNames: ",
    paste(errnames, collapse = ", ")
  )

  rescovs <- vapply(
    X         = rstruct,
    FUN       = paste,
    FUN.VALUE = "lhs ~~ rhs",
    collapse  = " ~~ "
  )
  mod <- paste0(
    "# Residual structure\n# ------------------\n",
    paste(rescovs, collapse = "\n")
  )
  mod
}

#' hemi_model
#' @rdname model_syntax
#' @keywords internal
hemi_model  <- function(rois, var_const = FALSE) {
  # input: roi names
  # output: lavaan structure model
  mod <- "# Hemisphere model\n# ----------------\n"
  for (roi in rois) {
    if (var_const) {
      mod <- paste0(mod, "lh_", roi, " ~~ cov_lat*rh_", roi, "\n")
    } else {
      mod <- paste0(mod, "lh_", roi, " ~~ cov_", roi, "*rh_", roi, "\n")
    }
  }
  mod
}

#' lat_index stuff
#' @rdname model_syntax
#' @keywords internal
lat_index <- function(rois) {
  # input: roi names
  # output: lateralisation index per region code
  nroi <- length(rois)

  mod <- "# Lateralization index code\n# -----------------------------------\n"

  for (i in 1:nroi) {
    mod <- paste0(mod, "LI_", rois[i], " := 1 - cov_", rois[i], " / (sqrt(v_", i, ") * sqrt(v_", i + nroi, "))\n")
  }
  mod
}

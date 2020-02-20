#' model_syntax
#' @rdname model_syntax
#' @keywords internal
efa_model <- function(dat, M = 3) {
  # input: dataset
  # output: efa model with M factors
  nm <- colnames(dat)
  mod <- paste(
    "# Exploratory factor model\n# ------------------------\n",
    paste0("efa('efa1')*F", 1:M, collapse = " +\n"),
    "=~",
    paste(nm, collapse = " + ")
  )
  # put relaxed constraints on the residual variances
  mod <- paste0(mod, "\n\n# Uniqueness constraints\n# ----------------------\n")
  for (p in 1:ncol(dat)) {
    var_p <- stats::var(dat[,p])
    nm_p  <- nm[p]
    mod <- paste0(mod, nm_p,
                  " ~~ lower(", -var_p * 0.05, ")*", nm_p,
                  " + upper(", var_p * 1.05, ")*", nm_p,
                  " + v_", p, "*", nm_p, "\n")
  }
  return(mod)
}

#' str_model
#' @rdname model_syntax
#' @keywords internal
str_model <- function(dat, rstruct) {
  # first, check if any of rstruct are not in names(dat)
  datanames <- names(dat)
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
    mod <- paste0(mod, "usum_", rois[i], " := v_", i, " + v_", i + nroi, "\n")
    mod <- paste0(mod, "LI_", rois[i], " := ",
                  "(usum_", rois[i], " - 2*cov_", rois[i], ")",
                  " / usum_", rois[i], "\n")
  }
  mod
}

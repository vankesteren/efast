covmat <- toeplitz(1/1:10)
diag(covmat[6:10, 1:5]) <-
  diag(covmat[1:5, 6:10]) <-
  diag(covmat[6:10, 1:5]) + 0.4
cov2cor(covmat)


diag2 <- function(mat, lower = TRUE) {
  P <- ncol(mat)
  diag(mat[1:(P/2), (P/2 + 1):P])
}

# minres

Asim <- cbind(
  c(rep(0.6, 7), rep(0,3)),
  c(rep(0,5), rep(0.6, 5))
)
AAtsim <- tcrossprod(Asim)
P <- nrow(Asim)
Bsim <- matrix(0, P, P/2)
diag(Bsim) <- diag(Bsim[P:1,(P/2):1]) <- sqrt(0.2)
BBtsim <- tcrossprod(Bsim)
Hsim <- diag(diag(AAtsim + BBtsim))
R <- AAtsim + BBtsim + diag(P) - Hsim

fa <- psych::fa(R, fm = "pa", nfactors = 2)

Imat <- diag(P)
diag(Imat[1:5, 6:10]) <- diag(Imat[6:10, 1:5]) <- 1

# EFAST-PAF, extension of algo 1 from https://sci-hub.tw/10.1007/BF02289468
Ri <- R
P  <- ncol(Ri)
M  <- 2
A <- matrix(0, P, M)
Str <- matrix(0, P, P)
H <- matrix(0, P, P)
for (i in 1:5) {
  # EFA matrix | Uniqueness, Structure
  eig <- eigen(R - diag(P) + H - Str)
  A   <- sqrt(eig$values[1:M]) * eig$vectors[,1:M, drop = FALSE]
  AAt <- tcrossprod(A)

  # Structure matrix | EFA matrix, uniqueness
  Res <- R - diag(P) + H - AAt
  Str <- Res * (Imat - diag(P))
  Str <- Str + diag(colSums(Str))

  # Communality matrix | Structure matrix, EFA matrix
  H <- diag(1 - diag(R - Str - AAt))

  # What to estimate EFA matrix from
  Fit <- crossprod(c(R - AAt - Str - diag(P) + H))
  cat(Fit, "\n")
}

(R - AAt) * Imat

# let's reproduce the fa object
residual <- R * Imat - H

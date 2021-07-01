# @param uv A data.frame object. The id, U, delta, A, and Cat for each 
#   participant with dimension {n x 5}
#
# @param cats An integer vector object. The response categories.
#
# @param txOpts A vector object. The treatment options. 
#
# @param itmax An integer object. The maximum number of iterations for the
#   Newton-Raphson algorithm.
#
# @param tol A numeric object. The value at which the Newton-Raphson is 
#   deemed to have converged.
#
#
# @returns a list object containing
#  \item{alpha} estimated alpha parameters
#  \item{beta} estimated beta parameters
#  \item{score} score matrix
#  \item{khat} Kaplan-Meier estimate
#  \item{respmat} Response indicator matrix
#
#' @include kaplanMeier.R
.ipw <- function(uv, cats, txOpts, itmax, tol) {

  # itmax must be a positive integer
  if (!is.numeric(x = itmax)) stop("itmax must be an integer object", call. = FALSE)
  if (!is.integer(x = itmax)) itmax <- as.integer(x = round(x = itmax, 
                                                            digits = 0L))
  if (itmax <= 0L) stop("inappropriate value of itmax given", call. = FALSE)

  # tol must be a positive numeric; flag if > 0.01 or < 0.000001
  if (!is.numeric(x = tol)) stop("tol must be a numeric object", call. = FALSE)
  if (tol < 0.0) stop("tol must be > 0", call. = FALSE)
  if (tol > 1e-2) warning("tol might be too large")
  if (tol < 1e-6) warning("tol might be too small")

  # number of response categories
  nCat <- length(x = cats)

  # number of alpha parameters
  nAlpha <- nCat - 1L

  # number of participants in the clinical trial
  n <- nrow(x = uv)

  # number of treatment options
  nTx <- length(x = txOpts)

  # number of beta parameters
  nBeta <- nTx - 1L

  # Kaplan-Meier estimates 
  # {n}
  khat <- .kaplanMeier(uv = uv, txOpts = txOpts)

  # initialize gradient matrix
  # {nAlpha + nBeta x nAlpha + nBeta}
  grad <- matrix(data = 0.0, nrow = nAlpha + nBeta, ncol = nAlpha + nBeta)

  # response indicator matrix; note that censored
  # subjects are NA and will be zeroed out by delta; category is an integer so
  # so no need to worry about equality test
  # {n x nAlpha}
  respmat <- outer(X = uv$Cat, Y = cats[-nCat], FUN = "<=")
  respmat[is.na(x = respmat)] <- TRUE

  # starting values for proportional odds model parameters
  # {nAlpha}
  alpha <- .logit(x = colMeans(x = respmat * uv$delta / khat))

  # {nBeta}
  beta <- rep(x = 0.0, times = nBeta)

  # create factor based design matrix to extract appropriate beta parameter
  # {n x nBeta}
  designA <- stats::model.matrix(object = ~as.factor(A), data = uv)
  designA <- unname(obj = designA[,-1L,drop = FALSE])

  # Newton-Raphson algorithm

  cvg <- FALSE

  for (iter in 1L:itmax) {
    
    # form the proportional odds model estimating function (score) and gradient

    # exipt_{ij} = expit(alpha_j + sum_k beta_k I(A_i = a_k)
    # {n x nAlpha}
    expit <- .expit(x = matrix(data = alpha, 
                               nrow = n, 
                               ncol = nAlpha, 
                               byrow = TRUE) + 
                        drop(x = designA %*% beta))

    # E{M_alpha_{j}} = sum_i Delta_i / K(U_i,A_i) (R_j - expit_{ij})
    # {nAlpha}
    M_alpha <- colSums(x = {respmat - expit} * uv$delta / khat)

    # E{M_beta_{k}} = sum_i Delta_i / K(U_i,A_i) I(A_i = a_k) sum_j (R_j - expit_{ij})
    # {nBeta}
    M_beta <- rowSums(x = t(x = designA) %*% {{respmat - expit} * uv$delta / khat})

    # E{M}
    # {nAlpha+nBeta}
    score <- c(M_alpha, M_beta)  

    # Derivatives

    # d/dalpha M_alpha
    # {n x nAlpha}
    pm1d <- uv$delta / khat * expit * {1.0 - expit}

    # E{d/dbeta M_beta}
    # {nBeta}
    tb <- colSums(x = rowSums(x = pm1d) * designA)

    # E{d M / d(alpha, beta)}
    diag(x = grad) <- c(colSums(x = pm1d), tb)

    # E{d/dalpha_j M_beta}
    # E{d/dbeta_k M_alpha}
    # {nBeta x nAlpha}
    tmp <- t(x = designA) %*% pm1d
    for (itx in 1L:nBeta) {
      grad[{1L:nAlpha},nAlpha+itx] <- tmp[itx,]
      grad[nAlpha+itx,{1L:nAlpha}] <- tmp[itx,]
    }

    # parameter update
    incr <- tryCatch(expr = solve(a = grad, b = score),
                     error = function(e) {
                               stop("unable to invert gradient in Newton-Raphson",
                                    e$message, call. = FALSE)
                              })

    alpha <- alpha + incr[1L:nAlpha]
    beta <- beta + incr[nCat:{nAlpha + nBeta}]

    # iteration counter
    iter <- iter + 1L

    if (max(abs(x = score)) < tol) {
      cvg <- TRUE
      break
    }

  }

  if (!cvg) warning("Newton-Raphson did not converge")

  # Return converged parameter values, score, KM estimates, and
  # response indicator matrix
  return( list("alpha" = alpha,
               "beta"  = beta,
               "score" = score,
               "khat" = khat,
               "respmat" = respmat) )
}

.logit <- function(x) { 
  log(x = x / {1.0 - x})
}

.expit <- function(x) {
  exp(x = x) / {1.0 + exp(x = x)}
}


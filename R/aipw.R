# include aumentation terms in estimator
#
# @param uv A data.frame object. The id, U, delta, A, and Cat for each 
#   participant with dimension {n x 5}
#
# @param ti A matrix object or NULL. If a matrix, the bases for the 
#   time-independent component with dimension {n x M+1}. NULL indicates
#   that the time-independent component is to be excluded from the AIPWCC
#   estimator.
#
# @param td A list object or NULL. If a list, each element is a matrix 
#   containing the bases for the time-dependent component with dimension 
#   {n_a x L}. Each element corresponds to a specific tx subgroup. NULL 
#   indicates that the time-dependent component is to be excluded from the
#   AIPWCC estimator.
#
# @param ipwObj The value object returned by .ipw(). A list containing
#   \item{alpha}{The estimated alpha parameters.}
#   \item{beta}{The estimated beta parameters.}
#   \item{score}{The score vector.}
#   \item{khat}{The Kaplan-Meier estimate}
#   \item{respmat}{The response indicator matrix}
#
#
# @param uniqueCensor A list object. Each element is a vector of the unique 
#   censoring times for the specific tx subgroup.
#
# @param txOpts A vector object. The treatment options. The order of the 
#   list elements of td and uniqueCensor correspond to the order of this object.
#
## @returns A list object. The
#   elements of the list correspond to the selected AIPWCC and/or IPWCC
#   estimators. For each estimator, two matrix objects are returned: $logOdds
#   contains the estimated beta parameters, their standard errors estimated
#   using the sandwich estimator, the 95\% confidence intervals,
#   and the p-values for the log odds ratio; $odds
#   contains the estimated odds ratio, their standard errors estimated
#   using the delta method, and the 95\% confidence intervals.
#
#' @importFrom stats qnorm pnorm model.matrix lm predict.lm
#' @include augment.R
.aipw <- function(uv, ti, td, ipwObj, uniqueCensor, txOpts) {

  # number of participants in the clinical trial
  n <- nrow(x = uv)

  # estimated alpha parameters
  # {nAlpha}
  alphaHat <- ipwObj$alpha
  nAlpha <- length(x = alphaHat)

  # estimated beta parameters
  # {nBeta}
  betaHat <- ipwObj$beta
  nBeta <- length(x = betaHat)

  # kaplan-meier estimate
  # {n}
  kHat <- ipwObj$khat

  # response indicator matrix
  # {n x nAlpha}
  R <- ipwObj$respmat

  # create factor based design matrix to extract appropriate beta parameter
  # {n x nBeta}
  designA <- stats::model.matrix(object = ~as.factor(A), data = uv)
  designA <- designA[,-1L,drop = FALSE]

  # expit(alpha_j - sum beta_k I(A = a_k))
  # {n x nAlpha}
  expit <- .expit(x = matrix(data = alphaHat, 
                             nrow = n,  
                             ncol = nAlpha,  
                             byrow = TRUE) + 
                      drop(x = designA %*% betaHat))

  # M(F; alpha, beta)
  # {n x nAlpha+nBeta}
  eM <- cbind(R - expit, designA * rowSums(x = R - expit)) * uv$delta / kHat

  # d M_{alpha} /d alpha)
  # {n x nAlpha}
  tmp <- - uv$delta / kHat * expit * {1.0 - expit}

  # (E{d M_alpha / d alpha})^{-1}
  # {nAlpha x nAlpha}
  B11inv <- diag(x = 1.0/colSums(x = tmp))

  # d M_{alpha} /d beta
  # {nBeta x nAlpha}
  B12T <- t(x = designA) %*% tmp

  # d M_{beta} /d alpha
  # {nAlpha x nBeta}
  B12 <- t(x = B12T)

  # d M_{beta} /d beta
  # {nBeta}
  B22 <- diag(x = rowSums(x = B12T), nrow = nBeta)

  ## beta components of [E{d M/d^T(alpha^T,beta^T)}]^{-1}

  # C A^{-1}
  # {nBeta x nAlpha}
  B12T_B11inv <- B12T %*% B11inv

  # [-C A^{-1}, I]
  # {nBeta x nAlpha+nBeta}
  B12T_B11inv_1 <- cbind(-B12T_B11inv, diag(x = nBeta))

  # m(F; alpha, beta)
  # {n x nBeta}
  score <- eM %*% t(x = B12T_B11inv_1)

  # (D - C A^{-1} B)^{-1}
  # {nBeta x nBeta}
  vee <- -{B22 - B12T_B11inv %*% B12}

  # {nBeta x nBeta}
  veeInv <- tryCatch(expr = solve(a = vee),
                     error = function(e) {
                               stop("unable to invert Fisher Information", 
                                    e$message, 
                                    call. = FALSE)
                             })

  # augmentation components
  #   $adj a matrix of dimension {n x nBeta}
  #   $design a list with matrices of dimension {nSubja x nTI + nTD*nTx}
  aug.results <- .augment(ti = ti,  
                          td = td, 
                          uv = uv,  
                          uniqueCensor = uniqueCensor, 
                          txOpts = txOpts,
                          score = score)

  # augmented score
  # {n x nBeta}
  if (!is.null(x = aug.results$adj)) {
    score <- score + aug.results$adj
  }

  # get the IPW SE {nBeta}
  betaIPW <- betaHat
  se_betaIPW <- veeInv %*% crossprod(x = score) %*% veeInv
  se_betaIPW <- sqrt(x = diag(x = se_betaIPW))

  if (!is.null(x = td) || !is.null(x = ti)) {
    # for each beta, regress the influence function on the design matrix 
    yHat <- matrix(data = 0.0, nrow = n, ncol = ncol(score))

    for (i in 1L:nBeta) {

      fit <- tryCatch(expr = stats::lm(formula = score[,i] ~ -1L + aug.results$design[[ i ]]),
                      error = function(e) {
                                message("error in linear regression step")
                                stop(e$message)
                              })

      yHat[,i] <- predict.lm(object = fit)

    }

    # one-step update and sandwich SE {nBeta}
    betaAIPW <- betaHat - drop(x = colSums(x = yHat) %*% veeInv)
    se_betaAIPW <- veeInv %*% crossprod(x = score - yHat) %*% veeInv
    se_betaAIPW <- sqrt(x = diag(x = se_betaAIPW))
  } else {
    betaAIPW <- NULL
    se_betaAIPW <- NULL
  }

  # prepare log odds ratio results for return

  # beta parameter estimates {nBeta x nEst}
  est <- cbind("IPWCC" = betaIPW, "AIPWCC" = betaAIPW)

  # sandwich estimator of standard errors {nBeta x nEst}
  se <- cbind("IPWCC" = se_betaIPW, "AIPWCC" = se_betaAIPW)

  # 95% confidence intervals {nBeta x nEst}
  lower <- est - stats::qnorm(p = 0.975)*se
  upper <- est + stats::qnorm(p = 0.975)*se

  # p-values {nBeta x nEst}
  pValue <- {1.0 - stats::pnorm(q = abs(x = est / se))}*2.0

  ipwResult <- list("logOdds" = cbind("beta" = est[,"IPWCC"],
                                      "se" = se[,"IPWCC"],
                                      "lower .95" = lower[,"IPWCC"],
                                      "upper .95" = upper[,"IPWCC"],
                                      "p-value" = pValue[,"IPWCC"]))

  rownames(x = ipwResult$logOdds) <- names(x = txOpts)[-1L]

  if (!is.null(x = ti) || !is.null(x = td)) {

    aipwResult <- list("logOdds" = cbind("beta" = est[,"AIPWCC"],
                                         "se" = se[,"AIPWCC"],
                                         "lower .95" = lower[,"AIPWCC"],
                                         "upper .95" = upper[,"AIPWCC"],
                                         "p-value" = pValue[,"AIPWCC"]))
    rownames(x = aipwResult$logOdds) <- names(x = txOpts)[-1L]
  }

  # prepare odds ratio results for return

  # estimated odds ratio {nBeta x nEst}
  est <- exp(x = est)

  # estimated standard error using delta method  {nBeta x nEst}
  dx <- diag(x = est[,"IPWCC"], ncol = nrow(x = est))
  se_IPW <- sqrt(x = diag(x = dx %*% veeInv %*% crossprod(x = score) %*% veeInv %*% dx))

  # 95% confidence interval using confidence interval of parameters
  lower <- exp(x = lower)
  upper <- exp(x = upper)

  ipwResult$odds <- cbind("est" = est[,"IPWCC"],
                          "se" = se_IPW,
                          "lower .95" = lower[,"IPWCC"],
                          "upper .95" = upper[,"IPWCC"])

  rownames(x = ipwResult$odds) <- names(x = txOpts)[-1L]

  result <- list("IPWCC" = ipwResult)


  if (!is.null(x = ti) || !is.null(x = td)) {
    dx <- diag(x = est[,"AIPWCC"], ncol = nrow(x = est))
    se_AIPW <- sqrt(x = diag(x = dx %*% veeInv %*% crossprod(x = score - yHat) %*% veeInv %*% dx))

    aipwResult$odds = cbind("est" = est[,"AIPWCC"],
                            "se" = se_AIPW,
                            "lower .95" = lower[,"AIPWCC"],
                            "upper .95" = upper[,"AIPWCC"])
    rownames(x = aipwResult$odds) <- names(x = txOpts)[-1L]
    result$AIPWCC <- aipwResult
  }

  
  # Return results
  return( result )
}

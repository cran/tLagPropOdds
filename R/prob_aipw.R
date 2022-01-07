# include aumentation terms in estimator for probability of falling into
# a specific category broken down by treatment
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
# @param ipwObj The value object returned by .prob_ipw(). A list containing
#   \item{piHat}{The estimated probabilities.}
#   \item{khat}{The Kaplan-Meier estimate}
#
# @param uniqueCensor A list object. Each element is a vector of the unique 
#   censoring times for the specific tx subgroup.
#
# @param txOpts A vector object. The treatment options. The order of the 
#   list elements of td and uniqueCensor correspond to the order of this object.
#
# @returns A list object. The
#   elements of the list correspond to the selected AIPWCC and/or IPWCC
#   estimators. For each estimator, a list of matrix objects is returned,
#   one for each treatment, 
#   containing the estimated probabilities, their standard errors, and
#   the 95\% confidence intervals.
#
#' @importFrom stats qnorm pnorm model.matrix lm predict.lm
#' @include augment.R
.prob_aipw <- function(uv, ti, td, ipwObj, uniqueCensor, cats, txOpts) {

  # number of participants in the clinical trial
  n <- nrow(x = uv)

  # number of categories
  nCat <- length(x = cats)

  # number of treatments
  nTx <- length(x = txOpts)

  # kaplan-meier estimate
  # {n}
  kHat <- ipwObj$khat

  # create factor based design matrix to extract appropriate beta parameter
  # {n x nTx}
  designA <- stats::model.matrix(object = ~as.factor(A)-1L, data = uv)
  designA <- unname(obj = designA)

  # number of participants in each treatment group
  na <- colSums(x = designA)
  # to avoid dividing by zero, set zeros to 1L
  na[na < 1L] <- 1L

  allResults <- list()
  allResults[[ "IPWCC" ]] <- NULL
  allResults[[ "AIPWCC" ]] <- NULL

  for (iCat in 1L:nCat) {  

    # estimated probabilities for each treatment
    # {nTx}
    piHat <- ipwObj$piHat[,iCat]

    # {n x nTx}
    piMatrix <- matrix(data = piHat, nrow = n, ncol = nTx, byrow = TRUE)

    # M(F; beta)
    # {n x nTx}
    tmp <- {uv$Cat == cats[iCat]}*1.0
    tmp[is.na(x = tmp)] <- 0L
    score <- designA * {tmp - piMatrix} * uv$delta / kHat

    # augmentation components
    #   $adj a matrix of dimension {n x nTx}
    #   $design a list {nTx} with matrices
    aug <- .augment(ti = ti,  
                    td = td, 
                    uv = uv,  
                    uniqueCensor = uniqueCensor, 
                    txOpts = txOpts,
                    score = score)

    score <- score + aug$adj

    # get the IPW SE {nTx}
    piIPW <- piHat

    se_piIPW <- crossprod(x = score)
    se_piIPW <- sqrt(x = diag(x = se_piIPW)) / na

    if (!is.null(x = td) || !is.null(x = ti)) {

      # for each beta, regress the influence function on the design matrix 
      yHat <- matrix(data = 0.0, nrow = n, ncol = nTx)
      for (i in 1L:nTx) {

        fit <- tryCatch(expr = stats::lm(formula = score[,i] ~ -1L + aug$design[[ i ]]),
                        error = function(e) {
                                  message("error in linear regression step")
                                  stop(e$message)
                                })

        yHat[,i] <- predict.lm(object = fit)

      }

      # one-step update and sandwich SE {nBeta}
      piAIPW <- piHat - colSums(x = yHat) / na
      se_piAIPW <- crossprod(x = score - yHat)
      se_piAIPW <- sqrt(x = diag(x = se_piAIPW)) / na
    } else {
      piAIPW <- NULL
      se_piAIPW <- NULL
    }

    # beta parameter estimates {nBeta x nEst}
    est <- cbind("IPWCC" = piIPW, "AIPWCC" = piAIPW)

    # sandwich estimator of standard errors {nBeta x nEst}
    se <- cbind("IPWCC" = se_piIPW, "AIPWCC" = se_piAIPW)

    # 95% confidence intervals {nBeta x nEst}
    lower <- est - stats::qnorm(p = 0.975)*se
    upper <- est + stats::qnorm(p = 0.975)*se

    allResults[[ "IPWCC" ]] <- rbind(allResults[[ "IPWCC" ]],
                                     data.frame("Cat" = cats[iCat], 
                                                "Tx" = txOpts,  
                                                "est" = est[,"IPWCC"],  
                                                "se" = se[,"IPWCC"],  
                                                "lower .95" = lower[,"IPWCC"],  
                                                "upper .95" = upper[,"IPWCC"]))

    rownames(x = allResults$IPWCC) <- NULL

    if (!is.null(x = ti) || !is.null(x = td)) {

      allResults[[ "AIPWCC" ]] <- rbind(allResults[[ "AIPWCC" ]],
                                        data.frame("Cat" = cats[iCat], 
                                                   "Tx" = txOpts,  
                                                   "est" = est[,"AIPWCC"],  
                                                   "se" = se[,"AIPWCC"],  
                                                   "lower .95" = lower[,"AIPWCC"],  
                                                   "upper .95" = upper[,"AIPWCC"]))
      rownames(x = allResults$AIPWCC) <- NULL
    }

  }

  # results are desired grouped by treatment rather than category
  ipw <- allResults$IPWCC
  aipw <- allResults$AIPWCC

  newIPW <- list()
  newAIPW <- list()

  for (i in 1L:nTx) {
    rw <- ipw$Tx == txOpts[i]
    newIPW[[ as.character(x = txOpts[i]) ]] <- ipw[rw,-2L,drop=FALSE]
    rownames(newIPW[[ as.character(x = txOpts[i]) ]]) <- NULL
    if (!is.null(x = aipw)) {
      newAIPW[[ as.character(x = txOpts[i]) ]] <- aipw[rw,-2L,drop=FALSE]
      rownames(newAIPW[[ as.character(x = txOpts[i]) ]]) <- NULL
    }
  }

  allResults$IPWCC <- newIPW
  if (!is.null(x = aipw)) allResults$AIPWCC <- newAIPW



  return( allResults )
}

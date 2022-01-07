# @param uv A data.frame object. The id, U, delta, A, and Cat for each 
#   participant with dimension {n x 5}
#
# @param cats An integer vector object. The response categories.
#
# @param txOpts A vector object. The treatment options. 
#
# @returns a list object containing
#  \item{beta} estimated beta parameters
#  \item{khat} Kaplan-Meier estimate
#
#' @include kaplanMeier.R
.prob_ipw <- function(uv, cats, txOpts) {

  # number of participants in the clinical trial
  n <- nrow(x = uv)

  # number of categories
  nCat <- length(x = cats)

  # number of treatments
  nTx <- length(x = txOpts)

  # Kaplan-Meier estimates 
  # {n}
  khat <- .kaplanMeier(uv = uv, txOpts = txOpts)

  # avoid dividing by zero
  khat[khat < 1e-12] <- 1L

  # response indicator matrix; note that censored
  # subjects are NA and will be zeroed out by delta; category is an integer so
  # no need to worry about equality test
  # {n x nCat}
  respmat <- outer(X = uv$Cat,  Y = cats, FUN = "==")
  respmat[is.na(x = respmat)] <- FALSE

  # create factor based design matrix to extract appropriate probability
  # removal of intecept allows for correct 0/1 for all treatments
  # {n x nTx}
  designA <- stats::model.matrix(object = ~as.factor(A)-1L, data = uv)
  designA <- unname(obj = designA)

  # number of participants in each treatment group
  na <- colSums(x = designA)
  # to avoid dividing by zero, set zeros to 1L
  na[na < 1L] <- 1L

  # {nTx x nCat}
  piHat <- matrix(data = NA, nrow = nTx, ncol = nCat)

  for (iCat in 1L:nCat) {

    # {n x nTx}
    score <- designA * {respmat[,iCat] * uv$delta / khat}

    piHat[,iCat] <- colSums(x = score) / na

  }

  # Return probabilities and KM estimates
  return( list("piHat" = piHat, # {nTx x nCat}
               "khat" = khat) ) # {n}
}

#  Create the martingale contributions for each treatment group. 
#  Calculate the adjustment term.
#
# @param td A list object or NULL. If a list, each element is a matrix 
#   containing the bases for the time-dependent component with dimension 
#   {n_a x L}. Each element corresponds to a specific tx subgroup. NULL 
#   indicates that the time-dependent component is to be excluded from the
#   AIPWCC estimator.
#
# @param ... Ignored. Included to require named inputs
#
# @param uv A data.frame object. The id, U, delta, A, and Cat for each 
#   participant with dimension {n x 5}
#
# @param uniqueCensor A list object. Each element is a vector of the unique 
#   censoring times for the specific tx subgroup.
#
# @param txOpts A vector object. The treatment options. The order of the 
#   list elements of td and uniqueCensor correspond to the order of this object.
#
# @param score A matrix object {n x nBeta}. The score vector for each
#   participant and each treatment option.
# 
#' @import methods
setGeneric(name = ".infl",
           def = function(td, ...) { standardGeneric(".infl") })

# anything not explicitly allowed is forbidden
setMethod(f = ".infl",
          signature = c(td = "ANY"),
          definition = function(td, ...) { stop("not allowed") })

# If only considering the IPWCC, calculate only the adjusment term.
setMethod(f = ".infl",
          signature = c(td = "NULL"),
          definition = function(td, ..., uv, uniqueCensor, txOpts, score) { 

            # number of treatment options
            nTx <- length(x = txOpts)

            # number of beta parameters
            # when being called by catProbs method, this is equivalent
            # to nTx
            nBeta <- ncol(x = score)

            # Initialize returned list for adjustment term
            adj.list <- list()

            for (i in 1L:nTx) {

              # when <5% of participants in tx subset are censored,
              # time-dependent component is not included
              if (any(is.na(x = uniqueCensor[[ i ]]))) {
                adj.list[[ i ]] <- NA
                next
              }

              # identify participants in treatment subgroup
              subja <- uv$A == txOpts[i]

              # number of participants in treatment subgroup
              nSubja <- sum(subja)

              # calculate dMC, Yt, and Ysum using only the participants
              # in the treatment subgroup. Returned object is a list containing
              #   $dMC  with dimension {nSubja x nUniqueCensor}
              #   $Yt   with dimension {nSubja x nUniqueCensor}
              #   $Ysum with dimension {nUniqueCensor}
              dmc <- .dMC(u = uv$U[subja], 
                          delta = uv$delta[subja],
                          uniqueCensor = uniqueCensor[[ i ]],  
                          txOpts = txOpts[i])

              # muHat(m,u,a; alpha, beta)
              # {nSubja x nUniqueCensor}
              adj.list[[ i ]] <- matrix(data = 0.0, nrow = nSubja, ncol = nBeta)

              # when called for catProbs, the only non-zero column
              # is that corresponding to the subset of a_i
              for (j in 1L:nBeta) {
                muKhat <- matrix(data = colSums(x = score[subja,j]*dmc$Yt) / 
                                        dmc$Ysum,
                                 nrow = nSubja, 
                                 ncol = length(x = uniqueCensor[[ i ]]), 
                                 byrow = TRUE)

                adj.list[[ i ]][,j] <- rowSums(x = muKhat*dmc$dMC)
              }

            }

            return( list("infl.list" = NULL, "adj.list" = adj.list) )
          })

# Conditions when time-dependent component is to be included in the AIPWCC
setMethod(f = ".infl",
          signature = c(td = "list"),
          definition = function(td, ..., uv, uniqueCensor, txOpts, score) { 

            # number of treatment options
            nTx <- length(x = txOpts)

            # number of beta parameters
            # when being called by catProbs method, this is equivalent
            # to nTx
            nBeta <- ncol(x = score)

            # Initialize returned lists for adjustment and influence
            adj.list <- list()
            infl.list <- list()

            for (i in 1L:nTx) {

              # when <5% of participants in tx subset are censored,
              # time-dependent component is not included
              if (any(is.na(x = uniqueCensor[[ i ]]))) {
                adj.list[[ i ]] <- NA
                infl.list[[ i ]] <- NA
                next
              }

              # identify participants in treatment subgroup
              subja <- uv$A == txOpts[i]

              # number of participants in treatment subgroup
              nSubja <- sum(subja)

              # calculate dMC, Yt, and Ysum using only the participants
              # in the treatment subgroup. Return object is a list containing
              #   $dMC  with dimension {nSubja x nUniqueCensor}
              #   $Yt   with dimension {nSubja x nUniqueCensor}
              #   $Ysum with dimension {nUniqueCensor}
              dmc <- .dMC(u = uv$U[subja], 
                          delta = uv$delta[subja],
                          uniqueCensor = uniqueCensor[[ i ]],  
                          txOpts = txOpts[i])

              # muHat(m,u,a; alpha, beta)
              # {nSubja x nUniqueCensor}
              adj.list[[ i ]] <- matrix(data = 0.0, nrow = nSubja, ncol = nBeta)

              # when called for catProbs, the only non-zero column
              # is that corresponding to the subset of a_i
              for (j in 1L:nBeta) {
                muKhat <- matrix(data = colSums(x = score[subja,j]*dmc$Yt) / 
                                        dmc$Ysum,
                                 nrow = nSubja, 
                                 ncol = length(x = uniqueCensor[[ i ]]), 
                                 byrow = TRUE)
                adj.list[[ i ]][,j] <- rowSums(x = muKhat*dmc$dMC)
              }

              # number of bases in time-dependent component
              nBasis <- dim(x = td[[ i ]])[3L]
              infl.list[[ i ]] <- matrix(data = 0.0, 
                                         nrow = nSubja, 
                                         ncol = nBasis)

              for (l in 1L:nBasis) {

                # h_{ell}(u,X_i,Lbar_i(u))
                # {nsubja x nUniqueCensor}
                thisL <- td[[ i ]][,,l]

                # muHat(h_{ell},u,a)
                # {nUniqueCensor}
                muHat <- colSums(x = dmc$Yt*thisL)/dmc$Ysum

                # h_{ell}(u,X_i,Lbar_i(u)) - muHat(h_{ell},u,a)
                # {nsubja x nUniqueCensor}
                l.lbar <- sweep(x = thisL, 
                                MARGIN = 2L,  
                                STATS = muHat,  
                                FUN = "-")

                infl.list[[ i ]][,l] <- rowSums(x = l.lbar*dmc$dMC)
              }
            }

            return(list("infl.list" = infl.list, "adj.list" = adj.list) )
          })

# general purpose function to calculate dMC, Yt, and Ysum for a specific
# treatment subgroup. Returns a list containing
#   $dMC  with dimension {nSubja x nUniqueCensor}
#   $Yt   with dimension {nSubja x nUniqueCensor}
#   $Ysum with dimension {nUniqueCensor}
.dMC <- function(u, delta, uniqueCensor, txOpts) {

  # I(U_i = u, Delta_i = 0)
  # {nSubja x nUniqueCensor}
  dNt <- {outer(X = u*{1L-delta}, Y = uniqueCensor+1e-8, FUN = "<") &
          outer(X = u*{1L-delta}, Y = uniqueCensor-1e-8, FUN = ">")} * 1.0

  # sum_{i=1}^{n} I(U_i = u, Delta_i = 0)
  # {nUniqueCensor}
  dNtsum <- colSums(x = dNt)

  # I(U_i >= u)
  # {nSubja x nUniqueCensor}
  Yt <- outer(X = u, Y = uniqueCensor-1e-8, FUN = ">")*1.0

  # sum_{i=1}^{n} I(U_i >= u)
  # {nUniqueCensor}
  Ysum <- colSums(x = Yt) 

  # dMHat_{ci}(u,a)
  # {nSubja x nUniqueCensor}
  dMC <- dNt - sweep(x = Yt, MARGIN = 2L, STATS = dNtsum/Ysum, FUN = "*")

  return( list("dMC" = dMC, "Yt" = Yt, "Ysum" = Ysum) )
}

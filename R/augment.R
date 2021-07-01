# Construct the augmentation terms for each treatment
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
#' @include infl.R
#' @import methods

setGeneric(name = ".augment",
           def = function(ti, td, ...) { standardGeneric(".augment") })

# anything not explicitly allowed is forbidden
setMethod(f = ".augment",
          signature = c(ti = "ANY",
                        td = "ANY"),
          definition = function(ti, td, ...) { stop("not allowed") })

# If only considering the IPWCC, no aumentation term needs to be calculated.
setMethod(f = ".augment",
          signature = c(ti = "NULL",
                        td = "NULL"),
          definition = function(ti, td, ...) { return( NULL ) })

# Conditions when only the time-independent component is to be included in the 
# AIPWCC.
setMethod(f = ".augment",
          signature = c(ti = "matrix",
                        td = "NULL"),
          definition = function(ti, td, ..., uv, uniqueCensor, txOpts, score) { 

            # total number of participants in the study
            n <- nrow(x = uv)

            # number of treatment options
            nTx <- length(x = txOpts)

            # number of beta parameters
            nBeta <- nTx - 1L

            # number of time-independent basis functions
            nTI <- ncol(x = ti)

            # martingale increments for the adjustment and
            # time-dependent augmentation terms
            # returns list containing
            #   infl.list list of matrices {nSubja x L}
            #   adj.list list of vectors {nSubja x nBeta}
            haug <- .infl(uv = uv, 
                          td = NULL,  
                          uniqueCensor = uniqueCensor,  
                          txOpts = txOpts,
                          score = score)

            # design matrix for each treatment subgroup
            des <- list()

            # adjusment term for dependent variable
            adj <- matrix(data = 0.0, nrow = n, ncol = nBeta)

            for (i in 1L:nBeta) {

              des[[ i ]] <- matrix(data = 0.0, nrow = n, ncol = 0L)

              for (j in 2L:nTx) {
                 # design matrices include only the M+1 time-independent covariates
                des[[ i ]] <- cbind(des[[ i ]], -ti * {{uv$A == txOpts[j]} - 
                                               mean(x = {uv$A == txOpts[j]})})
              }


              for (j in 1L:nTx) {
                if (any(is.na(x = haug$adj.list[[ j ]]))) next

                subja <- uv$A == txOpts[j]
                adj[subja,i] <- haug$adj.list[[ j ]][,i]
              }
            }

            return( list("design" = des, "adj" = adj) )

          })

# Conditions when only the time-dependent component is to be included in the 
# AIPWCC.
setMethod(f = ".augment",
          signature = c(ti = "NULL",
                        td = "list"),
          definition = function(ti, td, ..., uv, uniqueCensor, txOpts, score) { 

            # total number of participants in the study
            n <- nrow(x = uv)

            # number of treatment options
            nTx <- length(x = txOpts)

            # number of beta parameters
            nBeta <- nTx - 1L

            # martingale increments for the adjustment and
            # time-dependent augmentation terms
            # returns list containing
            #   infl.list list of matrices {nSubja x L}
            #   adj.list list of vectors {nSubja x nBeta}
            haug <- .infl(td = td,  
                          uv = uv, 
                          uniqueCensor = uniqueCensor,  
                          txOpts = txOpts,
                          score = score)

            # design matrix for each treatment subgroup
            des <- list()

            # adjusment term for dependent variable
            adj <- matrix(data = 0.0, nrow = n, ncol = nBeta)

            # design matrices include only the nTD time-dependent covariates
            # for each treatment option with >5% censoring.

            for (i in 1L:nBeta) {

              des[[ i ]] <- matrix(data = 0.0, nrow = n, ncol = 0L)

              for (j in 1L:length(x = txOpts)) {

                if (any(is.na(x = haug$infl.list[[ j ]]))) next

                # number of time-dependent basis functions
                nTD <- ncol(x = haug$infl.list[[ j ]])

                tmp <- matrix(data = 0.0, nrow = n, ncol = nTD)

                subja <- uv$A == txOpts[j]

                tmp[subja,] <- haug$infl.list[[ j ]]
               
                des[[ i ]] <- cbind(des[[ i ]], tmp)

                adj[subja,i] <- haug$adj.list[[ j ]][,i]
              }
            }

            return( list("design" = des, "adj" = adj) )

          })

# Conditions for full AIPWCC
setMethod(f = ".augment",
          signature = c(ti = "matrix",
                        td = "list"),
          definition = function(ti, td, ..., uv, uniqueCensor, txOpts, score) { 

            # total number of participants in the study
            n <- nrow(x = uv)

            # number of treatment options
            nTx <- length(x = txOpts)

            # number of beta parameters
            nBeta <- nTx - 1L

            # number of time-independent basis functions
            nTI <- ncol(x = ti)

            # martingale increments for the adjustment and
            # time-dependent augmentation terms
            # returns list containing
            #   infl.list list of matrices {nSubja x L}
            #   adj.list list of vectors {nSubja x nBeta}
            haug <- .infl(td = td,  
                          uv = uv, 
                          uniqueCensor = uniqueCensor,  
                          txOpts = txOpts,
                          score = score)

            # design matrix for each treatment subgroup
            des <- list()

            # adjusment term for dependent variable
            adj <- matrix(data = 0.0, nrow = n, ncol = nBeta)

            # design matrices include the nTI time-indenpent covariates and
            # nTD time-dependent covariates for each treatment option

            for (i in 1L:nBeta) {

              des[[ i ]] <- matrix(data = 0.0, nrow = n, ncol = 0L)

              for (j in 2L:nTx) {
                 # design matrices include only the M+1 time-independent covariates
                des[[ i ]] <- cbind(des[[ i ]], -ti * {{uv$A == txOpts[j]} - 
                                               mean(x = {uv$A == txOpts[j]})})
              }

              for (j in 1L:nTx) {

                if (any(is.na(x = haug$infl.list[[ j ]]))) next

                # number of time-dependent basis functions
                nTD <- ncol(x = haug$infl.list[[ j ]])

                tmp <- matrix(data = 0.0, nrow = n, ncol = nTD)

                subja <- uv$A == txOpts[j]

                tmp[subja,] <- haug$infl.list[[ j ]]

                des[[ i ]] <- cbind(des[[ i ]], tmp)

                adj[subja,i] <- haug$adj.list[[ j ]][,i]
              }
            }

            return( list("design" = des, "adj" = adj) )

          })

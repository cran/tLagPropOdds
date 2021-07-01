# Kaplan-Meier estimates of the censoring distribution evaluated at the U
# of each subject
#
# @param uv A data.frame object. The id, U, delta, A, and Cat for each 
#   participant with dimension {n x 5}
#
# @param txOpts A vector object. The treatment options. 
#
# @returns A vector object {n}. 
#
.kaplanMeier <- function(uv, txOpts) {

  # total number of participants in the study
  n <- nrow(x = uv)

  # number of treatment options
  nTx <- length(x = txOpts)

  # initialize returned vector
  Khat <- numeric(length = n)
    
  # loop through each treatment
  for (i in 1L:nTx) {

    # identify participants in treatment subgroup
    subja <- uv$A == txOpts[i]

    # censoring distribution evaluated at the U
    ss <- summary(object = survival::survfit(formula = Surv(U, {1L-delta}) ~ 1, 
                                             data = uv[subja,]), 
                  times = uv$U[subja])

    for (j in 1L:length(x = ss$time)) {
      # match evaluation times to participants
      use <- {uv$U > {ss$time[j] - 1e-8}} & {uv$U < {ss$time[j] + 1e-8}}
      Khat[use & subja] <- ss$surv[j]
    }

  }

  return( Khat )
}   

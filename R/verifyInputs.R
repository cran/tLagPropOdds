# @param data A data.frame object. A data.frame containing all observed data.
#    At a minimum, this data.frame must contain columns with headers 
#    "id", "U", "delta", "Cat" and "A". If the time-independent component of
#    the estimator is to be included, data.frame must also contain the 
#    bases of f(X). If the time-dependent component is included, data.frame
#    must also contain the bases of h(X,L) as well as the time intervals with
#    column headers {"tstart", "tstop"} or {"start","stop"}.
#
# @param ti A character or integer vector or NULL. The columns of data to be
#   included in the time-independent component of the estimator, 
#   f_m(X) m = 0, ..., M. If NULL, time-independent component is excluded 
#   from the AIPWCC estimator. If td = NULL, the IPW estimate is
#   calculated. If td != NULL, the IPW and AIPW2 estimates are calculated.  
#
# @param td A character or integer vector or NULL. The columns of data to be
#   included in the time-dependent component of the estimator, 
#   h_l(X,Lbar), l = 1, ..., L. If NULL, the time-dependent component is 
#   excluded from the AIPWCC estimator. If ti = NULL, the IPW estimate is
#   calculated. If ti != NULL, the IPW and AIPW1 estimates are calculated. 
#
# @returns a list containing
#   \describe{
#      \item{uv}{A data.frame with columns "id", "U", "delta", "Cat", and 
#                "A". {n x 5}}
#      \item{cats}{A vector object containing the possible outcome categories.}
#      \item{txOpts}{A vector of the treatment options}
#      \item{uniqueCensor}{A list with an element for each treatment option.
#                Each element contains a vector of the unique censoring times 
#                for the treatment subgroup.}
#      \item{ti}{A data.frame containing the M+1 bases of f(X). Note that
#                there is not an "id" column.}
#      \item{td}{A list with an element for each treatment option. Each
#                element contains the data.frame of bases for the time-dependent 
#                component rebinned according to the unique u values. Note that
#                there is not an "id" column.}
#   }
#
#' @importFrom dplyr distinct %>%
#' @importFrom survival coxph survfit Surv
.verifyInputs <- function(df, ti, td) {

  message("\ninputs processed as follows:")

  # returned object
  inputs <- list()

  ## general data

  # df must be a complete data.frame
  if (!is.data.frame(x = df)) df <- as.data.frame(x = df)

  # columns that make up the "uv" internal data.frame
  required <- c("id", "U", "delta", "Cat", "A")

  # data.frame must contain columns with headers
  # id, U, delta, Cat, A
  if (!all(required %in% colnames(x = df))) {
    stop("data.frame must contain column headers ", 
         paste(required, collapse=", "),
         call. = FALSE)
  }

  # create "uv" data.frame and remove duplicate rows
  df_uv <- df[,required]
  df_uv <- df_uv %>% dplyr::distinct()

  # note that we remove the requirement for Cat to be complete to
  # allow for censored individuals to be set as NA
  if (any(is.na(x = df[,c("id", "U", "delta", "A")]))) {
    stop("data must be complete", call. = FALSE)
  }

  # ensure that delta is integer or can be converted to integer without
  # loss of information
  if (is.factor(x = df_uv$delta)) {
    df_uv$delta <- .convertToInt(x = df_uv$delta)
  } else if (!is.integer(x = df_uv$delta)) {
    if (is.numeric(x = df_uv$delta)) {
      tmp <- as.integer(x = round(x = df_uv$delta, digits = 0L))
      if (!isTRUE(x = all.equal(target = tmp, current = df_uv$delta))) {
        stop("delta must be integer or factor", call. = FALSE)
      }
      df_uv$delta <- tmp
    } else {
      stop("delta must be integer or factor", call. = FALSE)
    }
  }

  # ensure that delta is binary 0/1
  if (!all(df_uv$delta %in% c(0L,1L))) {
    stop("delta must be integer 0/1", call. = FALSE)
  }

  # ensure that Cat is integer or can be converted to integer without
  # loss of information

  if (is.factor(x = df_uv$Cat)) {

    # this converts factors to integers. Censored cases are set as NA
    # all others are integer valued 1L:nCat
    df_uv$Cat <- .convertToInt(x = df_uv$Cat, delta = df_uv$delta)

  } else if (!is.integer(x = df_uv$Cat)) {

    # if numeric, convert to integer and ensure no loss of information

    if (is.numeric(x = df_uv$Cat)) {

      tmp <- as.integer(x = round(x = df_uv$Cat, digits = 0L))

      if (!isTRUE(x = all.equal(target = tmp, current = df_uv$Cat))) {
        stop("Cat must be integer or factor", call. = FALSE)
      }

      df_uv$Cat <- tmp

    } else {

      stop("Cat must be integer or factor", call. = FALSE)

    }
  }

  # generate message indicating number of categories
  # note that the NA of the censored cases is not included in this count

  cats <- sort(x = unique(x = df_uv$Cat[df_uv$delta == 1L]))
  ncats <- length(x = cats)

  if (ncats <= 0L) {

    stop("no outcome categories could be identified", call. = FALSE)

  } else if (ncats == 1L) {

    stop("only 1 outcome category could be identified", call. = FALSE)

  }

  message("  identified ", ncats, " outcome categories")

  inputs$cats <- cats

  # ensure that A is integer or can be converted to integer without
  # loss of information

  if (is.factor(x = df_uv$A)) {

    # store original treatment names for use in final printing
    origTxNames <- levels(x = df_uv$A)

    # this converts factors to integers. All cases are integer valued 1L:nTx
    df_uv$A <- .convertToInt(x = df_uv$A, 
                             delta = rep(x = 1.0, times = nrow(x = df_uv)))

  } else if (!is.integer(x = df_uv$A)) {

    # if numeric, convert to integer and ensure no loss of information

    if (is.numeric(x = df_uv$A)) {

      tmp <- as.integer(x = round(x = df_uv$A, digits = 0L))

      if (!isTRUE(x = all.equal(target = tmp, current = df_uv$A))) {
        stop("treatment (A) must be integer or factor", call. = FALSE)
      }

      df_uv$A <- tmp

      origTxNames <- sort(unique(x = df_uv$A))

    } else {

      stop("treatment (A) must be integer or factor", call. = FALSE)

    }
  } else {
    origTxNames <- sort(unique(x = df_uv$A))
  }    

  # return general data as 'uv'
  inputs$uv <- df_uv

  # treatment options

  txOpts <- sort(x = unique(x = df_uv$A))

  nTxOpts <- length(x = txOpts)

  names(x = txOpts) <- origTxNames

  message("  identified ", nTxOpts, " treatment options")

  # return treatment options as txOpts
  inputs$txOpts <- txOpts

  # unique censoring times for each treatment subgroup
  uniqueCensor <- list()

  isCovariateAdjusted <- TRUE

  for (i in 1L:nTxOpts) {

    # identify participants in treatment subgroup
    subja <- df_uv$A == txOpts[i]

    # report total number of participants in treatment
    nSubja <- sum(subja)
    message("    ", nSubja, " participants received tx ", txOpts[i])

    # fit a NULL model to obtain unique censoring times
    fit <- survival::coxph(Surv(U, {1L-delta}) ~ 1 , data = df_uv[subja,])

    ss <- survival::survfit(formula = fit)

    censoringTimes <- ss$time[ss$n.event != 0L]
    nct <- length(x = censoringTimes)

    if (nct == 0L) {

      # if no censoring in subgroup set uniqueCensor to NA.
      # if time-dependent basis provided, notify that it will not be
      # included for this subgroup

      if (!is.null(x = td)) {
        message("      no participants that received tx A = ", names(txOpts)[i], 
                " are censored;\n",
                "      **** time-dependent component for tx subset removed ****")
      }

      uniqueCensor[[ i ]] <- NA

    } else if (nct < as.integer(x = round(x = nSubja*0.05, digits = 0L))) {

      # if <5% censoring in subgroup set uniqueCensor to NA.
      # if time-dependent basis provided, notify that it will not be
      # included for this subgroup

      if (!is.null(x = td)) {

        message("      fewer than 5% of participants that received tx A = ", 
                names(txOpts)[i], 
                " are censored;\n",
                "      **** time-dependent component for tx subset removed ****")

      }

      uniqueCensor[[ i ]] <- NA

    } else {

      uniqueCensor[[ i ]] <- ss$time[ss$n.event != 0L]

      # flag that this is not just a covariate adjustment
      isCovariateAdjusted <- FALSE

    }
  }

  # if no censoring or <5% censoring, turn off time-dependent component
  if (isCovariateAdjusted) td <- NULL

  # generate message indicating what estimator is being calculated
  if (isCovariateAdjusted) {

    if (is.null(x = ti)) {

      message("  IPWCC estimator")

      type <- "IPWCC"

    } else {

      message("  IPWCC and partial, time-independent AIPWCC estimator")

      type <- c("IPWCC", "partial, time-independent AIPWCC")

    }

  } else if (is.null(x = ti) && is.null(x = td)) {

    message("  IPWCC estimator")

    type <- "IPWCC"

  } else if (is.null(x = td)) {

    message("  IPWCC and partial, time-independent AIPWCC estimators")

    type <- c("IPWCC", "partial, time-independent AIPWCC")

  } else if (is.null(x = ti)) {

    message("  IPWCC and partial, time-dependent AIPWCC estimators")

    type <- c("IPWCC", "partial, time-dependent AIPWCC")

  } else {

    message("  IPWCC and full AIPWCC estimators")

    type <- c("IPWCC", "full AIPWCC")

  }

  # return type for printing procedure
  inputs$type <- type

  # return list of unique censoring times for each treatment subgroup
  inputs$uniqueCensor <- uniqueCensor

  # time-independent component

  if (!is.null(x = ti)) {

    # if base is a numeric/integer object, pull appropriate column headers
    if (!is.character(x = ti)) {
      if (is.numeric(x = ti)) {

        # ensure that all elements of ti are > 0 and < # of columns
        if (max(ti) > ncol(x = df) || any(ti <= 0L)) {
          stop("inappropriate column index provided for input ti", 
               call. = FALSE)
        }

        ti <- colnames(x = df)[ti]

      } else {

        stop("if provided, ti must be a character or integer vector", 
             call. = FALSE)

      }
    }

    # ensure that none of the required headers are provided in ti
    if (any(required %in% ti)) {

      ti <- ti[!(ti %in% required)]

      if (length(x = ti) == 0L) {
        stop("inappropriate column index provided for input ti", 
             call. = FALSE)
      }

    }

    # ensure that the provided column headers are found in the data
    if (!all(ti %in% colnames(x = df))) {

      stop("ti basis function(s) ", 
           paste(ti[!(ti %in% colnames(x = df))], collapse = ", "),
           " not found in provided data.frame", call. = FALSE)

    }

    # extract time-independent bases and remove duplicate rows and id column
    df_ti <- df[,c("id",ti)]
    df_ti <- df_ti %>% dplyr::distinct()
    df_ti <- df_ti[, -1L, drop = FALSE]

    # test if f_0 was included in the model, if not include
    tstIntercept <- apply(X = df_ti, MARGIN = 2L, 
                          FUN = function(x){ 
                                  tst <- all(as.integer(x = round(x = x, digits = 0L)) == 1L)
                                  return( tst )
                                })

    if (!any(tstIntercept)) {

      message("  Intercept (f_0) included as a time-independent basis function")

      df_ti <- cbind("(Intercept)" = 1.0, df_ti)

    }

    if (nrow(x = df_ti) != nrow(x = df_uv)) {

      stop("number of rows for time-independent bases does not agree ",
           "with number of participants", call. = FALSE)

    }

    message("  ", ncol(x = df_ti), " bases in time-independent component")

    df_ti <- data.matrix(frame = df_ti)

  } else {

    df_ti <- NULL

  }

  inputs$ti <- df_ti

  #  time dependent bases

  if (!is.null(x = td)) {

    # if td is a numeric/integer object, pull appropriate column headers
    if (!is.character(x = td)) {
      if (is.numeric(x = td)) {

        # ensure that all elements of td are > 0 and < # of columns
        if ({max(td) > ncol(x = df)} || any(td <= 0L)) {
          stop("inappropriate column index provided for input td", 
               call. = FALSE)
        }

        td <- colnames(x = df)[td]

      } else {
        stop("if provided, td must be a character or integer vector", 
             call. = FALSE)
      }
    }

    # include 'stop' and 'start' in required column headers
    requiredTD <- c(required, "stop", "start")

    # internally data.frame must contain columns with headers
    # id, U, delta, Cat, A, start, stop (or tstart/tstop)
    if (!all(requiredTD %in% colnames(x = df))) {

      # allow for the original column headers generated by survival's tmerge()
      requiredTD <- c(required, "tstop", "tstart")

      if (!all(requiredTD %in% colnames(x = df))) {

        stop("if td != NULL, data must contain columns with headers ",
             "{'start', 'stop'} or {'tstart', 'tstop'}", call. = FALSE)

      } else {

        # rename columns to stop/start for internal use
        rename <- which(x = colnames(x = df) == "tstart")
        colnames(x = df)[rename] <- "start"

        rename <- which(x = colnames(x = df) == "tstop")
        colnames(x = df)[rename] <- "stop"

      }
        
    }

    # ensure that none of the required headers are provided in td
    if (any(requiredTD %in% td)) {

      td <- td[!(td %in% requiredTD)]

      if (length(x = td) == 0L) {
        stop("inappropriate column index provided for input td", 
             call. = FALSE)
      }

    }

    # ensure that the provided column headers are found in the data
    if (!all(td %in% colnames(x = df))) {

      stop("td basis function(s) ", 
           paste(td[!(td %in% colnames(x = df))], collapse = ", "),
           " not found in provided data.frame", call. = FALSE)

    }

    # extract time dependent bases
    df_td <- df[,c("id", "start", "stop", td)]

    # number of bases functions for time dependent component
    nBases <- length(x = td)
    message("  ", nBases, " bases in time-dependent component")

    # lists will contain unique censoring times and bases for each tx subgroup
    tdBases <- list()

    for (j in 1L:nTxOpts) {

      if (any(is.na(x = uniqueCensor[[ j ]]))) next

      # identify participants in treatment subgroup
      subja <- df_uv$A == txOpts[j]

      # number of participants in treatment subgroup
      nSubja <- sum(subja)

      # number of uniqueCensoring times in treatment subgroup
      nUniqueCensor <- length(x = uniqueCensor[[ j ]])

      # matrix to hold the bases of the time dependent bases
      # {nSubja x nUniqueCensor x nBases}

      tdOnU <- array(data = 0.0, dim = c(nSubja, nUniqueCensor, nBases))

      icnt <- 1L
      for (i in 1L:nrow(x = df_uv)) {

        if (!subja[i]) next

        # identify the rows of time dependent data for participant i
        idCol <- which(x = df_td$id == df_uv$id[i])
        nTimeIntervals <- length(x = idCol)

        # determine which value to use for each value of u

        cnt <- integer(length = nUniqueCensor) + nTimeIntervals + 1L

        for (k in 1L:nTimeIntervals) {
          cnt <- cnt - {uniqueCensor[[ j ]] < {df_td$stop[idCol[k]] + 1e-8}}
        }

        # nTimeIntervals+1 x nBases
        vals <- data.matrix(frame = rbind(df_td[idCol,], 
                                          df_td[idCol[nTimeIntervals],]))

        # nUniqueCensor x nBases
        tmp <- unname(obj = vals[cnt,][,td,drop=FALSE])

        tdOnU[icnt,,] <- tmp
        icnt <- icnt + 1L

      }

      dimnames(x = tdOnU) <- list(df_uv$id[subja], NULL, td)
      tdBases[[ j ]] <- tdOnU
    }
  } else {
    tdBases <- NULL
  }

  # return time-dependent dataset as td
  inputs$td <- tdBases

  return( inputs )

}

# general purpose function to convert factors to integers
# x is the factor variable
# delta is a 0/1 vector indicating which elements of x should be included
# in determining the level values
.convertToInt <- function(x, delta) {

  # participants that should be used in determining the levels
  use <- delta != 0L

  # default returned vector to NA
  tt <- rep(x = NA, times = length(x = x))

  # identify levels of factor
  levs <- levels(x = x)
  if (length(x = levs) == 0L) {
    stop("no outcome categories could be identified", call. = FALSE)
  }

  # initialize first level to be 1
  nL <- 1L
  for (i in 1L:length(x = levs)) {

    # identify those that have this level and should be included
    subUse <- {x == levs[i]} & use

    # if no one is in the subgroup cycle to next level without incrementing 
    # level count.
    if (sum(subUse) == 0L) {
      message("    no participants with outcome ", levs[i],
              "; category removed from analysis")
      next
    }

    # set this subgroup to current level count
    tt[subUse] <- nL

    # increment level count for next subgroup
    nL <- nL + 1L
  }

  return( tt )
}

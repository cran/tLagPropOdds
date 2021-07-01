#' Print Analysis Results
#'
#' Prints the key results.
#'
#' @param x A tLagObj object. The value returned by tLagPropOdds().
#'
#' @param ... Ignored. 
#'
#' @name print
#' @examples
#'
#' data(tLagData)
#'
#' # full AIPWCC estimator
#' res <- tLagPropOdds(data = tLagData, 
#'                     ti = "x",  
#'                     td = c("hospStatus", "daysOut"))
#'
#' print(x = res)
#' @method print tLagObj
#' @export
#' @importFrom R.utils printf

print.tLagObj <- function(x, ...) {



  headerLog = "    %4s  %6.4f  %6.4f  %10.4f  %10.4f  %7.4f\n"
  headerLogTitle = "    %4s  %6s  %6s  %10s  %10s  %7s\n"

  headerLog1 = "    %4s  %6.4f  %6.4f  %10.4f  %10.4f\n"
  headerLogTitle1 = "    %4s  %6s  %6s  %10s  %10s\n"

  for (j in 1L:length(x = x)) {

    cat("\n\n", attr(x = x, which = "type")[j], " estimator\n\n", sep="")

    cat("\n  log odds ratio\n\n")
    R.utils::printf(headerLogTitle, "Tx", "beta", "se", "lower .95", "upper .95", "p-value")
    for (i in 1L:nrow(x = x[[ j ]]$logOdds)) {
      R.utils::printf(headerLog, rownames(x[[ j ]]$logOdds)[i], 
                      x[[ j ]]$logOdds[i,1], 
                      x[[ j ]]$logOdds[i,2], 
                      x[[ j ]]$logOdds[i,3], 
                      x[[ j ]]$logOdds[i,4], 
                      x[[ j ]]$logOdds[i,5])
    }

    cat("\n  odds ratio\n\n")
    R.utils::printf(headerLogTitle1, "Tx", "est", "se", "lower .95", "upper .95")
    for (i in 1L:nrow(x = x[[ j ]]$odds)) {
      R.utils::printf(headerLog1, rownames(x[[ j ]]$odds)[i], 
                      x[[ j ]]$odds[i,1], 
                      x[[ j ]]$odds[i,2], 
                      x[[ j ]]$odds[i,3], 
                      x[[ j ]]$odds[i,4])
    }

  }

}

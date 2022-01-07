#' Print Analysis Results
#'
#' Prints the key results.
#'
#' @param x A catProbs object. The value returned by catProbs().
#'
#' @param ... Ignored. 
#'
#' @name print
#' @examples
#'
#' data(tLagData)
#'
#' # full AIPWCC estimator
#' res <- catProbs(data = tLagData, 
#'                 ti = "x",  
#'                 td = c("hospStatus", "daysOut"))
#'
#' print(x = res)
#' @method print catProbsObj
#' @export
#' @importFrom R.utils printf

print.catProbsObj <- function(x, ...) {

  headerLog = "    %4s  %6.4f  %6.4f  %10.4f  %10.4f\n"
  headerLogTitle = "    %4s  %6s  %6s  %10s  %10s\n"

  for (j in 1L:length(x = x)) {

    cat("\n\n", attr(x = x, which = "type")[j], " estimator\n\n", sep="")

    for (k in 1L:length(x = x[[ j ]])) {
      cat("\n  P(Cat = c | A = ", names(x = x[[ j ]])[k], ")\n")
      R.utils::printf(headerLogTitle, "c", "piHat", "se", "lower .95", "upper .95")
      for (i in 1L:nrow(x = x[[ j ]][[ k ]])) {
        R.utils::printf(headerLog,  
                        x[[ j ]][[ k ]][i,1L], 
                        x[[ j ]][[ k ]][i,2L], 
                        x[[ j ]][[ k ]][i,3L], 
                        x[[ j ]][[ k ]][i,4L], 
                        x[[ j ]][[ k ]][i,5L])
      }
    }

  }

  return( NULL )

}

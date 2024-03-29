#' Toy Dataset For Illustration
#'
#' These data are provided for the purposes of illustrating the use of
#' the software. Though the data were generated under a scenario similar
#' to a real-world COVID-19 therapeutics clinical trial, they should not be  
#' interpreted as representing true clinical trial data.
#' 
#' @usage data(tLagData)
#'
#' @format tLagData is a time-dependent data.frame containing the following 
#'   information for 602 participants ascertained at day 90 of a fictitious
#'   randomized clinical trial.
#'   \describe{
#'     \item{id:}{A unique participant identifier.}
#'     \item{A:}{The treatment received, where A=\{0,1\}.}
#'     \item{Cat:}{The ordered outcome category. There are 6 categories 
#'                 ascertained at day 90.
#'                 \describe{
#'                   \item{1:}{ at home and off oxygen, number of days >= 77;} 
#'                   \item{2:}{ at home and off oxygen, number of days 49-76;}
#'                   \item{3:}{ at home and off oxygen, number of days 1-48;}
#'                   \item{4:}{ not hospitalized and either at home on oxygen or not home;}
#'                   \item{5:}{ hospitalized for medical care or in hospice care; and}
#'                   \item{6:}{ dead.}
#'                 }
#'                 If participant is censored, Cat = NA.}
#'     \item{U:}{The time at which the outcome category was determined or
#'               the censoring time. For Cat = 1-5, U is the interim analysis
#'               time (90 days). For Cat = 6, U is the time of death.
#'               For Cat = NA, U is the censoring time.}
#'     \item{delta:}{The event indicator (1 if U is the time at which the
#'                   outcome category was determined;
#'                   0 if censored).}
#'     \item{x:}{A continuous baseline covariate.}
#'     \item{start:}{The lower bound of the time interval to which the
#'                   given covariate values pertain.}
#'     \item{stop:}{The upper bound of the time interval to which the
#'                   given covariate values pertain.}
#'     \item{hospStatus:}{A time-dependent indicator of hospital status, where 1
#'                        indicates that the participant was not in the hospital
#'                        during interval (start, stop]; 0 otherwise.}
#'     \item{daysOut:}{The expected number of continuous days out of hospital
#'                     at the time of the interim analysis (90 days).}
#'  }
#'
#' @name tLagData
#' @keywords datasets
NULL

#' @useDynLib AWKMT2cpp
#' @importFrom Rcpp sourceCpp
#' @import survival
NULL


#' Adaptively Weighted Kaplan-Meier Tests
#'
#' Performs the two-sample tests based on adaptively weighted differences between two Kaplan-Meier curves proposed by Uno, Tian, Claggett and Wei (2015).
#'
#'
#' @param indata          A data matrix (data frame). The 1st column is time-to-event variable, the 2nd column is event indicator (1=event, 0=censor), and the 3rd column is the treatment indicator (1=treat, 0=control). No missing values are allowed in this data matrix.
#' @param tau             A numeric value to specify the time interval of interest. The end of study time will be a general choice.
#' @param c_first         A first number in range to specify the search area of "c" for the versatile tests by Uno et al. (2015). Default is \code{0}.
#' @param c_last          A last number in range to specify the search area of "c" for the versatile tests by Uno et al. (2015). Default is \code{4}.
#' @param c_by            A number to specify the search area of "c" for the versatile tests by Uno et al. (2015). Default is \code{0.1}.
#' @param method          A name of the resampling method. It supports \code{"permutation"} (default) and \code{"perturbation"}.
#' @param nmethod         A number of iterations for the resampling. Recommended to specify at least \code{10000} (default) or larger.
#' @param seed            An integer value, used for the random number generation in the resampling procedures. Default is \code{1}.
#' @param test            Specify \code{"1_side"} for the one-sided test where the alternative hypothesis is that treatment group is superior to control group with respect to survival.
#'                        Specify \code{"2_side"} for the two-sided test where the alternative hypothesis is that treatment group is not equal to control group with respect to survival.
#'                        Default is \code{"1_side"}.
#' @return A list with components:
#' @return \item{test_method}{One-sided or two-sided test}
#' @return \item{resampling_method}{The resampling method.}
#' @return \item{crude_pvalue_V1}{Crude p-value of the test based on v1 in Uno et al. (2015).}
#' @return \item{crude_pvalue_V2}{Crude p-value of the test based on v2 in Uno et al. (2015).}
#' @return \item{bona_fide_pvalue_V1}{Bona-fide p-value of the test based on v1 in Uno et al. (2015).}
#' @return \item{bona_fide_pvalue_V2}{Bona-fide p-value of the test based on v2 in Uno et al. (2015).}
#'
#' @references
#'  Uno H, Tian L, Claggett B, Wei LJ. A versatile test for equality of two survival functions based on weighted differences of Kaplan-Meier curves.
#'  Statistics in Medicine 2015, 34, 3680-3695.
#' @export
AWKMT2cpp <- function(indata, tau, c_first=0, c_last=4, c_by=0.1, method="permutation", nmethod=10000, seed=1, test="1_side"){
  rank2 <- function(x){
    rank(x, ties.method="min")
  }

  if(! method %in% c("permutation", "perturbation")){
    stop("Only permutation and perturbation method are supported")
  }

  if(! test %in% c("1_side", "2_side")){
    stop("Only 1_side and 2_side test are supported")
  }

  #-- Get tau1 --
  tau1 = 0

  #-- Get tau2 --
  tau2 = tau

  #-- Get unique_time, n_times --
  t_idx   = unique(sort(c(indata$time, 0, tau1, tau2, tau)))
  n_times = length(t_idx)

  #-- Get crange --
  crange   = seq(c_first, c_last, by=c_by)

  x0 <- subset(indata, arm == 0)$time
  x1 <- subset(indata, arm == 1)$time

  status0 <- subset(indata, arm == 0)$status
  status1 <- subset(indata, arm == 1)$status

  if(! is.null(seed) ) set.seed(seed)

  res <- AWKMT2_eigen(t_idx = as.numeric(t_idx),
                      status0 = as.numeric(status0),
                      status1 = as.numeric(status1),
                      x0 = as.numeric(x0),
                      x1 = as.numeric(x1),
                      tau1 = as.numeric(tau1), tau2 = as.numeric(tau2),
                      crange = as.numeric(crange), test = test, type = method, nmethod = nmethod)

  #-- Get crude p-value --
  cc_T1 = rbind(res$V1, res$V1_per)
  cc_T2 = rbind(res$V2, res$V2_per)

  #-- T1 --
  dd_T1    = apply(-cc_T1, 2, rank2)-1
  d2_T1    = dd_T1[1, ]
  crude_T1 = min(d2_T1/nmethod)

  #-- T2 --
  dd_T2    = apply(-cc_T2, 2, rank2)-1
  d2_T2    = dd_T2[1, ]
  crude_T2 = min(d2_T2/nmethod)

  #-- Get the null distribution of pb to choose the threshold value for claiming a statistical significance based on pb --
  #-- T1 --
  ee_T1        = apply(- res$V1_per, 2, rank2)-1
  ours_bona_T1 = apply(ee_T1, 1, min)/(nmethod-1)

  aa_crude_pvalue_T1  = rep(crude_T1, nmethod)
  ours_bona_pvalue_T1 = mean(ours_bona_T1 < aa_crude_pvalue_T1)

  #-- T2 --
  ee_T2 = apply(-res$V2_per, 2, rank2)-1
  ours_bona_T2 = apply(ee_T2, 1, min)/(nmethod-1)

  aa_crude_pvalue_T2  = rep(crude_T2, nmethod)
  ours_bona_pvalue_T2 = mean(ours_bona_T2 < aa_crude_pvalue_T2)

  #-- Output --
  output = list()
  output$test_method                = test
  output$resampling_method          = method
  output$crude_pvalue_V1               = crude_T1
  output$crude_pvalue_V2               = crude_T2
  output$bona_fide_pvalue_V1           = ours_bona_pvalue_T1
  output$bona_fide_pvalue_V2           = ours_bona_pvalue_T2
  output
}

#' Conducts a DiProPerm test
#'
#' @name DiProPerm
#'
#' @description This package conducts a Direction-Projection-Permutation (DiProPerm) test.
#' For more details see Wei et al, 2016 \href{https://www.tandfonline.com/doi/full/10.1080/10618600.2015.1027773}{here}.
#' DiProPerm assesses whether a binary linear classifier can detect a difference between two high-dimensional distributions.
#'
#' @param X An \code{nxp} data matrix.
#' @param y A vector of \code{n} binary class labels -1 and 1.
#' @param B The number of permutations for the DiProPerm test. The default is 1000.
#' @param classifier A string designating the binary linear classifier. classifier="dwd", distance weighted discrimination, is the default. classifier="dwd" implements a generalized DWD model from the \code{\link[DWDLargeR]{genDWD}} function in the \code{DWDLargeR} package.
#' The penalty parameter, \code{C}, in the \code{genDWD} function is calculated using the \code{\link[DWDLargeR]{penaltyParameter}} function in \code{DWDLargeR}. More details on the algorithm used to calculate the DWD solution can be found \href{https://www.tandfonline.com/doi/full/10.1080/10618600.2017.1366915}{here}.
#' One other option for the binary classifier is the "md", mean difference direction.
#' @param univ.stat A string indicating the univariate statistic used for the projection step. univ.stat="md", mean difference, is the default.
#' @param balance A logical indicator for whether a balanced permutation design should be implemented.  The default is TRUE.
#' @param alpha An integer indicating the level of significance. The default is 0.05.
#' @param cores An integer indicating the number of cores to be used for parallel processing. The default is 2. Note, parallel processing is only available on MacOS and Ubuntu operating systems at this time. Windows users will default to using 1 core.
#'
#' @return A list containing:
#' \item{\code{X}}{The observed \code{nxp} data matrix.}
#' \item{\code{y}}{The observed vector of \code{n} binary class labels -1 and 1.}
#' \item{\code{obs_teststat}}{The observed univariate test statistic.}
#' \item{\code{xw}}{Projection scores used to compute the specified univariate statistic.}
#' \item{\code{w}}{The loadings of the binary classification.}
#' \item{\code{Z}}{The \code{Z} score of the observed test statistic.}
#' \item{\code{cutoff_value}}{The cutoff value to achieve an alpha level of significance.}
#' \item{\code{pvalue}}{The pvalue from the permutation test.}
#' \item{\code{perm_dist}}{A list containing the permuted projection scores and permuted class labels for each permutation.}
#' \item{\code{perm_stats}}{A \code{B} dimensional vector of univariate test statistics.}
#'
#' @export
#'
#' @import DWDLargeR
#' @importFrom usethis use_package
#'
#' @author Andrew G. Allmon, J.S. Marron, Michael G. Hudgens
#' @examples
#' \donttest{data(mushrooms)
#' X <- Matrix::t(mushrooms$X)
#' y <- mushrooms$y
#' dpp <- DiProPerm(X=X,y=y,B=10)}
#' \dontshow{data(mushrooms)
#' X.temp <- SparseM::as.matrix(mushrooms$X)[,1:50]
#' X <- Matrix::t(X.temp)
#' y <- mushrooms$y[1:50]
#' dpp <- DiProPerm(X=X,y=y,B=100,classifier="md")}
#'
#' @references
#' {Lam, X. Y., Marron, J. S., Sun, D., & Toh, K.-C. (2018). Fast Algorithms for
#'  Large-Scale Generalized Distance Weighted Discrimination. Journal of
#'  Computational and Graphical Statistics, 27(2), 368–379.
#'  \url{https://doi.org/10.1080/10618600.2017.1366915}}
#'
#' {Wei, S., Lee, C., Wichers, L., & Marron, J. S. (2016). Direction-Projection-Permutation for High-Dimensional Hypothesis Tests. Journal of Computational and Graphical Statistics, 25(2), 549–569. \url{https://doi.org/10.1080/10618600.2015.1027773}}
DiProPerm <- function(X,y,B=1000,classifier="dwd",univ.stat="md",balance=TRUE,alpha=0.05,cores=2) {

  X.t <- Matrix::t(X)

  if(is.na(match("matrix.csr",class(X)))){X.t <- SparseM::as.matrix.csr(X.t)}

########### Step 1: Calculate hyperplane using classifier direction #####################
  ######## DWD direction ##############
  if (classifier=="dwd") {

    C = DWDLargeR::penaltyParameter(X.t,y,expon=1,rmzeroFea = 0)
    # solve the generalized DWD model
    result = DWDLargeR::genDWD(X.t,y,C=C,expon=1,rmzeroFea = 0) ## Iain uses C=0.1 in his example
    w.obs <- result$w / norm_vec(result$w)
  }

  ######### Mean Difference Direction ########
  if (classifier=="md") {
    X.temp <- SparseM::as.matrix(X)
    w.md <- apply(X.temp[y==-1,],2,mean)-apply(X.temp[y==1,],2,mean)
    w.obs <- w.md / norm_vec(w.md)
  }

########### Step 2: Caluculate two sample univariate meand difference statistic #####################
  xw.obs <- X %*% w.obs
  if (univ.stat=="md") {
    obs.teststat <- abs(mean(xw.obs[y==1])-mean(xw.obs[y==-1])) ## Very close to obs.test.md which is in DPP python package
  }

########### Step 3: Do a permutation test #####################

  ## This set is the best and fastest: takes about 30 seconds when C is from the penalty param function from Marron et al. ##
  if (classifier=="dwd") {
    perm_list <- RepParallel(B,n.cores=cores,dwd_scores(X.t,n=nrow(X),balance))
    permdist_rsamp <- unlist(lapply(perm_list,sumfun_diffmean))
  }

  if (classifier=="md") {
    perm_list <- RepParallel(B,n.cores=cores,md_scores(X.temp,n=nrow(X),balance))
    permdist_rsamp <- unlist(lapply(perm_list,sumfun_diffmean))
  }

  pval <- mean(permdist_rsamp>=obs.teststat)

  pct <- stats::quantile(permdist_rsamp, probs = (1-alpha))

  Z <- (obs.teststat-mean(permdist_rsamp))/stats::sd(permdist_rsamp)

  return(list(X=X,y=y,obs_teststat=obs.teststat,xw=xw.obs,w=w.obs,Z=Z,cutoff_value=pct[[1]],pvalue=pval,perm_dist=perm_list,perm_stats=permdist_rsamp))
}

#source("R/dpp_functions.R")

## Data must be in pxn format for Marron's code ##
#data(mushrooms)
#X <- t(mushrooms$X)
#y <- mushrooms$y
#set.seed(21083)
#perm <- DiProPerm(X=X,y=y,B=10)



#unlist(lapply(perm,sumfun_diffmean))
#X.t <- t(X)
##X.t <- as(Matrix(X.t,sparse=TRUE),"matrix.csr")

#C = penaltyParameter(X.t,y,expon=1)
# solve the generalized DWD model
#result = genDWD(X.t,y,C=C,expon=1) ## Iain uses C=0.1 in his example
#w.obs <- result$w / norm_vec(result$w)
#xw.obs <- X %*% w.obs
#obs.teststat <- abs(mean(xw.obs[y==1])-mean(xw.obs[y==-1])) ## Very close to obs.test.md which is in DPP python package

#dwd_scores()
#### Possible Code to use exact distribution ####
#### For Exact Distribution #####

#ind_comb <- combinations(nrow(X), sum(y==1))
#nrow(ind_comb) ## count combinations

#head(ind_comb) ## look at the first few
#permdist_diffmean <- apply(ind_comb,
#                           MARGIN=1,function(x) sumfun_diffmean(simfun(x)))
#2*mean(permdist_diffmean>=sumfun_diffmean(ants))
#################################################


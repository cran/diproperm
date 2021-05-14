#' Conducts a DiProPerm test
#'
#' @name DiProPerm
#'
#' @description This package conducts a Direction-Projection-Permutation (DiProPerm) test.
#' DiProPerm is a two-sample hypothesis test for comparing two high-dimensional
#' distributions. The DiProPerm test is exact, i.e., the type I error is guaranteed
#' to be controlled at the nominal level for any sample size. For more details see Wei et al. (2016).
#'
#' @param X An \code{nxp} data matrix.
#' @param y A vector of \code{n} binary class labels -1 and 1.
#' @param B The number of permutations for the DiProPerm test. The default is 1000.
#' @param classifier A string designating the binary linear classifier. classifier="dwd", distance weighted discrimination (DWD), is the default. classifier="dwd" implements a generalized DWD model from the \code{\link[DWDLargeR]{genDWD}} function in the \code{DWDLargeR} package.
#' The penalty parameter, \code{C}, in the \code{genDWD} function is calculated using the \code{\link[DWDLargeR]{penaltyParameter}} function in \code{DWDLargeR}. The \code{genDWD} and \code{penaltyParameter} functions have several arguments which are set to recommended default values. More details on the algorithm used to calculate the DWD solution can be found in Lam et al. (2018).
#' Other options for the binary classifier include the "md", mean difference direction, and "svm", support vector machine. The "svm" option uses the default implementation from \code{\link[e1071]{svm}}.
#' @param univ.stat A string indicating the univariate statistic used for the projection step. univ.stat="md", the mean difference, is the default.
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
#' \dontshow{data(mushrooms)
#' X <- Matrix::t(mushrooms$X)
#' X <- X[1:50,]
#' y <- mushrooms$y[1:50]
#' dpp <- DiProPerm(X=X,y=y,B=100,classifier="svm")}
#'
#' @references
#' {Lam, X. Y., Marron, J. S., Sun, D., & Toh, K.-C. (2018). Fast Algorithms for
#'  Large-Scale Generalized Distance Weighted Discrimination. Journal of
#'  Computational and Graphical Statistics, 27(2), 368–379.
#'  \doi{10.1080/10618600.2017.1366915}}
#'
#' {Wei, S., Lee, C., Wichers, L., & Marron, J. S. (2016). Direction-Projection-Permutation for High-Dimensional Hypothesis Tests. Journal of Computational and Graphical Statistics, 25(2), 549–569. \doi{10.1080/10618600.2015.1027773}}
DiProPerm <- function(X,y,B=1000,classifier="dwd",univ.stat="md",balance=TRUE,alpha=0.05,cores=2) {

  #Xbefore <<- FALSE
  #ybefore <<- FALSE

  #if(eval(is.element('X',ls()),envir = .GlobalEnv)){Xbefore <<- TRUE}
  #if(eval(is.element('y',ls()),envir = .GlobalEnv)){ybefore <<- TRUE}

  X <<- X
  y <<- y

  ## Perform error checks
  if(B<1000) {message("Setting the number of permutations, B, below 1000 is not recommended.")}
  if(dim(table(y)) != 2){stop("Error: y vector must be 2 dimensions with '-1' and '1' for class labels.")}
  if(sum(is.element(c(1,-1),y)) != 2){stop("y vector must have class labels of '-1' and '1'.")}
  if(is.na(match("matrix",class(X))) & is.na(match("matrix.csr",class(X))) ) {stop("Error: Please make sure your input data, X, is a data matrix.")}

  X.t <- Matrix::t(X)

  if(is.na(match("matrix.csr",class(X)))){X.t <- SparseM::as.matrix.csr(X.t)}

########### Step 1: Calculate hyperplane using classifier direction #####################
  ######## DWD direction ##############
  if (classifier=="dwd") {

    C = quiet(DWDLargeR::penaltyParameter(X.t,y,expon=1,rmzeroFea = 0))
    # solve the generalized DWD model
    result.out = utils::capture.output(DWDLargeR::genDWD(X.t,y,C=C,expon=1,rmzeroFea = 0))
    result = quiet(DWDLargeR::genDWD(X.t,y,C=C,expon=1,rmzeroFea = 0)) ## Iain uses C=0.1 in his example
    cat(result.out[1],result.out[3:8],sep = "\n") ## Exclude the primfeas, dualfeas, and relative gap information from output ##
    w.obs <- result$w / norm_vec(result$w)
  }

  ######### Mean Difference Direction ########
  if (classifier=="svm") {
    ## y must be in factor format for svm to do classification
    y.temp <- as.factor(y)
    result <- e1071::svm(X,y.temp, kernel = "linear")
    w.svm <- Matrix::as.matrix(drop(t(result$coefs)%*%X[result$index,]))
    w.obs <- w.svm[1,] / norm_vec(w.svm[1,])
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

  if (classifier=="svm") {
    perm_list <- RepParallel(B,n.cores=cores,svm_scores(X,n=nrow(X),balance))
    permdist_rsamp <- unlist(lapply(perm_list,sumfun_diffmean))
  }

  if (classifier=="md") {
    perm_list <- RepParallel(B,n.cores=cores,md_scores(X.temp,n=nrow(X),balance))
    permdist_rsamp <- unlist(lapply(perm_list,sumfun_diffmean))
  }

  pval <- paste0(round(mean(permdist_rsamp>=obs.teststat),3))

  if(pval=="0") {pval="<0.001"}

  pct <- stats::quantile(permdist_rsamp, probs = (1-alpha))

  Z <- (obs.teststat-mean(permdist_rsamp))/stats::sd(permdist_rsamp)

  odat <- X
  ydat <- y

  #if(Xbefore==FALSE){rm(X,envir = .GlobalEnv)}
  #if(ybefore==FALSE){rm(y,envir = .GlobalEnv)}

  return(list(X=odat,y=ydat,obs_teststat=obs.teststat,xw=xw.obs,w=w.obs,Z=Z,cutoff_value=pct[[1]],pvalue=pval,perm_dist=perm_list,perm_stats=permdist_rsamp))
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


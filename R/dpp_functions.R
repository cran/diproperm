#########################################
### Author: Drew Allmon
## Purpose: Functions for creating DiProPerm package
#########################################
#devtools::install_github("elbamos/clusteringdatasets")
#devtools::install_github("grayclhn/dbframe-R-library")
usethis::use_package("Matrix")
usethis::use_package("SparseM")
usethis::use_package("DWDLargeR")
usethis::use_package("e1071")
usethis::use_package("parallel")
#usethis::use_package("devtools")
#usethis::use_package("pracma")
usethis::use_package("ggplot2")
usethis::use_package("lemon")
usethis::use_package("gridExtra")
usethis::use_package("dplyr")
#usethis::use_package("methods")
usethis::use_package("sampling")

## Parallel version of replicate taken from dbframe on github ##
RepParallel <- function(B, n.cores=2 , expr, simplify = "array",...) {

  if(Sys.info()["sysname"] == 'Windows') {
      n.cores=1
      message("Note: Parallel processing is not available for Windows machines at this time. Thus, only 1 core will be used.")
    }

  answer <- parallel::mclapply(integer(B), eval.parent(substitute(function(...) expr)),mc.cores = n.cores,...)
  if (!identical(simplify, FALSE) && length(answer))
    return(simplify2array(answer, higher = (simplify == "array")))
  else return(answer)
}

## Calculates the norm of a vector as done in python DPP package
norm_vec <- function(a) sqrt(sum(a^2))

## Suppresses output from DiProPerm iterations ##
quiet <- function(b) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(b))
}

## Conduct a balanced permutation ##
dwd_rsamp <- function(balance,n) {

  if(balance==TRUE) {
    ## Initalize Vecotr to hold the class labels
    simy <- vector(length = n)

    ## Reassign strata values to match what can be read into balancedstratification function
    strata <- y
    strata[strata==-1] <- 2

    n1 <- table(strata)[[1]]  ## Number in class 1
    n2 <- table(strata)[[2]]  ## Number in class -1

    # matrix of balancing variables
    index=cbind(seq(1:length(strata)))

    # Vector of inclusion probabilities.
    pik=rep(n1/n,times=n)

    ## Permform balanced permutation
    s=sampling::balancedstratification(index,strata,pik,comment = FALSE)

    # the sample is
    c1 <- (1:length(pik))[s==1]
    c2 <- (1:length(pik))[s==0]

    ## Reassign previous class labels to work with DWD
    simy[c1] <- 1
    simy[c2] <- -1

    return(simy)
  }

  if(balance==FALSE) {

    ## Performs random permutation
    return(sample(y,prob = rep(0.5,n)))
  }
  #perm <- randperm(y)
  #return(perm)
}

## Calculates the DWD scores ##
dwd_scores <- function(X.t,n,balance) {
  set.seed(NULL)
  perm_y <- dwd_rsamp(balance,n)

  ## Calculate the penalty parameter to be used for DiProPerm ##
  C = quiet(DWDLargeR::penaltyParameter(X.t,perm_y,expon=1,rmzeroFea = 0))

  # solve the generalized DWD model
  result = quiet(DWDLargeR::genDWD(X.t,perm_y,C=C,expon=1,rmzeroFea = 0)) ## Iain uses C=0.1 in his example

  ## Calculate Permuted Scores ##
  w <- result$w / norm_vec(result$w) ## Loadings of Separating Hyperplane
  xw <- X %*% w  ## Projected scores onto hyperplane

  return(list(data.frame(xw,perm_y)))
}

## Calculates the SVM scores ##
svm_scores <- function(Xtemp,n,balance) {
  set.seed(NULL)
  perm_y <- dwd_rsamp(balance,n)
  perm_y_temp <- as.factor(perm_y)

  # solve the SVM model
  result = e1071::svm(Xtemp,perm_y_temp, kernel="linear")

  w.svm <- Matrix::as.matrix(drop(t(result$coefs)%*%Xtemp[result$index,]))

  ## Calculate Permuted Scores ##
  w <- w.svm[1,] / norm_vec(w.svm[1,]) ## Loadings of Separating Hyperplane
  xw <- X %*% w  ## Projected scores onto hyperplane

  return(list(data.frame(xw,perm_y)))
}

## Calculates the mean difference direction
md_scores <- function(X.temp,n,balance) {
  set.seed(NULL)
  perm_y <- dwd_rsamp(balance,n)

  w.md <- abs(apply(X.temp[perm_y==-1,],2,mean)-apply(X.temp[perm_y==1,],2,mean))
  w <- w.md / norm_vec(w.md)
  xw <- X.temp %*% w

  return(list(data.frame(xw,perm_y)))
}

## Calculates the mean difference univariate statistic ##
sumfun_diffmean <- function(dat) {
  with(dat,abs(mean(xw[perm_y==1])-mean(xw[perm_y==-1])))
  #with(dat,abs(mean(xw[y==1])-mean(xw[y==-1])))
}

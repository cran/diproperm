#' Plots diagnostics from DiProPerm test
#'
#' @name plotdpp
#'
#' @description This function plots the diagnostics of a DiProPerm test including the projection scores for the observed data, projection scores for the permutations with the smallest and largest univariate statistic values, and permutation distribution for \code{B} univariate statistics.
#'
#' @param dpp A DiProPerm object.
#' @param plots A string designating the desired plots to be displayed:
#' "obs" displays the projection scores for the observed data,
#' "min" displays the projection scores for the permutation with the smallest univariate statistic value,
#' "max" displays the projection scores for the permutation with the largest univariate statistic value,
#' "permdist" displays the permutation distribution for \code{B} univariate statistics, and
#' "all" displays all 4 diagnostic plots in one plot.
#' Additionally, one can specify "perm1" to display the projection scores for the first permutation and "perm2" to display the projection scores for the second permutation.
#' @param w An integer indicating the width of the jitter. The default is 0.001.
#' @param h An integer indicating the height of the jitter.  The default is 0.001.
#'
#' @return A ggplot
#'
#' @export
#'
#' @author Andrew G. Allmon, J.S. Marron, Michael G. Hudgens
#' @examples
#' \donttest{data(mushrooms)
#' X <- Matrix::t(mushrooms$X)
#' y <- mushrooms$y
#' dpp <- DiProPerm(X=X,y=y,B=10)
#' plotdpp(dpp)}
#' \dontshow{data(mushrooms)
#' X.temp <- SparseM::as.matrix(mushrooms$X)[,1:50]
#' X <- Matrix::t(X.temp)
#' y <- mushrooms$y[1:50]
#' dpp <- DiProPerm(X=X,y=y,B=100,classifier="md")
#' plotdpp(dpp)}
#' \dontshow{data(mushrooms)
#' X <- Matrix::t(mushrooms$X)
#' X <- X[1:50,]
#' y <- mushrooms$y[1:50]
#' dpp <- DiProPerm(X=X,y=y,B=100,classifier="svm")
#' plotdpp(dpp)}

plotdpp <- function(dpp,plots="all",w=0.001,h=0.001) {

  pdata <- as.data.frame(dpp$perm_stats)
  names(pdata) <- c("perm_stats")

  xwdata <- as.data.frame(dpp$xw[,1])
  names(xwdata) <- c("xw")

  permdata1 <- as.data.frame(dpp$perm_dist[[1]])
  permdata2 <- as.data.frame(dpp$perm_dist[[2]])

  minindex <- match(min(dpp$perm_stats),dpp$perm_stats)
  maxindex <- match(max(dpp$perm_stats),dpp$perm_stats)

  permmin <- as.data.frame(dpp$perm_dist[[minindex]])
  permmax <- as.data.frame(dpp$perm_dist[[maxindex]])

  names(permdata1) <- c("xw","permy")
  names(permdata2) <- c("xw","permy")

  names(permmin) <- c("xw","permy")
  names(permmax) <- c("xw","permy")

  ## Determines the height and width of the jitter plot ##
  jitter <- ggplot2::position_jitter(width = w, height = h)

  ## Permutation Distribution Plot ##
  p.diag1 <- ggplot2::ggplot(pdata,ggplot2::aes(x=perm_stats)) +
    ggplot2::geom_density(alpha=0.2,fill="yellow") + ggplot2::geom_vline(xintercept=dpp$obs_teststat,colour="brown",linetype="dashed") +
    #ggplot2::geom_point(position = jitter, ggplot2::aes(y=0.05,color="darkgreen")) +
    ggplot2::scale_color_manual(values = c("darkgreen"),guide=F) +
    ggplot2::labs(y="Density",x="Permutation Statistics") +
    ggplot2::theme_minimal()

  t1 <- ggplot2::ggplot_build(p.diag1)
  minx <- min(t1$data[[1]]["x"])
  maxx <- max(dpp$obs_teststat,dpp$perm_stats)
  miny <- min(t1$data[[1]]["y"])
  maxy <- max(t1$data[[1]]["y"])

  p.diagt1 <- p.diag1 + ggplot2::geom_point(position = jitter, ggplot2::aes(y=maxy*0.25,color="darkgreen")) +
              ggplot2::geom_label(
              label=paste("Z:", round(dpp$Z,3), ", Obs. Test Stat:", round(dpp$obs_teststat,3), ", Cutoff Value:", round(dpp$cutoff_value,3),", Pvalue:",dpp$pvalue),
              x=maxx*0.95,
              y=maxy*0.95,
              vjust=1,
              hjust=1,
              label.padding = ggplot2::unit(0.2, "lines"), # Rectangle size around label
              label.size = 0.5,
              size=3,
              color = "black",
              fill="#69b3a2")

  p.diagt11 <- p.diag1 + ggplot2::geom_point(position = jitter, ggplot2::aes(y=maxy*0.25,color="darkgreen")) +
              ggplot2::geom_label(
                label=paste("Z:", round(dpp$Z,3),"\nObs. Test Stat:", round(dpp$obs_teststat,3), "\nCutoff Value:", round(dpp$cutoff_value,3),"\n","Pvalue:",dpp$pvalue,""),
                x=maxx*0.95,
                y=maxy*0.95,
                vjust=1,
                hjust=1,
                label.padding = ggplot2::unit(0.55, "lines"), # Rectangle size around label
                label.size = 1,
                size=3.25,
                color = "black",
                fill="#69b3a2")

  ## Permutation 1 Plot ##
  p.diag2 <- ggplot2::ggplot(permdata1,ggplot2::aes(x=xw,fill=as.character(permy))) +
    ggplot2::geom_density(alpha=0.2) +
    #ggplot2::geom_point(position = jitter, ggplot2::aes(y=0.05,shape=as.character(permy),color=as.character(permy))) +
    ggplot2::scale_fill_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values=c("red", "blue")) +
    ggplot2::scale_shape_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c(1,4)) +
    ggplot2::scale_color_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c("darkred","darkblue")) +
    ggplot2::scale_x_continuous(limits = c(min(permdata1$xw)-3,max(permdata1$xw)+3)) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse=TRUE),shape=ggplot2::guide_legend(reverse = TRUE),color=ggplot2::guide_legend(reverse=TRUE)) +
    ggplot2::labs(y="Density",x="Projection Score (Perm=1)") +
    ggplot2::theme_minimal()

  t2 <- ggplot2::ggplot_build(p.diag2)
  minx <- min(t2$data[[1]]["x"])
  maxx <- max(t2$data[[1]]["x"])
  miny <- min(t2$data[[1]]["y"])
  maxy <- max(t2$data[[1]]["y"])

  p.diagt2 <- p.diag2 + ggplot2::geom_point(position = jitter, ggplot2::aes(y=maxy*0.25,shape=as.character(permy),color=as.character(permy)))

  ## Permutation 2 Plot ##
  p.diag3 <- ggplot2::ggplot(permdata2,ggplot2::aes(x=xw,fill=as.character(permy))) +
    ggplot2::geom_density(alpha=0.2) +
    #ggplot2::geom_point(position = jitter, ggplot2::aes(y=0.05,shape=as.character(permy),color=as.character(permy))) +
    ggplot2::scale_fill_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values=c("red", "blue")) +
    ggplot2::scale_shape_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c(1,4)) +
    ggplot2::scale_color_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c("darkred","darkblue")) +
    ggplot2::scale_x_continuous(limits = c(min(permdata2$xw)-3,max(permdata2$xw)+3)) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse=TRUE),shape=ggplot2::guide_legend(reverse = TRUE),color=ggplot2::guide_legend(reverse=TRUE)) +
    ggplot2::labs(y="Density",x="Projection Score (Perm=2)") +
    ggplot2::theme_minimal()

  t3 <- ggplot2::ggplot_build(p.diag3)
  minx <- min(t3$data[[1]]["x"])
  maxx <- max(t3$data[[1]]["x"])
  miny <- min(t3$data[[1]]["y"])
  maxy <- max(t3$data[[1]]["y"])

  p.diagt3 <- p.diag3 + ggplot2::geom_point(position = jitter, ggplot2::aes(y=maxy*0.25,shape=as.character(permy),color=as.character(permy)))


  ## Observed Data Plot ##
  p.diag4 <- ggplot2::ggplot(xwdata,ggplot2::aes(x=xw,fill=as.character(y))) +
    ggplot2::geom_density(alpha=0.2) +
    #ggplot2::geom_point(position = jitter, ggplot2::aes(y=0.05,shape=as.character(y),color=as.character(y))) +
    ggplot2::scale_fill_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values=c("red", "blue")) +
    ggplot2::scale_shape_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c(1,4)) +
    ggplot2::scale_color_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c("darkred","darkblue")) +
    ggplot2::scale_x_continuous(limits = c(min(xwdata$xw)-3,max(xwdata$xw)+3)) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse=TRUE),shape=ggplot2::guide_legend(reverse = TRUE),color=ggplot2::guide_legend(reverse=TRUE)) +
    ggplot2::labs(y="Density",x="Projection Score (Observed)") +
    ggplot2::theme_minimal()

  t4 <- ggplot2::ggplot_build(p.diag4)
  minx <- min(t4$data[[1]]["x"])
  maxx <- max(t4$data[[1]]["x"])
  miny <- min(t4$data[[1]]["y"])
  maxy <- max(t4$data[[1]]["y"])

  p.diagt4 <- p.diag4 + ggplot2::geom_point(position = jitter, ggplot2::aes(y=maxy*0.25,shape=as.character(y),color=as.character(y)))

  ## Minimum Permutation Plot ##
  p.diag5 <- ggplot2::ggplot(permmin,ggplot2::aes(x=xw,fill=as.character(permy))) +
    ggplot2::geom_density(alpha=0.2) +
    #ggplot2::geom_point(position = jitter, ggplot2::aes(y=0.05,shape=as.character(permy),color=as.character(permy))) +
    ggplot2::scale_fill_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values=c("red", "blue")) +
    ggplot2::scale_shape_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c(1,4)) +
    ggplot2::scale_color_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c("darkred","darkblue")) +
    ggplot2::scale_x_continuous(limits = c(min(permmin$xw)-3,max(permmin$xw)+3)) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse=TRUE),shape=ggplot2::guide_legend(reverse = TRUE),color=ggplot2::guide_legend(reverse=TRUE)) +
    ggplot2::labs(y="Density",x="Projection Score (Perm=Min)") +
    ggplot2::theme_minimal()

  t5 <- ggplot2::ggplot_build(p.diag5)
  minx <- min(t5$data[[1]]["x"])
  maxx <- max(t5$data[[1]]["x"])
  miny <- min(t5$data[[1]]["y"])
  maxy <- max(t5$data[[1]]["y"])

  p.diagt5 <- p.diag5 + ggplot2::geom_point(position = jitter, ggplot2::aes(y=maxy*0.25,shape=as.character(permy),color=as.character(permy)))


  ## Maximum Permutation Plot ##
  p.diag6 <- ggplot2::ggplot(permmax,ggplot2::aes(x=xw,fill=as.character(permy))) +
    ggplot2::geom_density(alpha=0.2) +
    #ggplot2::geom_point(position = jitter, ggplot2::aes(y=0.05,shape=as.character(permy),color=as.character(permy))) +
    ggplot2::scale_fill_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values=c("red", "blue")) +
    ggplot2::scale_shape_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c(1,4)) +
    ggplot2::scale_color_manual(name="Class Labels",limits=c("1","-1"),labels=c("1","-1"),values = c("darkred","darkblue")) +
    ggplot2::scale_x_continuous(limits = c(min(permmax$xw)-3,max(permmax$xw)+3)) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse=TRUE),shape=ggplot2::guide_legend(reverse = TRUE),color=ggplot2::guide_legend(reverse=TRUE)) +
    ggplot2::labs(y="Density",x="Projection Score (Perm=Max)") +
    ggplot2::theme_minimal()

  t6 <- ggplot2::ggplot_build(p.diag6)
  minx <- min(t6$data[[1]]["x"])
  maxx <- max(t6$data[[1]]["x"])
  miny <- min(t6$data[[1]]["y"])
  maxy <- max(t6$data[[1]]["y"])

  p.diagt6 <- p.diag6 + ggplot2::geom_point(position = jitter, ggplot2::aes(y=maxy*0.25,shape=as.character(permy),color=as.character(permy)))



  ## Hides all but one legend and arranges ggplots into one plot ##
  if(plots=="all") {
    nol <- ggplot2::theme(legend.position = "hidden")
    lemon::grid_arrange_shared_legend(p.diagt4,gridExtra::arrangeGrob(p.diagt5+nol,p.diagt6+nol,ncol = 2),p.diagt1,ncol=1,nrow = 3,position = "top")
  }

  if(plots=="obs") {print(p.diagt4)}
  if(plots=="perm1") {print(p.diagt2)}
  if(plots=="perm2") {print(p.diagt3)}
  if(plots=="permdist") {print(p.diagt11)}
  if(plots=="min") {print(p.diagt5)}
  if(plots=="max") {print(p.diagt6)}

}




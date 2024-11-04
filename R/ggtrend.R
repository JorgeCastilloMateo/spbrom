#' @title Plot trend
#'
#' @description This function generates a ggplot graph with a weighted average
#' of indicators with an stationary decay of 1/t. Being t (time) the x-axis.
#'
#' @usage ggtrend(vol3D, yint, ychar, xchar, xnumbreaks, xlabbreaks)
#'
#' @param vol3D Input 3D volume.
#' @param yint y-axis limits.
#' @param ychar y-axis label.
#' @param xchar x-axis label.
#' @param xnumbreaks x-axis breaks.
#' @param xlabbreaks x-axis label breaks.
#'
#' @details The weighted trend will be summarized over the first dimension of the input 3D volume.
#'
#' @return A ggplot2 object depicting the t-times weighted trend.
#' 
#' @author Zeus Gracia-Tabuenca
#'
#' @examples
#' \donttest{
#' Y365 <- readRDS("data/Y365.rds")
#' ggtrend(vol3D = Y365,
#'         yint = c(0, 3),
#'         ychar = expression(t %*% hat(p)[t]),
#'         xchar = "t (year)",
#'         xnumbreaks = c(1, 21, 41, 61),
#'         xlabbreaks = c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"))
#' }
#'
#' @importFrom ggplot2 ggplot geom_hline geom_line geom_smooth theme_bw ylab xlab ylim scale_x_continuous
#' @export

ggtrend <- function(vol3D, yint, ychar, xchar, xnumbreaks, xlabbreaks){
  
  # Set records
  I365 <- apply(X = vol3D,
                MARGIN = c(2,3),
                FUN =  function(x) c(1,as.numeric(diff(cummax(x))>0)))
  I365[is.na(I365)] <- 0
  
  # Figure
  g.out <- ggplot(data = data.frame(p=apply(I365,1,mean), t=1:dim(I365)[1]),
                  mapping =  aes(x=t, y=p*t)) +
    geom_hline(yintercept=1,
               color = "gray") +
    geom_line() +
    geom_smooth(method = "loess",
                formula = "y ~ x",
                se=F,
                linetype = "dashed",
                color = "black") +
    theme_bw() +
    ylab(ychar) +
    xlab(xchar) +
    ylim(yint[1],yint[2]) +
    scale_x_continuous(breaks=xnumbreaks,
                       labels=xlabbreaks)
  
  # Output
  return(g.out)
}

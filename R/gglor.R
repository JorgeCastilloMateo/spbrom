#' @title Plot log odds ratio
#'
#' @description This function generates a ggplot graph with a log odds ratio
#' between a 3D volumes and its lag volume and conditioned by its lag-two. 
#' Computation is over the fist dimension (time) the x-axis.
#'
#' @usage function(vol3D, cond, corr, yint, ychar, xchar, xnumbreaks, xlabbreaks)
#'
#' @param vol3D Input 3D volume.
#' @param cond Condition of the lag-two volume. Options: marginal (.), positive (1), or negative (0).
#' @param corr Continuity correction for the odds ratio.
#' @param yint y-axis limits.
#' @param ychar y-axis label.
#' @param xchar x-axis label.
#' @param xnumbreaks x-axis breaks.
#' @param xlabbreaks x-axis label breaks.
#'
#' @details The log odds ratio is computed over the first dimension of the input 3D volume.
#'
#' @return A ggplot2 object depicting the log odds ratio over the first dimension.
#' 
#' @author Zeus Gracia-Tabuenca
#'
#' @examples
#' \donttest{
#' Y365 <- readRDS("data/Y365.rds")
#' gglor(vol3D = Y365,
#'       cond = ".",
#'       corr = 0.5,
#'       yint = c(1.1, 3.6),
#'       ychar = expression(LOR[t]),
#'       xchar = "t (year)",
#'       xnumbreaks = c(1, 21, 41, 61),
#'       xlabbreaks = c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"))
#' }
#'
#' @importFrom ggplot2 ggplot geom_line geom_smooth stat_smooth theme theme_bw ylab xlab ylim scale_x_continuous
#' @export

gglor <- function(vol3D, cond, corr, yint, ychar, xchar, xnumbreaks, xlabbreaks){
  
  # Set records
  I365 <- apply(X = vol3D,
                MARGIN = c(2,3),
                FUN =  function(x) c(1,as.numeric(diff(cummax(x))>0)))
  I365[is.na(I365)] <- 0
  
  # Apply lag1
  I365.lag1 <- lag3d(I365)
  # Apply lag2
  I365.lag2 <- lag3d(I365.lag1)
  
  # Concatenate volumes
  allvol <- simplify2array(list(I365, I365.lag1, I365.lag2))
  
  # Combinations
  comb01 <- expand.grid(rep(list(0:1),3))
  n.ls <- vector("list",nrow(comb01))
  names(n.ls) <- paste0(1:length(n.ls))
  for(ii in 1:nrow(comb01)){
    # Get vector
    vec01 <- unlist(comb01[ii,])
    # Variable name
    names(n.ls)[ii] <- paste0("n",paste0(vec01, collapse = ""))
    # Create variable
    auxvol <- apply(allvol, 1:3, function(x) as.numeric(all(x==vec01)))
    n.ls[[ii]] <- apply(auxvol, 1, function(x) sum(x, na.rm = T))
  }
  
  # Extract OR numerator/denominator based on 'cond'
  condvec <- as.character(cond)
  if(nchar(condvec)!=1) stop("'Cond' length is not correct. Revisit.")
  if(is.na(match(condvec, c(".","0","1")))) stop("Incorrect 'cond' input. Only '.01' characters are allowed.")
  
  # Add conditional based on case
  if(condvec=="."){
    n11 <- n.ls$n110 + n.ls$n111
    n10 <- n.ls$n100 + n.ls$n101
    n01 <- n.ls$n010 + n.ls$n011
    n00 <- n.ls$n000 + n.ls$n001
  }
  if(condvec=="0"){
    n11 <- n.ls$n110
    n10 <- n.ls$n100
    n01 <- n.ls$n010
    n00 <- n.ls$n000
  }
  if(condvec=="1"){
    n11 <- n.ls$n111
    n10 <- n.ls$n101
    n01 <- n.ls$n011
    n00 <- n.ls$n001
  }
  
  # Add continuity correction
  n11 <- n11 + corr
  n10 <- n10 + corr
  n01 <- n01 + corr
  n00 <- n00 + corr
  
  # Compute odds ratios by first dimension.
  or.df <- data.frame(OR=(n11*n00)/(n10*n01), t=1:dim(I365)[1])
  
  # Figure
  g.out <- ggplot(data = or.df[-1,],
                  mapping =  aes(x=t, y=log(OR))) +
    geom_line() +
    geom_smooth(method = "loess",
                formula = y ~ x,
                se = F,
                colour="black",
                linetype="dashed") +
    stat_smooth(method = "lm",
                formula = y~I(log(x-1)),
                se = F,
                colour="gray",
                linetype="solid") +
    theme_bw() +
    theme(legend.position = "none") +
    ylab(ychar) +
    xlab(xchar) +
    ylim(yint[1],yint[2]) +
    scale_x_continuous(breaks=xnumbreaks,
                       labels=xlabbreaks)
  
  # Output
  return(g.out)
}

lag3d <- function(vol3d){
  # Transpose
  vol3d_dim <- dim(vol3d)
  aux <- apply(vol3d, 3, function(x) t(x))
  # Apply lag1
  aux2 <- rbind(rep(0,ncol(aux)), aux[-nrow(aux),])
  # Transpose again
  aux3 <- array(data = apply(X = array(data = aux2, dim = vol3d_dim[c(2,1,3)]),
                             MARGIN = 3,
                             FUN = function(x) t(x)),
                dim = vol3d_dim)
  return(aux3)
}

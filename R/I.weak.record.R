#' @title Record Indicators
#' 
#' @description This function returns the record indicators of the values in a 
#'   vector. The record indicator for each value in a vector is a binary 
#'   variable that takes the value 1 if the corresponding value in the vector 
#'   is a record and 0 otherwise. Indicators for \eqn{r}-tied records have a 
#'   value \eqn{1 / r}.
#'  
#' @param X The vector
#' @return The binary vector of the same length as \code{X}, 
#'   indicating the record occurrence.
#'   
#' @author Jorge Castillo-Mateo
#' 
#' @export
I.weak.record <- function(X) {
  
  X[is.na(X)] <- -Inf
  
  T <- length(X)
  I <- c(1, rep(0, T - 1))
  count <- 1
  R <- X[1]
  
  for (t in 2:T) {
    if (X[t] > R) {
      count <- 1
      R <- X[t]
      I[t] <- 1
    } else if (X[t] == R) {
      count <- count + 1
      I[t] <- 1 / count
    }
  }
  
  return(I)
}

#' Convergence diagnostics for the models
#'
#' @description
#'   This function is used to check the convergence of the models.
#'
#' @param object A matrix of models (same model different chains) by columns or
#'   a list with each list a chain and each element of the chain several models
#'   by columns associated with different folds.
#' @return Diagnostics matrix.
#'
#' @author Jorge Castillo-Mateo
#'
#' @export
MCMCconvergence <- function(object) {

  if (is.matrix(object)) {

    n <- length(object[,1]$params)

    output <- matrix(nrow = 2, ncol = n,
                     dimnames = list(
                       "metric" = c("PSRF", "ESS"),
                       "params" = names(object[,1]$params)
                     )
    )

    for (i in 1:n) {
      if (!is.matrix(object[,1]$params[[i]]) ||
          ncol(object[,1]$params[[i]]) == 1) {
        output[1,i] <- coda::gelman.diag(coda::mcmc.list(
          coda::mcmc(object[,1]$params[[i]]),
          coda::mcmc(object[,2]$params[[i]])))$psrf[1]
        output[2,i] <-
          coda::effectiveSize(object[,1]$params[[i]]) +
          coda::effectiveSize(object[,2]$params[[i]])
      } else if (ncol(object[,1]$params[[i]]) <= 100) {
        const <- which(apply(object[,1]$params[[i]], 2, stats::sd) == 0 |
                         apply(object[,2]$params[[i]], 2, stats::sd) == 0)
        if (length(const) > 0) {
          object[,1]$params[[i]] <- object[,1]$params[[i]][,-const]
          object[,2]$params[[i]] <- object[,2]$params[[i]][,-const]
        }
        output[1,i] <- coda::gelman.diag(coda::mcmc.list(
          coda::mcmc(object[,1]$params[[i]]),
          coda::mcmc(object[,2]$params[[i]])))$mpsrf
        output[2,i] <- min(
          coda::effectiveSize(object[,1]$params[[i]]) +
          coda::effectiveSize(object[,2]$params[[i]]))
      } else {
        N <- ncol(object[,1]$params[[i]])
        output[1,i] <- 1
        output[2,i] <- 100000
        for (j in 1:N) {
          output[1,i] <- max(output[1,i],
            coda::gelman.diag(coda::mcmc.list(
              coda::mcmc(object[,1]$params[[i]][,j]),
              coda::mcmc(object[,2]$params[[i]][,j])))$psrf[1])
          output[2,i] <- min(output[2,i],
              coda::effectiveSize(object[,1]$params[[i]][,j]) +
              coda::effectiveSize(object[,2]$params[[i]][,j]))
        }
      }
    }

  } else {

    K <- ncol(object[[1]])
    n <- length(object[[1]][,1]$params)

    output <- array(dim = c(K, n, 2),
      dimnames = list(
        "fold" = 1:K,
        "params" = names(object[[1]][,1]$params),
        "metric" = c("PSRF", "ESS")
      )
    )

    for (k in 1:K) {
      for (i in 1:n) {
        if (!is.matrix(object[[1]][,k]$params[[i]]) ||
            ncol(object[[1]][,k]$params[[i]]) == 1) {
          output[k,i,1] <- coda::gelman.diag(coda::mcmc.list(
            coda::mcmc(object[[1]][,k]$params[[i]]),
            coda::mcmc(object[[2]][,k]$params[[i]])))$psrf[1]
          output[k,i,2] <-
            coda::effectiveSize(object[[1]][,k]$params[[i]]) +
            coda::effectiveSize(object[[2]][,k]$params[[i]])
        } else if (ncol(object[[1]][,k]$params[[i]]) <= 100) {
          const <- which(apply(object[[1]][,k]$params[[i]], 2, stats::sd) == 0 |
                         apply(object[[2]][,k]$params[[i]], 2, stats::sd) == 0)
          if (length(const) > 0) {
            object[[1]][,k]$params[[i]] <- object[[1]][,k]$params[[i]][,-const]
            object[[2]][,k]$params[[i]] <- object[[2]][,k]$params[[i]][,-const]
          }
          output[k,i,1] <- coda::gelman.diag(coda::mcmc.list(
            coda::mcmc(object[[1]][,k]$params[[i]]),
            coda::mcmc(object[[2]][,k]$params[[i]])))$mpsrf
          output[k,i,2] <- min(
            coda::effectiveSize(object[[1]][,k]$params[[i]]) +
            coda::effectiveSize(object[[2]][,k]$params[[i]]))
        } else {
          N <- ncol(object[[1]][,1]$params[[i]])
          output[k,i,1] <- 1
          output[k,i,2] <- 100000
          for (j in 1:N) {
            output[k,i,1] <- max(output[k,i,1],
              coda::gelman.diag(coda::mcmc.list(
                coda::mcmc(object[[1]][,k]$params[[i]][,j]),
                coda::mcmc(object[[2]][,k]$params[[i]][,j])))$psrf[1])
            output[k,i,2] <- min(output[k,i,2],
              coda::effectiveSize(object[[1]][,k]$params[[i]][,j]) +
              coda::effectiveSize(object[[2]][,k]$params[[i]][,j]))
          }
        }
      }
    }

  }

  return(output)
}

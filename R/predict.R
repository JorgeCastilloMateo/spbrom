#' Predict Method for Bivariate Record Occurrence Model Fits
#' @aliases predict.default predict.brom predict.rom predict
#'
#' @description
#'   Predicted values based on \code{brom} (bivariate record occurrence model)
#'   object.
#'
#' @details
#'   Predicted values are on the scale of the response (probabilities).
#'
#' @param object Object of class inheriting from "brom"
#' @param newdata An optional list with \code{x} in which to look for variables
#'   with which to predict, \code{coords} and \code{site} in which to predict,
#'   and \code{day} or \code{year}.
#'   If omitted, the fitted values are used.
#' @param newcoords,newextra.coords New coordinates for the sites in
#'   \code{newdata$site} or \code{newsite}.
#' @param newmatrix.mean.A New matrix containing the variables
#'   in the mean of A (if sp.A was TRUE).
#' @param newsite,newyear,newday A vector where for each observation includes
#'   the site, year or day.
#' @param type the type of prediction required. For \eqn{K}-fold
#'   cross-validation \code{"KFCV"} uses observed indicators from \code{data}
#'   and samples weak records from their prior probability.
#' @return A matrix where rows are MCMC posterior predictive samples and cols
#'   are observations.
#'
#' @rdname predict
#' @method predict brom
#' @export predict.brom
#' @export
predict.brom <- function(
    object,
    newdata,
    newcoords,
    newextra.coords,
    newmatrix.mean.A = matrix(1, nrow = nrow(newcoords)),
    newsite,
    newyear,
    newday,
    type = c("KFCV", "sample")) {

  type <- match.arg(type)

  if (missing(newsite)) newsite <- newdata$site
  if (missing(newyear)) newyear <- newdata$year
  if (missing(newday))  newday  <- newdata$day
  newsite <- match(newsite, sort(unique(newsite)))
  newyear <- match(newyear, sort(unique(newyear)))
  newday  <- match(newday,  sort(unique(newday)))
  SS <- length(unique(object$date$site))
  TT <- length(unique(object$date$year))
  LL <- length(unique(object$date$day))
  newSS <- length(unique(newsite))
  newNN <- nrow(newdata)

  tt <- stats::terms(object)
  Y <- stats::model.response(stats::model.frame(tt, newdata))
  Terms <- stats::delete.response(tt)
  m <- stats::model.frame(Terms, newdata)
  m0 <- m1 <- m
  m0[!(m0[,"lag.tx"] %in% 0:1),"lag.tx"] <- 0
  m0[!(m0[,"lag.tn"] %in% 0:1),"lag.tn"] <- 0
  m1[!(m1[,"lag.tx"] %in% 0:1),"lag.tx"] <- 1
  m1[!(m1[,"lag.tn"] %in% 0:1),"lag.tn"] <- 1
  X  <- stats::model.matrix(Terms, m )[,colnames(object$x)]
  X0 <- stats::model.matrix(Terms, m0)[,colnames(object$x)]
  X1 <- stats::model.matrix(Terms, m1)[,colnames(object$x)]

  if (type == "KFCV") {
    if (is.null(object$extra.coords)) {
      M <- 1
      dAux <- array(dist1(rbind(newcoords, object$coords)),
                    dim = c(newSS + SS, newSS + SS, 1))
    } else {
      M <- 1 + ncol(object$extra.coords)
      dAux <- array(dim = c(newSS + SS, newSS + SS, M))
      dAux[,,1] <- dist1(rbind(newcoords, object$coords))
      for (mm in 2:M) {
        dAux[,,mm] <- dist1(rbind(
          newextra.coords[,mm-1, drop = FALSE],
          object$extra.coords[,mm-1, drop = FALSE]
        ))
      }
    }
  }

  indLag1  <- grep("lag.tx", colnames(X))
  indLag2  <- grep("lag.tn", colnames(X))

  if (object$scale.cov) {
    X  <- sweep(X,  2, attr(object$x, "scaled:center"), FUN = '-')
    X0 <- sweep(X0, 2, attr(object$x, "scaled:center"), FUN = '-')
    X1 <- sweep(X1, 2, attr(object$x, "scaled:center"), FUN = '-')
    X  <- sweep(X,  2, attr(object$x, "scaled:scale"), FUN = '/')
    X0 <- sweep(X0, 2, attr(object$x, "scaled:scale"), FUN = '/')
    X1 <- sweep(X1, 2, attr(object$x, "scaled:scale"), FUN = '/')

    if ("(Intercept)" %in% colnames(X)) {
      X[,"(Intercept)"] <- 1
      X0[,"(Intercept)"] <- 1
      X1[,"(Intercept)"] <- 1
    }
  }

  BB <- nrow(object$params$beta1)

  weak1 <- which(!(m[,"lag.tx"] %in% 0:1))
  weak2 <- which(!(m[,"lag.tn"] %in% 0:1))
  Nweak <- c(length(weak1), length(weak2))
  prob1 <- m[weak1,"lag.tx"]
  prob2 <- m[weak2,"lag.tn"]

  if (type == "KFCV") {

    if (!object$sp.A) {
      object$params$hpA <- as.matrix(NA)
      object$x.A <- as.matrix(NA)
    }

    Xb <- predictBglm(
      BB, Nweak, newNN,
      SS, newSS, TT, LL,
      X, X0, X1,
      object$date$year - 1, object$date$day - 1,
      newsite - 1, newyear - 1, newday - 1,
      dAux,
      object$params$beta1, object$params$beta2,
      object$params$w1tls, object$params$w2tls,
      object$params$a, object$sp.A, object$params$hpA,
      object$x.A, newmatrix.mean.A, ncol(newmatrix.mean.A),
      as.matrix(object$params$decay), M,
      weak1 - 1, weak2 - 1,
      indLag1 - 1, indLag2 - 1,
      prob1, prob2)

  } else {

    Xb <- predGlmBerKFCV(
      BB, Nweak, newNN,
      X, X0, X1,
      object$params$beta1, object$params$beta2,
      weak1 - 1, weak2 - 1,
      indLag1 - 1, indLag2 - 1,
      prob1, prob2)

    if (object$sp.A) {
      for (ii in 1:SS) {
        s <- 1:(TT*LL) + (ii-1) * (TT*LL)
        Xb[,s,1] <- Xb[,s,1] +
          object$params$a[,ii] * object$params$w1tls[,s]
        Xb[,s,2] <- Xb[,s,2] +
          object$params$a[,ii + 2*SS] * object$params$w1tls[,s] +
          object$params$a[,ii + SS] * object$params$w2tls[,s]
      }
    } else {
      Xb[,,1] <- Xb[,,1] +
        object$params$a[,1] * object$params$w1tls
      Xb[,,2] <- Xb[,,2] +
        object$params$a[,3] * object$params$w1tls +
        object$params$a[,2] * object$params$w2tls
    }
  }

  Y <- stats::pnorm(Xb)

  return(Y)
}



#' @rdname predict
#' @method predict rom
#' @export predict.rom
#' @export
predict.rom <- function(
    object,
    newdata,
    newcoords,
    newextra.coords,
    newsite,
    newyear,
    newday,
    type = c("KFCV", "sample")) {

  type <- match.arg(type)

  if (missing(newsite)) newsite <- newdata$site
  if (missing(newyear)) newyear <- newdata$year
  if (missing(newday))  newday  <- newdata$day
  newsite <- match(newsite, sort(unique(newsite)))
  newyear <- match(newyear, sort(unique(newyear)))
  newday  <- match(newday,  sort(unique(newday)))
  SS <- length(unique(object$date$site))
  TT <- length(unique(object$date$year))
  LL <- length(unique(object$date$day))
  newSS <- length(unique(newsite))
  newNN <- nrow(newdata)

  tt <- stats::terms(object)
  Y <- stats::model.response(stats::model.frame(tt, newdata))
  Terms <- stats::delete.response(tt)
  m <- stats::model.frame(Terms, newdata)
  m0 <- m1 <- m
  m0[!(m0[,"lag"] %in% 0:1),"lag"] <- 0
  m1[!(m1[,"lag"] %in% 0:1),"lag"] <- 1
  X  <- stats::model.matrix(Terms, m )[,colnames(object$x)]
  X0 <- stats::model.matrix(Terms, m0)[,colnames(object$x)]
  X1 <- stats::model.matrix(Terms, m1)[,colnames(object$x)]

  if (type == "KFCV") {
    if (is.null(object$extra.coords)) {
      M <- 1
      dAux <- array(dist1(rbind(newcoords, object$coords)),
                    dim = c(newSS + SS, newSS + SS, 1))
    } else {
      M <- 1 + ncol(object$extra.coords)
      dAux <- array(dim = c(newSS + SS, newSS + SS, M))
      dAux[,,1] <- dist1(rbind(newcoords, object$coords))
      for (mm in 2:M) {
        dAux[,,mm] <- dist1(rbind(
          newextra.coords[,mm-1, drop = FALSE],
          object$extra.coords[,mm-1, drop = FALSE]
        ))
      }
    }
  }

  indLag  <- grep("lag", colnames(X))

  if (object$scale.cov) {
    X  <- sweep(X,  2, attr(object$x, "scaled:center"), FUN = '-')
    X0 <- sweep(X0, 2, attr(object$x, "scaled:center"), FUN = '-')
    X1 <- sweep(X1, 2, attr(object$x, "scaled:center"), FUN = '-')
    X  <- sweep(X,  2, attr(object$x, "scaled:scale"), FUN = '/')
    X0 <- sweep(X0, 2, attr(object$x, "scaled:scale"), FUN = '/')
    X1 <- sweep(X1, 2, attr(object$x, "scaled:scale"), FUN = '/')

    if ("(Intercept)" %in% colnames(X)) {
      X[,"(Intercept)"] <- 1
      X0[,"(Intercept)"] <- 1
      X1[,"(Intercept)"] <- 1
    }
  }

  BB <- nrow(object$params$beta)

  weak <- which(!(m[,"lag"] %in% 0:1))
  Nweak <- c(length(weak), 0)
  prob <- m[weak,"lag"]

  if (type == "KFCV") {

    if (is.null(object$params$hpA))
      object$params$hpA <- as.matrix(NA)

    Xb <- predictBglm(
      BB, Nweak, newNN,
      SS, newSS, TT, LL,
      X, X0, X1,
      object$date$year - 1, object$date$day - 1,
      newsite - 1, newyear - 1, newday - 1,
      dAux,
      object$params$beta, object$params$beta,
      object$params$wtls, object$params$wtls,
      cbind(object$params$a, object$params$a, object$params$a),
      object$sp.A, object$params$hpA,
      as.matrix(NA), as.matrix(NA), 0,
      as.matrix(object$params$decay), M,
      weak - 1, weak - 1,
      indLag - 1, indLag - 1,
      prob, prob)[,,1]

  } else {

    Xb <- predGlmBerKFCV(
      BB, Nweak, newNN,
      X, X0, X1,
      object$params$beta, object$params$beta,
      weak - 1, weak - 1,
      indLag - 1, indLag - 1,
      prob, prob)[,,1]

    if (object$sp.A) {
      for (ii in 1:SS) {
        s <- 1:(TT*LL) + (ii-1) * (TT*LL)
        Xb[,s] <- Xb[,s] +
          object$params$a[,ii] * object$params$wtls[,s]
      }
    } else {
      Xb <- Xb + object$params$a * object$params$wtls
    }
  }

  Y <- stats::pnorm(Xb)

  return(Y)
}

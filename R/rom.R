#' Univariate Record Occurrence Model Fitting
#'
#' @description
#'   This function fits the Univariate Record Occurrence Model in a Bayesian
#'   framework using MCMC.
#'
#' @details
#'   The fitting algorithm is a data-augmentation Metropolis-within-Gibbs MCMC
#'   algorithm.
#'
#'   \code{formula} requires one variable, \code{lag}, that is the record
#'   indicator of the previous day.
#'
#'   \code{data} requires column names \code{trend} (annual trend), and
#'   \code{lag}.
#'
#'   The response variable must be \code{Y}. It can be 0's and 1's indicating
#'   the occurrence of a record, but values \eqn{1/r} (\eqn{r=2,3,\ldots}) are
#'   allowed indicating the probability of a tied record being a strong record.
#'   (If all values are 0's and 1's the function could return an error.)
#'
#' @note It is necessary to have data ordered as using if c() over a
#'   3-dimensional array where the first dimension is year, the second is day,
#'   and the third is site. This is, first all data for the first site, then
#'   for the first day and it progresses by years. (This will also be necessary
#'   for prediction.)
#'
#' @param formula an object of class "formula": a symbolic description of the
#'   model to be fitted. (Has some requirements with \code{trend}, and
#'   \code{lag.tn}).
#' @param data a data frame (\eqn{N \times p}) containing the variables in the
#'   model.
#' @param coords a matrix or data frame (\eqn{n \times 2}) containing the
#'   coordinates of the sites.
#' @param extra.coords a matrix or data frame (\eqn{n \times k}) containing
#'   spatial covariates to be included in the spatial correlation function.
#' @param sp.A a boolean indicating a spatially-varying coregionalization.
#' @param site,year,day vectors (\eqn{N \times 1}) containing the site, year,
#'   and day the measurements were taken.
#' @param scale.cov a boolean indicating whether or not the columns of the
#'   design matrix must be scaled to have zero mean and unit variance.
#' @param decay.prior a character string indicating the prior for the decay
#'   parameter, gamma, uniform in an interval, or fixed.
#' @param inits a vector of initial values (if NULL, random values).
#' @param prior a list containing the parameters of the prior distributions:
#'   the precision for beta; the scale for a; for the main decay if fixed the
#'   fixed values, 1 + 1 for each column in \code{extra.coords}, if gamma(a, b)
#'   first all a's then all b's (again 1 + 1 each), if unif(a, b) the same;
#'   both hyperparameters of the gamma priors, the fixed decay for the
#'   spatially-varying a's.
#' @param n.sims,n.thin,n.burnin,n.report
#'   (i) Number of iterations not discarded.
#'   (ii) Thinning rate.
#'   (iii) Number of iterations discarded at the beginning.
#'   (iv) Report the number of iterations rate.
#' @return A \code{rom} list with many elements. Samples from the model
#'   parameters are in \code{params} where rows are MCMC simulations and cols
#'   are parameters.
#' @examples
#' # rom(tx ~ trend + lag, data, coords)
#'
#' @export
rom <- function(
    formula,
    data,
    coords,
    extra.coords,
    sp.A = FALSE,
    site,
    year,
    day,
    scale.cov = TRUE,
    decay.prior = c("fixed", "unif", "gamma"),
    inits = NULL,
    prior = list("beta" = 1e-02,
                 "a" = 1 / 5^2,
                 "decay" = c(2, 100),
                 "sigma" = c(2, 1),
                 "decay_A" = 0.01),
    n.sims = 1000,
    n.thin = 1,
    n.burnin = 1000,
    n.report = 100
) {

  if (missing(site)) site <- data$site
  if (missing(year)) year <- data$year
  if (missing(day))  day  <- data$day
  site <- match(site, sort(unique(site)))
  year <- match(year, sort(unique(year)))
  day  <- match(day,  sort(unique(day)))
  SS <- length(unique(site))
  TT <- length(unique(year))
  LL <- length(unique(day))

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- al <- mf[c(1L, m)]
  al[[1L]] <- quote(stats::alias)
  names(al)[2] <- "object"
  al <- eval(al, parent.frame())
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")

  mf0 <- mf1 <- mf
  mf0[!(mf0[,"lag"] %in% 0:1),"lag"] <- 0
  mf1[!(mf1[,"lag"] %in% 0:1),"lag"] <- 1

  if (is.null(attr(al$Complete, "dimnames"))) {
    X  <- model.matrix(mt, mf )
    X0 <- model.matrix(mt, mf0)
    X1 <- model.matrix(mt, mf1)
  } else {
    X  <- model.matrix(mt, mf )[, attr(al$Complete, "dimnames")[[2]]]
    X0 <- model.matrix(mt, mf0)[, attr(al$Complete, "dimnames")[[2]]]
    X1 <- model.matrix(mt, mf1)[, attr(al$Complete, "dimnames")[[2]]]
  }

  decay.prior <- match.arg(decay.prior)
  if (missing(extra.coords)) {
    M <- 1
    dAux <- array(dist1(coords), dim = c(SS, SS, 1))
    drange <- as.matrix(range(dAux[dAux != 0]))
  } else {
    M <- 1 + ncol(extra.coords)
    dAux <- array(dim = c(SS, SS, M))
    dAux[,,1] <- dist1(coords)
    for (mm in 2:M) {
      dAux[,,mm] <- dist1(extra.coords[,mm-1, drop = FALSE])
    }
    drange <- apply(dAux, 3, function(x) range(x[x != 0]))
  }

  if (decay.prior == "fixed") {
    decay.prior <- 1
  } else if (decay.prior == "unif") {
    decay.prior <- 2
  } else { # decay.prior == "gamma"
    decay.prior <- 3
  }

  indLag  <- grep("lag", colnames(X))

  if (scale.cov) {
    X <- scale(X)
    X0 <- sweep(X0, 2, attr(X, "scaled:center"), FUN = '-')
    X1 <- sweep(X1, 2, attr(X, "scaled:center"), FUN = '-')
    X0 <- sweep(X0, 2, attr(X, "scaled:scale"), FUN = '/')
    X1 <- sweep(X1, 2, attr(X, "scaled:scale"), FUN = '/')

    if ("(Intercept)" %in% colnames(X)) {
      X[,"(Intercept)"] <- 1
      X0[,"(Intercept)"] <- 1
      X1[,"(Intercept)"] <- 1
    }
  }

  k <- ncol(X)
  N <- nrow(X)
  abLim <- rep(-1, N)
  abLim[Y == 1] <- 1
  weak <- which(!(Y %in% 0:1))
  Nweak <- rep(NA, 2)
  Nweak[1] <- length(weak)

  day1 <- day[weak] + 1
  year1 <- year[weak]

  Weak <- match(apply(apply(cbind(day1, year1, site[weak]), 2, format, width = 3), 1, paste, collapse = ""),
                apply(apply(cbind(day, year, site), 2, format, width = 3), 1, paste, collapse = ""))

  weakl0 <- which(day == min(day) & !(mf[,"lag"] %in% 0:1))
  probl0 <- mf[weakl0,"lag"]
  Nweak[2] <- length(weakl0)

  ### MODEL
  if (sp.A) {
    stop("Model not implemented.")
    # beta + wtls + a11s + mu prec decay 11 + decay
    keep <- matrix(nrow = n.sims / n.thin, ncol = k + N + SS + 3 + M)
    if (is.null(inits)) {
      inits <- rnorm(k + N + SS + 3 + M)
      inits[1:SS + k + N] <- exp(inits[1:SS + k + N])
      inits[2 + k + N + SS] <- exp(inits[2 + k + N + SS])
      inits[3 + k + N + SS] <- prior$decay_A
      if (decay.prior == 1) {
        inits[1:M + k + N + SS + 3] <- prior$decay[1:M]
      } else if (decay.prior == 2) {
        inits[1:M + k + N + SS + 3] <- runif(M, prior$decay[1], prior$decay[2])
      } else {
        inits[1:M + k + N + SS + 3] <- runif(M, 3 / drange[2,], 3 / drange[1,])
      }
    } else if ((k + N + SS + 3 + M) != length(inits)) {
      stop("'inits' does not have the proper length")
    }

    ## code must be udated, build function uglm2 from bglm2
    #keep <- uglm2(
    #  N, k, Nweak, X, X0, X1,
    #  cbind(inits[1:k], inits[1:k + k]),
    #  cbind(inits[1:N + 2 * k], inits[1:N + 2 * k + N]),
    #  inits[1:SS + 2 * k + 2 * N],
    #  inits[1:SS + 2 * k + 2 * N + SS],
    #  inits[1:SS + 2 * k + 2 * N + 2 * SS],
    #  inits[1:9 + 2 * k + 2 * N + 3 * SS],
    #  inits[1:M + 2 * k + 2 * N + 3 * SS + 9], M,
    #  dAux, decay.prior,
    #  SS, TT, LL,
    #  site - 1, year - 1, day - 1,
    #  weak1 - 1, weak2 - 1,
    #  Weak1 - 1, Weak2 - 1,
    #  !is.na(Weak1), !is.na(Weak2),
    #  indLag1 - 1, indLag2 - 1,
    #  Y, abLim,
    #  weak1l0 - 1, weak2l0 - 1,
    #  prob1l0, prob2l0,
    #  prior$beta * diag(k),
    #  0, prior$beta,
    #  prior$sigma[1], prior$sigma[2],
    #  prior$decay[1:M], prior$decay[1:M + M],
    #  keep,
    #  n.sims, n.thin, n.burnin, n.report
    #)

    colnames(keep) <- c(
      colnames(X),
      paste0("wt", 1:TT, "l", rep(1:LL, each = TT), "s", rep(1:SS, each = TT * LL)),
      paste0("a11s", 1:SS), "la11", "preca11", "decaya11",
      paste0("decay", 1:M))

    keep.list <- list()
    keep.list$params$beta  <- keep[,1:k]
    keep.list$params$wtls  <- keep[,1:N + k]
    keep.list$params$a     <- keep[,1:SS + k + N]
    keep.list$params$hpA   <- keep[,1:3 + k + N + SS]
    keep.list$params$decay <- keep[,1:M + k + N + SS + 3]
    keep.list$date$day  <- day
    keep.list$date$year <- year
    keep.list$date$site <- site
    keep.list$coords <- coords
    if (!missing(extra.coords))
      keep.list$extra.coords <- extra.coords
    keep.list$sp.A <- sp.A
    keep.list$y <- Y
    keep.list$x <- X
    keep.list$call <- cl
    keep.list$terms <- mt
    keep.list$model <- mf
    keep.list$scale.cov <- scale.cov
  } else {
    # beta + wtls + a11 + decay
    keep <- matrix(nrow = n.sims / n.thin, ncol = k + N + 1 + M)
    if (is.null(inits)) {
      inits <- rnorm(k + N + 1 + M)
      inits[1 + k + N] <- abs(inits[1 + k + N])
      if (decay.prior == 1) {
        inits[1:M + k + N + 1] <- prior$decay[1:M]
      } else if (decay.prior == 2) {
        inits[1:M + k + N + 1] <- runif(M, prior$decay[1:M], prior$decay[1:M + M])
      } else {
        inits[1:M + k + N + 1] <- runif(M, 3 / drange[2,], 3 / drange[1,])
      }
    } else if ((k + N + 1 + M) != length(inits)) {
      stop("'inits' does not have the proper length")
    }

    keep <- uglm(
      N, k, Nweak, X, X0, X1,
      inits[1:k],
      inits[1:N + k],
      inits[1 + k + N],
      inits[1:M + k + N + 1], M,
      dAux, decay.prior,
      SS, TT, LL,
      site - 1, year - 1, day - 1,
      weak - 1, Weak - 1,
      !is.na(Weak), indLag - 1,
      Y, abLim,
      weakl0 - 1, probl0,
      prior$beta * diag(k),
      0, prior$beta,
      prior$a,
      prior$decay[1:M], prior$decay[1:M + M],
      keep,
      n.sims, n.thin, n.burnin, n.report
    )

    colnames(keep) <- c(
      colnames(X),
      paste0("wt", 1:TT, "l", rep(1:LL, each = TT), "s", rep(1:SS, each = TT * LL)),
      "a11", paste0("decay", 1:M))

    keep.list <- list()
    keep.list$params$beta  <- keep[,1:k]
    keep.list$params$wtls  <- keep[,1:N + k]
    keep.list$params$a     <- keep[,1 + k + N]
    keep.list$params$decay <- keep[,1:M + k + N + 1]
    keep.list$date$day  <- day
    keep.list$date$year <- year
    keep.list$date$site <- site
    keep.list$coords <- coords
    if (!missing(extra.coords))
      keep.list$extra.coords <- extra.coords
    keep.list$sp.A <- sp.A
    keep.list$y <- Y
    keep.list$x <- X
    keep.list$call <- cl
    keep.list$terms <- mt
    keep.list$model <- mf
    keep.list$scale.cov <- scale.cov
  }

  keep.list <- structure(keep.list, class = "rom")

  return(keep.list)
}

#' Metrics for Validation
#'
#' @description
#'   This function is used to validate the predicted probabilities of record
#'   with the actual values.
#'
#' @details
#'   Values different from 0 or 1 (tied records) are removed from the metrics.
#'
#' @param Y The observed values (vector)
#' @param p The estimated probabilities (matrix with replicates by rows)
#' @param year The year associated with the values (vector)
#' @param metric A character string \code{"DIC"}, \code{"CV"} or \code{"pCV"}
#'   (CV is the posterior mean of the metric and pCV uses the posterior mean
#'   of the probabilities of record)
#' @param D1,D2 Two shorter year periods on which to calculate the metrics, in
#'   addition to the period of all years
#' @return If \code{metric = "CV"}; the metrics (vector):
#'   \item{AUC}{Area under the ROC curve}
#'   \item{AUC1}{AUC in \code{D1}}
#'   \item{AUC2}{AUC in \code{D2}}
#'   \item{J}{Jaccard index}
#'   \item{J1}{J in \code{D1}}
#'   \item{J2}{J in \code{D2}}
#'   \item{AUPRC}{Area under the Precission-Recall curve}
#'   \item{AUPRC1}{AUPRC in \code{D1}}
#'   \item{AUPRC2}{AUPRC in \code{D2}}
#'   \item{BS}{Brier Score}
#'   \item{BS1}{BS in \code{D1}}
#'   \item{BS2}{BS in \code{D2}}
#'
#'   If \code{metric = "DIC"}; a vector with DIC, D (deviance), and penalty (pD)
#'
#' @author Jorge Castillo-Mateo
#'
#' @export
metrics <- function(Y, p, year,
                    metric = c("DIC", "CV", "pCV"),
                    D1 = 2:31, D2 = 32:64) {

  metric <- match.arg(metric)
  weak <- !(Y %in% 0:1)
  Y    <- Y[!weak]
  p    <- p[,!weak]

  if (metric == "DIC") {

    p <- t(p)
    # D <- -2 * mean(colSums(Y * log(p) + (1 - Y) * log(1 - p)))
    D <- -2 * mean(colSums(log(abs((1 - Y) - p))))

    p <- rowMeans(p)
    pD <- D + 2 * sum(log(abs((1 - Y) - p)))

    return(c("DIC" = D + pD, "D" = D, "pD" = pD))

  } else if (metric == "CV") {

    year <- year[!weak]

    metric <- rep(NA, 12)
    names(metric) <- c(
      "AUC", "AUC1", "AUC2",
      "J", "J1", "J2",
      "AUPRC", "AUPRC1", "AUPRC2",
      "BS", "BS1", "BS2")

    # AUC
    roc <- function(x, Y) suppressMessages(pROC::auc(Y, x))
    metric[1] <- mean(apply(p, 1, roc, Y = Y))
    metric[2] <- mean(apply(p[,year %in% D1], 1, roc, Y = Y[year %in% D1]))
    metric[3] <- mean(apply(p[,year %in% D2], 1, roc, Y = Y[year %in% D2]))

    # Jaccard index
    A1 <- Y == 1
    A0 <- Y == 0
    TP <- rowSums(p[,A1])
    FP <- rowSums(p[,A0])
    FN <- rowSums(1 - p[,A1])
    metric[4] <- mean(TP / (TP + FP + FN))
    TP <- rowSums(p[,A1 & year %in% D1])
    FP <- rowSums(p[,A0 & year %in% D1])
    FN <- rowSums(1 - p[,A1 & year %in% D1])
    metric[5] <- mean(TP / (TP + FP + FN))
    TP <- rowSums(p[,A1 & year %in% D2])
    FP <- rowSums(p[,A0 & year %in% D2])
    FN <- rowSums(1 - p[,A1 & year %in% D2])
    metric[6] <- mean(TP / (TP + FP + FN))

    # AUPRC
    prroc <- function(x, Y) PRROC::pr.curve(scores.class0 = x, weights.class0 = Y)$auc.integral
    metric[7] <- mean(apply(p, 1, prroc, Y = Y))
    metric[8] <- mean(apply(p[,year %in% D1], 1, prroc, Y = Y[year %in% D1]))
    metric[9] <- mean(apply(p[,year %in% D2], 1, prroc, Y = Y[year %in% D2]))

    # Brier Score
    bs <- function(x, Y) mean((Y - x)^2)
    metric[10] <- mean(apply(p, 1, bs, Y = Y))
    metric[11] <- mean(apply(p[,year %in% D1], 1, bs, Y = Y[year %in% D1]))
    metric[12] <- mean(apply(p[,year %in% D2], 1, bs, Y = Y[year %in% D2]))

    return(metric)

  } else {

    year <- year[!weak]

    metric <- rep(NA, 12)
    names(metric) <- c(
      "AUC", "AUC1", "AUC2",
      "J", "J1", "J2",
      "AUPRC", "AUPRC1", "AUPRC2",
      "BS", "BS1", "BS2")

    p <- colMeans(p)

    # AUC
    metric[1] <- suppressMessages(pROC::auc(Y, p))
    metric[2] <- suppressMessages(pROC::auc(Y[year %in% D1], p[year %in% D1]))
    metric[3] <- suppressMessages(pROC::auc(Y[year %in% D2], p[year %in% D2]))

    # Jaccard index
    A1 <- Y == 1
    A0 <- Y == 0
    TP <- sum(p[A1])
    FP <- sum(p[A0])
    FN <- sum(1 - p[A1])
    metric[4] <- TP / (TP + FP + FN)
    TP <- sum(p[A1 & year %in% D1])
    FP <- sum(p[A0 & year %in% D1])
    FN <- sum(1 - p[A1 & year %in% D1])
    metric[5] <- TP / (TP + FP + FN)
    TP <- sum(p[A1 & year %in% D2])
    FP <- sum(p[A0 & year %in% D2])
    FN <- sum(1 - p[A1 & year %in% D2])
    metric[6] <- TP / (TP + FP + FN)

    # AUPRC
    metric[7] <- PRROC::pr.curve(scores.class0 = p, weights.class0 = Y)$auc.integral
    metric[8] <- PRROC::pr.curve(scores.class0 = p[year %in% D1],
                                 weights.class0 = Y[year %in% D1])$auc.integral
    metric[9] <- PRROC::pr.curve(scores.class0 = p[year %in% D2],
                                 weights.class0 = Y[year %in% D2])$auc.integral

    # Brier Score
    metric[10] <- mean((Y - p)^2)
    metric[11] <- mean((Y[year %in% D1] - p[year %in% D1])^2)
    metric[12] <- mean((Y[year %in% D2] - p[year %in% D2])^2)

    return(metric)
  }
}

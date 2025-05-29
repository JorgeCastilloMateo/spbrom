##################################
### Section 4 - MODEL CHECKING ###
##################################

# INPUT: pp is the array of probabilities
#   obtained in 008_DIC.R (any model)

### Posterior predictive distribution
I <- array(dim = c(1000, 231840, 2))
I[,,1] <- matrix(rbinom(1000 * 231840, 1, pp[,,1]), 1000, 231840)
I[,,2] <- matrix(rbinom(1000 * 231840, 1, pp[,,2]), 1000, 231840)
Ipre <- array(I, dim = c(1000, 63, 92, 40, 2))
Iobs <- array(M5[,1]$y, dim = c(63, 92, 40, 2))

# Preparation
# Months
month <- list(
  jun = 1:30,
  jul = 31:61,
  aug = 62:92)

# Features
# Nall
Nallobs <- matrix(nrow =  3 * 40, ncol = 3)
Nallpre <- array(dim = c(3 * 40, 1000, 3))
j <- 0; t <- 1:63
for (i in 1:40) {
  for (m in 1:3) {
    j <- j + 1
    l <- month[[m]]
    L <- length(l)
    Nallobs[j, 1]   <- sum(Iobs[t, l, i, 1]) / L + 1
    Nallobs[j, 2]   <- sum(Iobs[t, l, i, 2]) / L + 1
    Nallobs[j, 3]   <- sum(Iobs[t, l, i, 1] * Iobs[t, l, i, 2]) / L + 1
    Nallpre[j, , 1] <- apply(Ipre[, t, l, i, 1], 1, sum) / L + 1
    Nallpre[j, , 2] <- apply(Ipre[, t, l, i, 2], 1, sum) / L + 1
    Nallpre[j, , 3] <- apply(Ipre[, t, l, i, 1] * Ipre[, t, l, i, 2], 1, sum) / L + 1
  }
}

# N21st
N21stobs <- matrix(nrow =  3 * 40, ncol = 3)
N21stpre <- array(dim = c(3 * 40, 1000, 3))
j <- 0; t <- 41:63
for (i in 1:40) {
  for (m in 1:3) {
    j <- j + 1
    l <- month[[m]]
    L <- length(l)
    N21stobs[j, 1]   <- sum(Iobs[t, l, i, 1]) / L
    N21stobs[j, 2]   <- sum(Iobs[t, l, i, 2]) / L
    N21stobs[j, 3]   <- sum(Iobs[t, l, i, 1] * Iobs[t, l, i, 2]) / L
    N21stpre[j, , 1] <- apply(Ipre[, t, l, i, 1], 1, sum) / L
    N21stpre[j, , 2] <- apply(Ipre[, t, l, i, 2], 1, sum) / L
    N21stpre[j, , 3] <- apply(Ipre[, t, l, i, 1] * Ipre[, t, l, i, 2], 1, sum) / L
  }
}

# R
E <- sum(1 / 55:64)
Robs <- matrix(nrow =  3 * 40, ncol = 3)
Rpre <- array(dim = c(3 * 40, 1000, 3))
j <- 0; t <- 54:63
for (i in 1:40) {
  for (m in 1:3) {
    j <- j + 1
    l <- month[[m]]
    L <- length(l)
    Robs[j, 1]   <- sum(Iobs[t, l, i, 1]) / (L * E)
    Robs[j, 2]   <- sum(Iobs[t, l, i, 2]) / (L * E)
    Robs[j, 3]   <- sum(Iobs[t, l, i, 1] * Iobs[t, l, i, 2]) / L
    Rpre[j, , 1] <- apply(Ipre[, t, l, i, 1], 1, sum) / (L * E)
    Rpre[j, , 2] <- apply(Ipre[, t, l, i, 2], 1, sum) / (L * E)
    Rpre[j, , 3] <- apply(Ipre[, t, l, i, 1] * Ipre[, t, l, i, 2], 1, sum) / L
  }
}

# ERS
ERSobs <- matrix(nrow =  3 * 63, ncol = 3)
ERSpre <- array(dim = c(3 * 63, 1000, 3))
j <- 0
for (t in 1:63) {
  for (m in 1:3) {
    j <- j + 1
    l <- month[[m]]
    L <- length(l)
    ERSobs[j, 1]   <- (t + 1) * sum(Iobs[t, l, , 1]) / (40 * L)
    ERSobs[j, 2]   <- (t + 1) * sum(Iobs[t, l, , 2]) / (40 * L)
    ERSobs[j, 3]   <- (t + 1) * sum(Iobs[t, l, , 1] * Iobs[t, l, , 2]) / (40 * L)
    ERSpre[j, , 1] <- (t + 1) * apply(Ipre[, t, l, , 1], 1, sum) / (40 * L)
    ERSpre[j, , 2] <- (t + 1) * apply(Ipre[, t, l, , 2], 1, sum) / (40 * L)
    ERSpre[j, , 3] <- (t + 1) * apply(Ipre[, t, l, , 1] * Ipre[, t, l, , 2], 1, sum) / (40 * L)
  }
}

# PIT
PIThist <- function(I, p, J = 10,
                    title = NULL,
                    picture.name = "photo.pdf",
                    save = TRUE) {

  func <- function(u, I, p) {

    p0 <- rowMeans(I <  p)
    p1 <- rowMeans(I <= p)

    ifelse(u <= p0, 0, ifelse(u > p1, 1, (u - p0) / (p1 - p0)))
  }

  f <- rep(0, J); f[J] <- 1
  for (j in 1:(J - 1)) {
    f[j] <- mean(func(j / J, I, p))
  }

  f <- diff(c(0, f))

  df <- data.frame("PIT" = 1:J / J, "Relative Frequency" = J * f)

  gg <- ggplot2::ggplot(df, ggplot2::aes(PIT, Relative.Frequency)) +
    ggplot2::geom_bar(stat = 'identity', just = 1, width = 1 / J, col = "black", fill = "white") +
    ggplot2::theme_classic() +
    ggplot2::labs(x="PIT", y="Relative Frequency", title=title) +
    ggplot2::ylim(c(0, 1.75))

  if (save) {
    ggplot2::ggsave(picture.name, gg, width = 4, height = 4)
  }

  gg
}

i <- 1
PIThist(Nallobs[,i], Nallpre[,,i],
        title = bquote(bar(N)["64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_Nall_max.pdf")
PIThist(N21stobs[,i], N21stpre[,,i],
        title = bquote(bar(N)["42:64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_N21st_max.pdf")
PIThist(Robs[,i], Rpre[,,i],
        title = bquote(R["55:64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_R_max.pdf")
PIThist(ERSobs[,i], ERSpre[,,i],
        title = expression(t %*% bar(ERS)["t,month"](D)),
        picture.name = "inst/img/SUPP_PIT_ERS_max.pdf")

i <- 2
PIThist(Nallobs[,i], Nallpre[,,i],
        title = bquote(underline(N)["64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_Nall_min.pdf")
PIThist(N21stobs[,i], N21stpre[,,i],
        title = bquote(underline(N)["42:64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_N21st_min.pdf")
PIThist(Robs[,i], Rpre[,,i],
        title = bquote(underline(R)["55:64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_R_min.pdf")
PIThist(ERSobs[,i], ERSpre[,,i],
        title = expression(t %*% underline(ERS)["t,month"](D)),
        picture.name = "inst/img/SUPP_PIT_ERS_min.pdf")

i <- 3
PIThist(Nallobs[,i], Nallpre[,,i],
        title = bquote(bar(underline(N))["64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_Nall_joint.pdf")
PIThist(N21stobs[,i], N21stpre[,,i],
        title = bquote(bar(underline(N))["42:64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_N21st_joint.pdf")
PIThist(Robs[,i], Rpre[,,i],
        title = bquote(bar(underline(R))["55:64,month"](s)),
        picture.name = "inst/img/SUPP_PIT_R_joint.pdf")
PIThist(ERSobs[,i], ERSpre[,,i],
        title = expression(t %*% bar(underline(ERS))["t,month"](D)),
        picture.name = "inst/img/SUPP_PIT_ERS_joint.pdf")



# Observed vs predicted
OPplot <- function(I, p,
                   lower.lim = 0,
                   upper.lim = NULL,
                   title = NULL,
                   picture.name = "photo.pdf",
                   save = TRUE) {

  df <- data.frame(
    cbind(I, rowMeans(p),
          t(apply(p, 1, quantile, prob = c(0.05, 0.95), type = 2))
    )
  )
  colnames(df) <- c("Xobs", "Xpre", "q05", "q95")

  df$col <- apply(df[,c(1,3,4)], 1, function(x) (x[1] - x[2]) * (x[1] - x[3]) < 0)
  print(mean(df$col))

  if (is.null(upper.lim)) upper.lim <- ceiling(max(df[,1:4]))

  gg <- ggplot2::ggplot(df) +
    ggplot2::geom_abline(intercept=0, slope=1) +
    ggplot2::geom_linerange(ggplot2::aes(y=Xobs, xmin=q05, xmax=q95, color=col), alpha=0.1) +
    ggplot2::geom_point(ggplot2::aes(x=Xpre, y=Xobs, color=col), alpha=0.8) +
    ggplot2::scale_color_manual(values=c("red", "black")) +
    ggplot2::labs(x='Predicted', y='Observed', title=title) +
    ggplot2::xlim(c(lower.lim, upper.lim)) +
    ggplot2::ylim(c(lower.lim, upper.lim)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none") +
    ggplot2::coord_fixed()

  if (save) {
    ggplot2::ggsave(picture.name, gg, width = 4, height = 4)
  }

  gg
}

i <- 1
OPplot(Nallobs[,i], Nallpre[,,i],
       lower.lim = 3.5,
       upper.lim = 8.5,
       title = bquote(bar(N)["64,month"](s)),
       picture.name = "inst/img/SUPP_OP_Nall_max.pdf")
OPplot(N21stobs[,i], N21stpre[,,i],
       title = bquote(bar(N)["42:64,month"](s)),
       picture.name = "inst/img/SUPP_OP_N21st_max.pdf")
OPplot(Robs[,i], Rpre[,,i],
       upper.lim = 8,
       title = bquote(bar(R)["55:64,month"](s)),
       picture.name = "inst/img/SUPP_OP_R_max.pdf")
OPplot(ERSobs[,i], ERSpre[,,i],
       upper.lim = 11,
       title = expression(t %*% bar(ERS)["t,month"](D)),
       picture.name = "inst/img/SUPP_OP_ERS_max.pdf")

i <- 2
OPplot(Nallobs[,i], Nallpre[,,i],
       lower.lim = 3.5,
       upper.lim = 8.5,
       title = bquote(underline(N)["64,month"](s)),
       picture.name = "inst/img/SUPP_OP_Nall_min.pdf")
OPplot(N21stobs[,i], N21stpre[,,i],
       title = bquote(underline(N)["42:64,month"](s)),
       picture.name = "inst/img/SUPP_OP_N21st_min.pdf")
OPplot(Robs[,i], Rpre[,,i],
       upper.lim = 8,
       title = bquote(underline(R)["55:64,month"](s)),
       picture.name = "inst/img/SUPP_OP_R_min.pdf")
OPplot(ERSobs[,i], ERSpre[,,i],
       upper.lim = 11,
       title = expression(t %*% underline(ERS)["t,month"](D)),
       picture.name = "inst/img/SUPP_OP_ERS_min.pdf")

i <- 3
OPplot(Nallobs[,i], Nallpre[,,i],
       lower.lim = 1,
       title = bquote(bar(underline(N))["64,month"](s)),
       picture.name = "inst/img/SUPP_OP_Nall_joint.pdf")
OPplot(N21stobs[,i], N21stpre[,,i],
       upper.lim = 1.5,
       title = bquote(bar(underline(N))["42:64,month"](s)),
       picture.name = "inst/img/SUPP_OP_N21st_joint.pdf")
OPplot(Robs[,i], Rpre[,,i],
       upper.lim = 0.75,
       title = bquote(bar(underline(N))["55:64,month"](s)),
       picture.name = "inst/img/SUPP_OP_R_joint.pdf")
OPplot(ERSobs[,i], ERSpre[,,i],
       upper.lim = 4.5,
       title = expression(t %*% bar(underline(ERS))["t,month"](D)),
       picture.name = "inst/img/SUPP_OP_ERS_joint.pdf")

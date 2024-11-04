##################################################
### Section 4 - RESULTS FULL MODEL (INFERENCE) ###
##################################################

library("sf")
library("sp")
library("stars")
library("smoothr")
library("tidyverse")
library("geomtextpath")
library("rnaturalearth")
library("rnaturalearthdata")

#' @param Z The values
#' @param coords_limit Boundaries of the map (typical xlim and ylim, better to
#'   leave them by default)
#' @param ref number: reference white color in the blue to red continuous color
#'   scale (not factor)
#' @param zlim limits of the Z if continuous (not factor)
#' @param picture.name name of the .png to be saved
#' @param title map title
#' @param legend.name legend title
#' @param save saves the map in the directory
#' @param contour one of c("none", "auto", "smooth")
#' @param breaks the contour levels will be drawn at these values * ref
#' @param label label for each break
#' @param smoothness level of smoothness if contour = "smooth"
#' @param threshold minimum length of the contour line to be drawn in meters
#' @param threshold2 if the length of the line in meters is smaller, the label
#'   is not drawn
#' @param dist vector of distance to the coast (if < 2 km the pixel is not
#'   considered in the contour line)
#' @return A map
mapSpain <- function(Z,
                     Zobs,
                     grid,
                     stations,
                     coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
                     ref = .5,
                     zlim = NULL,
                     picture.name = "photo.png",
                     title = expression(p["tl"](s)),
                     legend.name = "",
                     save = TRUE,
                     cross = FALSE,
                     contour = c("none", "auto", "smooth"),
                     breaks = 8:12 / 10,
                     label = "",
                     n.contour = 1,
                     smoothness = 5,
                     threshold = 1e+05,
                     threshold2 = 4 * threshold,
                     dist) {

  contour <- match.arg(contour)

  coords_limit <- st_transform(
    as(
      SpatialPointsDataFrame(coords = coords_limit,
                             data = coords_limit,
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
      'sf'),
    2062)

  # background
  world_map <- ne_countries(scale = "large", returnclass = 'sf')
  european_union <- c("Andorra","France","Morocco","Portugal","Spain","Gibraltar","Algeria")
  background <-
    world_map %>%
    filter(name %in% european_union) %>%
    st_transform(2062)

  map <- ggplot(data = background) +
    geom_sf(fill="antiquewhite") +
    xlab("Longitude") + ylab("Latitude") + ggtitle(title) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6,angle=90),
          axis.title=element_text(size=10,face="bold")) +
    geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[,1], y = st_coordinates(grid)[,2], fill = colMeans(Z, na.rm = TRUE)))

  if (!missing(Zobs)) {
    map <- map +
      geom_point(data = stations,
        aes(x = st_coordinates(stations)[, 1], y = st_coordinates(stations)[, 2], fill = Zobs),
        color = "black", pch=21, size=3)
  }

  map <- map +
    scale_fill_gradient2(midpoint = ref, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), space = "Lab", limits = zlim, name = legend.name)

  if (cross) {
    grid2 <- grid[apply(Z, 2, quantile, prob = 0.05) < 1 & 1 < apply(Z, 2, quantile, prob = 0.95), ]
    map <- map +
      geom_point(data = grid2, aes(x = st_coordinates(grid2)[, 1], y = st_coordinates(grid2)[, 2]), shape = 4, size = 0.8)
  }

  if (contour == "auto") {
    map <- map +
      geom_textcontour(data = grid,
                       mapping = ggplot2::aes(x = sf::st_coordinates(grid)[,1], y = sf::st_coordinates(grid)[,2], z = colMeans(Z)),
                       breaks = ref * breaks)
  } else if (contour == "smooth") {
    raster <- st_rasterize(st_sf(geometry = st_geometry(grid),
                                 value = colMeans(Z, na.rm = TRUE)))
    raster$value[raster$value == 0] <- NA
    contour1 <- stars::st_contour(stars::st_as_stars(raster), contour_lines = TRUE, breaks = ref * breaks)
    contour2 <- smoothr::smooth(contour1, method = "ksmooth", smoothness = smoothness)
    units(threshold) <- units(threshold2) <- "m"
    contour3 <- contour2[st_length(contour2) > threshold,]
    label <- label[match(round(contour3$value, 3), round(ref * breaks, 3))]
    label[st_length(contour3) < threshold2] <- ""

    map <- map +
      geom_textsf(data = contour3, label = label,
                  color = "black", linecolor = "black",
                  size = 2, linewidth = 0.275)
  }

  map <- map + coord_sf(xlim = st_coordinates(coords_limit)[,1], ylim = st_coordinates(coords_limit)[,2])# +
  #ggplot2::geom_point(ggplot2::aes(x = X, y = Y), data = data.frame(coords * 1000), color = "black")

  if (save) {
    ggplot2::ggsave(picture.name, map, width = 8.27 / 2, height = 11.69 / 4)
  }

  map
}

# ERS
ggERS <- function(y, xchar, ychar, ylim, xnumbreaks, xlabbreaks,
                  title = NULL,
                  picture.name = "photo.png",
                  save = TRUE) {

  n <- ncol(y) + 1
  df <- data.frame(
    y = c(1, 2:n * colMeans(y)),
    t = 1:n,
    CI1 = c(1, 2:n * apply(y, 2, quantile, prob = 0.05)),
    CI2 = c(1, 2:n * apply(y, 2, quantile, prob = 0.95))
  )

  gg <- ggplot(data = df,
               mapping =  aes(x = t, y = y)) +
    geom_hline(yintercept=1,
               color = "gray") +
    theme_bw() +
    theme(legend.position="none") +
    ylab(ychar) +
    xlab(xchar) +
    ylim(ylim) +
    scale_x_continuous(breaks=xnumbreaks,
                       labels=xlabbreaks) +
    geom_ribbon(ggplot2::aes(
      ymin = CI1, ymax = CI2),
      alpha = 0.2) +
    geom_line(size = 0.2)

  if (!is.null(title)) {
    gg <- gg + ggtitle(title)
  }

  if (save) {
    ggplot2::ggsave(picture.name, gg, width = 8.27 / 2, height = 11.69 / 4)
  }

  gg
}

# meter iteracion b-esima del MCMC
core <- function(TT, LL, SS, newSS,
                 sigma, mu,
                 prob_l0,
                 beta1, beta2,
                 w1tls, w2tls, decay,
                 a11s0, a22s0, a21s0,
                 dd, ddx,
                 grid) {

  ### re-scale regression coefficients
  beta1[-1] <- beta1[-1] / sigma
  beta2[-1] <- beta2[-1] / sigma
  beta1[1]  <- beta1[1] - sum(beta1[-1] * mu)
  beta2[1]  <- beta2[1] - sum(beta2[-1] * mu)

  ### matrices
  R11       <- exp(- decay[1] * dd[1:newSS, 1:newSS] -
                     decay[2] * ddx[1:newSS, 1:newSS])
  R12       <- exp(- decay[1] * dd[1:newSS, newSS + 1:SS] -
                     decay[2] * ddx[1:newSS, newSS + 1:SS])
  R22       <- exp(- decay[1] * dd[newSS + 1:SS, newSS + 1:SS] -
                     decay[2] * ddx[newSS + 1:SS, newSS + 1:SS])
  R22inv    <- solve(R22)
  R12R22inv <- R12 %*% R22inv
  R         <- R11 - R12R22inv %*% t(R12)

  MVN <- mvrnormArma(2 * TT * LL, Sigma = R)

  # l = 0
  prob_l0 <- log(prob_l0 / (1 - prob_l0)) + rnorm(2 * TT, sd = 0.2)
  prob_l0[prob_l0 == -Inf] <- rnorm(sum(prob_l0 == -Inf), mean = -8)
  prob_l0 <- 1 / (1 + exp(-prob_l0))

  #
  X  <- cbind(1, NA, NA, NA, log1p(grid$dist), NA, NA, NA)
  Xb <- matrix(nrow = newSS, ncol = 2)
  p  <- array(dim = c(TT, LL, newSS, 2, 2))

  for (t in 1:TT){

    X[,2]   <- qnorm(1 / (1 + t))
    X[,3:4] <- matrix(rbinom(2 * newSS, size = 1, prob = prob_l0[t,]), ncol = 2, byrow = TRUE)
    X[,6:8] <- X[1, 2] * X[,3:5]

    # l in JJA
    for (l in 1:LL) {

      w1tls0 <-
        R12R22inv %*% w1tls[t + (l - 1) * TT + 0:(SS - 1) * (LL * TT)] +
        MVN[2 * (l + (t - 1) * LL) - 1,]

      w2tls0 <-
        R12R22inv %*% w2tls[t + (l - 1) * TT + 0:(SS - 1) * (LL * TT)] +
        MVN[2 * (l + (t - 1) * LL),]

      Xb[,1] <- X %*% beta1 + a11s0 * w1tls0
      Xb[,2] <- X %*% beta2 + a21s0 * w1tls0 + a22s0 * w2tls0

      p[t, l, , , 1] <- pnorm(Xb)

      X[,3:4] <- p[t, l, , , 2] <- matrix(rbinom(2 * newSS, size = 1, prob = p[t, l,,,1]), ncol = 2)

      X[,6:7] <- X[1, 2] * X[,3:4]

    }
  }

  return(p)
}

library(Rcpp)
library(RcppArmadillo)
sourceCpp("inst/scripts/MVN.cpp")

pp <- array(dim = c(1000, 63, 92, 790, 2, 2)) # B T L s xn pI
set.seed(998877)
for (b in 1:500) {
  print(b)
  pp[b,,,,,] <- core(63, 92, 40, 790,
    sigma = attr(M5[,1]$x, "scaled:scale")[-1],
    mu    = attr(M5[,1]$x, "scaled:center")[-1],
    prob_l0 = cbind(apply(Ix[-1,151,], 1, mean),
                    apply(In[-1,151,], 1, mean)),
    beta1 = M5[,1]$params$beta1[b,],
    beta2 = M5[,1]$params$beta2[b,],
    w1tls = M5[,1]$params$w1tls[b,],
    w2tls = M5[,1]$params$w2tls[b,],
    decay = M5[,1]$params$decay[b,],
    a11s0 = a11s0[b,],
    a22s0 = a22s0[b,],
    a21s0 = a21s0[b,],
    dd  = spbrom:::dist1(rbind(sf::st_coordinates(grid), sf::st_coordinates(stations)) / 1000),
    ddx = spbrom:::dist1(as.matrix(c(grid$logSDx, log(stations$SDx)))),
    grid = grid)
  pp[b + 500,,,,,] <- core(63, 92, 40, 790,
    sigma = attr(M5[,2]$x, "scaled:scale")[-1],
    mu    = attr(M5[,2]$x, "scaled:center")[-1],
    prob_l0 = cbind(apply(Ix[-1,151,], 1, mean),
                    apply(In[-1,151,], 1, mean)),
    beta1 = M5[,2]$params$beta1[b,],
    beta2 = M5[,2]$params$beta2[b,],
    w1tls = M5[,2]$params$w1tls[b,],
    w2tls = M5[,2]$params$w2tls[b,],
    decay = M5[,2]$params$decay[b,],
    a11s0 = a11s0[b + 500,],
    a22s0 = a22s0[b + 500,],
    a21s0 = a21s0[b + 500,],
    dd  = spbrom:::dist1(rbind(sf::st_coordinates(grid), sf::st_coordinates(stations)) / 1000),
    ddx = spbrom:::dist1(as.matrix(c(grid$logSDx, log(stations$SDx)))),
    grid = grid)
}



# Jaccard index
t <- 55:64
Z <- apply(pp[,t-1,,,1,1] * pp[,t-1,,,2,1], c(1,4), mean) / apply(pp[,t-1,,,1,1] * pp[,t-1,,,2,1] + (1 - pp[,t-1,,,1,1]) * pp[,t-1,,,2,1] + pp[,t-1,,,1,1] * (1 - pp[,t-1,,,2,1]), c(1,4), mean)
ref <-  mean(1 / t * 1 / t) / mean(1 / t * 1 / t + (1 - 1 / t) * 1 / t + 1 / t * (1 - 1 / t))

mapSpain(Z,
         #Zobs = apply(Ix[t,JJA,] * In[t,JJA,], 3, mean) / apply((Ix[t,JJA,] + In[t,JJA,]) > 0, 3, mean),
         grid = grid,
         stations = stations,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 0.15,
         zlim = c(0, .25),
         title = expression(paste("Jaccard index (", T[x], ", ", T[n], ") (2014-2023)")),
         legend.name = "",
         save = TRUE,
         picture.name = "inst/img/MAIN_mapJ.png",
         contour = "smooth",
         breaks = c(.15, .2) / 0.15,
         label = c(.15, .2),
         smoothness = 3,
         threshold = 1e+05 / 2,
         threshold2 = 1e+05,
         dist = grid$dist)

Z <- matrix(nrow = 100, ncol = 63)
ref <- rep(NA, 63)
for (t in 2:64) {
  Z[,t-1] <- rowMeans(apply(pp[,t-1,,,1,1] * pp[,t-1,,,2,1], c(1,3), mean) / apply(pp[,t-1,,,1,1] * pp[,t-1,,,2,1] + (1 - pp[,t-1,,,1,1]) * pp[,t-1,,,2,1] + pp[,t-1,,,1,1] * (1 - pp[,t-1,,,2,1]), c(1,3), mean))
  ref[t-1] <-  mean(1 / t * 1 / t) / mean(1 / t * 1 / t + (1 - 1 / t) * 1 / t + 1 / t * (1 - 1 / t))
}

# Block average of Jaccard index
ggJ <- function(y, ref,
                xchar, ychar, ylim, xnumbreaks, xlabbreaks,
                title = NULL,
                picture.name = "photo.png",
                save = TRUE) {

  n <- ncol(y) + 1
  df <- data.frame(
    y = colMeans(y),
    ref = ref,
    t = 2:n,
    CI1 = apply(y, 2, quantile, prob = 0.05),
    CI2 = apply(y, 2, quantile, prob = 0.95)
  )

  gg <- ggplot(data = df,
               mapping = aes(x = t, y = y)) +
    geom_line(aes(x = t, y = ref),
              color = "gray") +
    theme_bw() +
    theme(legend.position="none") +
    ylab(ychar) +
    xlab(xchar) +
    ylim(ylim) +
    scale_x_continuous(breaks=xnumbreaks,
                       labels=xlabbreaks) +
    geom_ribbon(ggplot2::aes(
      ymin = CI1, ymax = CI2),
      alpha = 0.2) +
    geom_line(size = 0.2)

  if (!is.null(title)) {
    gg <- gg + ggtitle(title)
  }

  if (save) {
    ggplot2::ggsave(picture.name, gg, width = 8.27 / 2, height = 11.69 / 4)
  }

  gg
}

ggJ(Z, ref,
    "t (year)",
    expression(paste("Jaccard index (", T[x], ", ", T[n], ")")),
    c(0, 0.55),
    c(1, 21, 41, 61),
    c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
    picture.name = "inst/img/MAIN_tsJ.pdf")





# Persistence
Z <- matrix(NA, nrow = 100, ncol = 790)
t <- 44 # (2003)
for (b in 1:100) {
  for (i in 1:790) {
    if (sum(pp[b,t-1,-92,i,1,2] == 1) > 0) {
      Z[b,i] <-
        mean(pp[b,t-1,-1,i,2,1][pp[b,t-1,-92,i,1,2] == 1]) /
        mean(pp[b,t-1,-1,i,2,1][pp[b,t-1,-92,i,1,2] == 0])
    }
  }
}
range(sort(colMeans(Z)))

mapSpain(Z,
         grid = grid,
         stations = stations,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 5,
         zlim = c(2, 10.5),
         title = expression(frac(
           paste("P(", underline(I)["44,l"](s), "= 0 | ", bar(I)["44,l-1"](s), "= 1)"),
           paste("P(", underline(I)["44,l"](s), "= 0 | ", bar(I)["44,l-1"](s), "= 0)"))),
         legend.name = "",
         save = TRUE,
         picture.name = "inst/img/MAIN_map_persistence_2003.png",
         contour = "smooth",
         breaks = c(4,6,8) / 5,
         label = c(4,6,8),
         smoothness = 5,
         threshold = 2e+05 / 2,
         threshold2 = 2e+05,
         dist = grid$dist)

Z <- matrix(NA, nrow = 100, ncol = 790)
t <- 63 # (2022)
for (b in 1:100) {
  for (i in 1:790) {
    if (sum(pp[b,t-1,-92,i,1,2] == 1) > 0) {
      Z[b,i] <-
        mean(pp[b,t-1,-1,i,2,1][pp[b,t-1,-92,i,1,2] == 1]) /
        mean(pp[b,t-1,-1,i,2,1][pp[b,t-1,-92,i,1,2] == 0])
    }
  }
}
range(sort(colMeans(Z)))

mapSpain(Z,
         grid = grid,
         stations = stations,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 5,
         zlim = c(2, 10.5),
         title = expression(frac(
           paste("P(", underline(I)["63,l"](s), "= 0 | ", bar(I)["63,l-1"](s), "= 1)"),
           paste("P(", underline(I)["63,l"](s), "= 0 | ", bar(I)["63,l-1"](s), "= 0)"))),
         legend.name = "",
         save = TRUE,
         picture.name = "inst/img/MAIN_map_persistence_2022.png",
         contour = "smooth",
         breaks = c(4,6,8) / 5,
         label = c(4,6,8),
         smoothness = 5,
         threshold = 2e+05 / 2,
         threshold2 = 2e+05,
         dist = grid$dist)

Z <- array(NA, dim = c(100, 63, 790))
for (t in 2:64) {
  print(t)
  for (b in 1:100) {
    for (i in 1:790) {
      if (sum(pp[b,t-1,-92,i,1,2] == 1) > 0) {
        Z[b,t-1,i] <-
          mean(pp[b,t-1,-1,i,2,1][pp[b,t-1,-92,i,1,2] == 1]) /
          mean(pp[b,t-1,-1,i,2,1][pp[b,t-1,-92,i,1,2] == 0])
      }
    }
  }
}

ggJ(apply(Z, c(1, 2), mean, na.rm = TRUE), rep(1, 63),
    "t (year)",
    expression(frac(
      paste("P(", underline(I)["tl"](s), "= 0 | ", bar(I)["t,l-1"](s), "= 1)"),
      paste("P(", underline(I)["tl"](s), "= 0 | ", bar(I)["t,l-1"](s), "= 0)"))),
    c(2, 20.5),
    c(1, 21, 41, 61),
    c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
    picture.name = "inst/img/MAIN_ts_persistence.pdf")



### calendar heatmap
library("lubridate")
library("ggTimeSeries")

# Tx 1999-2023
calendar <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("1999-06-01"):as.Date("1999-08-31"),
      as.Date("2000-06-01"):as.Date("2000-08-31"),
      as.Date("2001-06-01"):as.Date("2001-08-31"),
      as.Date("2002-06-01"):as.Date("2002-08-31"),
      as.Date("2003-06-01"):as.Date("2003-08-31"),
      as.Date("2004-06-01"):as.Date("2004-08-31"),
      as.Date("2005-06-01"):as.Date("2005-08-31"),
      as.Date("2006-06-01"):as.Date("2006-08-31"),
      as.Date("2007-06-01"):as.Date("2007-08-31"),
      as.Date("2008-06-01"):as.Date("2008-08-31"),
      as.Date("2009-06-01"):as.Date("2009-08-31"),
      as.Date("2010-06-01"):as.Date("2010-08-31"),
      as.Date("2011-06-01"):as.Date("2011-08-31"),
      as.Date("2012-06-01"):as.Date("2012-08-31"),
      as.Date("2013-06-01"):as.Date("2013-08-31"),
      as.Date("2014-06-01"):as.Date("2014-08-31"),
      as.Date("2015-06-01"):as.Date("2015-08-31"),
      as.Date("2016-06-01"):as.Date("2016-08-31"),
      as.Date("2017-06-01"):as.Date("2017-08-31"),
      as.Date("2018-06-01"):as.Date("2018-08-31"),
      as.Date("2019-06-01"):as.Date("2019-08-31"),
      as.Date("2020-06-01"):as.Date("2020-08-31"),
      as.Date("2021-06-01"):as.Date("2021-08-31"),
      as.Date("2022-06-01"):as.Date("2022-08-31"),
      as.Date("2023-06-01"):as.Date("2023-08-31"))),
    value = c(
      apply(pp[,39,,,1,1], 2, mean),
      apply(pp[,40,,,1,1], 2, mean),
      apply(pp[,41,,,1,1], 2, mean),
      apply(pp[,42,,,1,1], 2, mean),
      apply(pp[,43,,,1,1], 2, mean),
      apply(pp[,44,,,1,1], 2, mean),
      apply(pp[,45,,,1,1], 2, mean),
      apply(pp[,46,,,1,1], 2, mean),
      apply(pp[,47,,,1,1], 2, mean),
      apply(pp[,48,,,1,1], 2, mean),
      apply(pp[,49,,,1,1], 2, mean),
      apply(pp[,50,,,1,1], 2, mean),
      apply(pp[,51,,,1,1], 2, mean),
      apply(pp[,52,,,1,1], 2, mean),
      apply(pp[,53,,,1,1], 2, mean),
      apply(pp[,54,,,1,1], 2, mean),
      apply(pp[,55,,,1,1], 2, mean),
      apply(pp[,56,,,1,1], 2, mean),
      apply(pp[,57,,,1,1], 2, mean),
      apply(pp[,58,,,1,1], 2, mean),
      apply(pp[,59,,,1,1], 2, mean),
      apply(pp[,60,,,1,1], 2, mean),
      apply(pp[,61,,,1,1], 2, mean),
      apply(pp[,62,,,1,1], 2, mean),
      apply(pp[,63,,,1,1], 2, mean))),
  "date",
  "value"
) +
  scale_fill_gradient2(
    midpoint = 0.1,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,1),
    name = "ERS"
  ) +
  theme_minimal()

ggplot2::ggsave("inst/img/SUPP_calendar_Tx.pdf", calendar,
                width = 11.69, height = 8.27)

# Tn 1999-2023
calendar <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("1999-06-01"):as.Date("1999-08-31"),
      as.Date("2000-06-01"):as.Date("2000-08-31"),
      as.Date("2001-06-01"):as.Date("2001-08-31"),
      as.Date("2002-06-01"):as.Date("2002-08-31"),
      as.Date("2003-06-01"):as.Date("2003-08-31"),
      as.Date("2004-06-01"):as.Date("2004-08-31"),
      as.Date("2005-06-01"):as.Date("2005-08-31"),
      as.Date("2006-06-01"):as.Date("2006-08-31"),
      as.Date("2007-06-01"):as.Date("2007-08-31"),
      as.Date("2008-06-01"):as.Date("2008-08-31"),
      as.Date("2009-06-01"):as.Date("2009-08-31"),
      as.Date("2010-06-01"):as.Date("2010-08-31"),
      as.Date("2011-06-01"):as.Date("2011-08-31"),
      as.Date("2012-06-01"):as.Date("2012-08-31"),
      as.Date("2013-06-01"):as.Date("2013-08-31"),
      as.Date("2014-06-01"):as.Date("2014-08-31"),
      as.Date("2015-06-01"):as.Date("2015-08-31"),
      as.Date("2016-06-01"):as.Date("2016-08-31"),
      as.Date("2017-06-01"):as.Date("2017-08-31"),
      as.Date("2018-06-01"):as.Date("2018-08-31"),
      as.Date("2019-06-01"):as.Date("2019-08-31"),
      as.Date("2020-06-01"):as.Date("2020-08-31"),
      as.Date("2021-06-01"):as.Date("2021-08-31"),
      as.Date("2022-06-01"):as.Date("2022-08-31"),
      as.Date("2023-06-01"):as.Date("2023-08-31"))),
    value = c(
      apply(pp[,39,,,2,1], 2, mean),
      apply(pp[,40,,,2,1], 2, mean),
      apply(pp[,41,,,2,1], 2, mean),
      apply(pp[,42,,,2,1], 2, mean),
      apply(pp[,43,,,2,1], 2, mean),
      apply(pp[,44,,,2,1], 2, mean),
      apply(pp[,45,,,2,1], 2, mean),
      apply(pp[,46,,,2,1], 2, mean),
      apply(pp[,47,,,2,1], 2, mean),
      apply(pp[,48,,,2,1], 2, mean),
      apply(pp[,49,,,2,1], 2, mean),
      apply(pp[,50,,,2,1], 2, mean),
      apply(pp[,51,,,2,1], 2, mean),
      apply(pp[,52,,,2,1], 2, mean),
      apply(pp[,53,,,2,1], 2, mean),
      apply(pp[,54,,,2,1], 2, mean),
      apply(pp[,55,,,2,1], 2, mean),
      apply(pp[,56,,,2,1], 2, mean),
      apply(pp[,57,,,2,1], 2, mean),
      apply(pp[,58,,,2,1], 2, mean),
      apply(pp[,59,,,2,1], 2, mean),
      apply(pp[,60,,,2,1], 2, mean),
      apply(pp[,61,,,2,1], 2, mean),
      apply(pp[,62,,,2,1], 2, mean),
      apply(pp[,63,,,2,1], 2, mean))),
  "date",
  "value"
) +
  scale_fill_gradient2(
    midpoint = 0.1,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,1),
    name = "ERS"
  ) +
  theme_minimal()

ggplot2::ggsave("inst/img/SUPP_calendar_Tn.pdf", calendar,
                width = 11.69, height = 8.27)

# Tx x Tn 1999-2023
calendar <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("1999-06-01"):as.Date("1999-08-31"),
      as.Date("2000-06-01"):as.Date("2000-08-31"),
      as.Date("2001-06-01"):as.Date("2001-08-31"),
      as.Date("2002-06-01"):as.Date("2002-08-31"),
      as.Date("2003-06-01"):as.Date("2003-08-31"),
      as.Date("2004-06-01"):as.Date("2004-08-31"),
      as.Date("2005-06-01"):as.Date("2005-08-31"),
      as.Date("2006-06-01"):as.Date("2006-08-31"),
      as.Date("2007-06-01"):as.Date("2007-08-31"),
      as.Date("2008-06-01"):as.Date("2008-08-31"),
      as.Date("2009-06-01"):as.Date("2009-08-31"),
      as.Date("2010-06-01"):as.Date("2010-08-31"),
      as.Date("2011-06-01"):as.Date("2011-08-31"),
      as.Date("2012-06-01"):as.Date("2012-08-31"),
      as.Date("2013-06-01"):as.Date("2013-08-31"),
      as.Date("2014-06-01"):as.Date("2014-08-31"),
      as.Date("2015-06-01"):as.Date("2015-08-31"),
      as.Date("2016-06-01"):as.Date("2016-08-31"),
      as.Date("2017-06-01"):as.Date("2017-08-31"),
      as.Date("2018-06-01"):as.Date("2018-08-31"),
      as.Date("2019-06-01"):as.Date("2019-08-31"),
      as.Date("2020-06-01"):as.Date("2020-08-31"),
      as.Date("2021-06-01"):as.Date("2021-08-31"),
      as.Date("2022-06-01"):as.Date("2022-08-31"),
      as.Date("2023-06-01"):as.Date("2023-08-31"))),
    value = c(
      apply(pp[,39,,,1,1] * pp[,39,,,2,1], 2, mean),
      apply(pp[,40,,,1,1] * pp[,40,,,2,1], 2, mean),
      apply(pp[,41,,,1,1] * pp[,41,,,2,1], 2, mean),
      apply(pp[,42,,,1,1] * pp[,42,,,2,1], 2, mean),
      apply(pp[,43,,,1,1] * pp[,43,,,2,1], 2, mean),
      apply(pp[,44,,,1,1] * pp[,44,,,2,1], 2, mean),
      apply(pp[,45,,,1,1] * pp[,45,,,2,1], 2, mean),
      apply(pp[,46,,,1,1] * pp[,46,,,2,1], 2, mean),
      apply(pp[,47,,,1,1] * pp[,47,,,2,1], 2, mean),
      apply(pp[,48,,,1,1] * pp[,48,,,2,1], 2, mean),
      apply(pp[,49,,,1,1] * pp[,49,,,2,1], 2, mean),
      apply(pp[,50,,,1,1] * pp[,50,,,2,1], 2, mean),
      apply(pp[,51,,,1,1] * pp[,51,,,2,1], 2, mean),
      apply(pp[,52,,,1,1] * pp[,52,,,2,1], 2, mean),
      apply(pp[,53,,,1,1] * pp[,53,,,2,1], 2, mean),
      apply(pp[,54,,,1,1] * pp[,54,,,2,1], 2, mean),
      apply(pp[,55,,,1,1] * pp[,55,,,2,1], 2, mean),
      apply(pp[,56,,,1,1] * pp[,56,,,2,1], 2, mean),
      apply(pp[,57,,,1,1] * pp[,57,,,2,1], 2, mean),
      apply(pp[,58,,,1,1] * pp[,58,,,2,1], 2, mean),
      apply(pp[,59,,,1,1] * pp[,59,,,2,1], 2, mean),
      apply(pp[,60,,,1,1] * pp[,60,,,2,1], 2, mean),
      apply(pp[,61,,,1,1] * pp[,61,,,2,1], 2, mean),
      apply(pp[,62,,,1,1] * pp[,62,,,2,1], 2, mean),
      apply(pp[,63,,,1,1] * pp[,63,,,2,1], 2, mean))),
  "date",
  "value"
) +
  scale_fill_gradient2(
    midpoint = 0.05,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,.5),
    name = "ERS"
  ) +
  theme_minimal()

ggplot2::ggsave("inst/img/SUPP_calendar_Tx_Tn_scale.pdf", calendar,
                width = 11.69, height = 8.27)

# Tx - Tn - Tx Tn - 2003, 2013, 2022, 2023
calendar <- list()
calendar[[1]] <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("2003-06-01"):as.Date("2003-08-31"),
      as.Date("2013-06-01"):as.Date("2013-08-31"))),
    value = c(
      apply(pp[,43,,,1,1], 2, mean),
      apply(pp[,53,,,1,1], 2, mean))),
  "date",
  "value"
) + xlab("") + ylab(expression(T[x])) +
  scale_fill_gradient2(
    midpoint = 0.1,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,1),
    name = "ERS"
  ) + theme_minimal() +
  theme(plot.margin = unit(c(0,-0.7,0,0), 'lines'),
    axis.title.y = element_text(angle = 0, vjust = 1.2), axis.text.x = element_blank(),
    axis.ticks = element_blank(), strip.background = element_blank(), legend.position="none")

calendar[[2]] <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("2022-06-01"):as.Date("2022-08-31"),
      as.Date("2023-06-01"):as.Date("2023-08-31"))),
    value = c(
      apply(pp[,62,,,1,1], 2, mean),
      apply(pp[,63,,,1,1], 2, mean))),
  "date",
  "value"
) + xlab("") + ylab("") +
  scale_fill_gradient2(
    midpoint = 0.1,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,1),
    name = "ERS"
  ) + theme_minimal() +
  theme(
    legend.key.size = unit(.5, 'cm'),
    legend.key.height = unit(.5, 'cm'),
    legend.key.width = unit(.5, 'cm'),
    legend.title = element_text(size=6),
    legend.text = element_text(size=6),
    plot.margin = unit(c(0,0,0,-0.7), 'lines'),
    axis.text = element_blank(), axis.ticks = element_blank(),
    strip.background = element_blank())

calendar[[3]] <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("2003-06-01"):as.Date("2003-08-31"),
      as.Date("2013-06-01"):as.Date("2013-08-31"))),
    value = c(
      apply(pp[,43,,,2,1], 2, mean),
      apply(pp[,53,,,2,1], 2, mean))),
  "date",
  "value"
) + xlab("") + ylab(expression(T[n])) +
  scale_fill_gradient2(
    midpoint = 0.1,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,1),
    name = "ERS"
  ) +
  theme_minimal() +
  theme(plot.margin = unit(c(0,-0.7,0,0), 'lines'),
    axis.title.y = element_text(angle = 0, vjust = 1.2), axis.text.x = element_blank(),
    axis.ticks = element_blank(), strip.background = element_blank(), legend.position="none")

calendar[[4]] <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("2022-06-01"):as.Date("2022-08-31"),
      as.Date("2023-06-01"):as.Date("2023-08-31"))),
    value = c(
      apply(pp[,62,,,2,1], 2, mean),
      apply(pp[,63,,,2,1], 2, mean))),
  "date",
  "value"
) + xlab("") + ylab("") +
  scale_fill_gradient2(
    midpoint = 0.1,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,1),
    name = "ERS"
  ) +
  theme_minimal() +
  theme(
    legend.key.size = unit(.5, 'cm'),
    legend.key.height = unit(.5, 'cm'),
    legend.key.width = unit(.5, 'cm'),
    legend.title = element_text(size=6),
    legend.text = element_text(size=6),
    plot.margin = unit(c(0,0,0,-0.7), 'lines'),
    axis.text = element_blank(), axis.ticks = element_blank(),
    strip.background = element_blank())

calendar[[5]] <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("2003-06-01"):as.Date("2003-08-31"),
      as.Date("2013-06-01"):as.Date("2013-08-31"))),
    value = c(
      apply(pp[,43,,,1,1] * pp[,43,,,2,1], 2, mean),
      apply(pp[,53,,,1,1] * pp[,53,,,2,1], 2, mean))),
  "date",
  "value"
) + xlab("") + ylab(expression(T[x] %*% T[n])) +
  scale_fill_gradient2(
    midpoint = 0.05,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,.5),
    name = "ERS"
  ) +
  theme_minimal() +
  theme(plot.margin = unit(c(0,-0.7,0,0), 'lines'),
    axis.title.y = element_text(angle = 0, vjust = 1.2), axis.ticks = element_blank(),
    strip.background = element_blank(), legend.position="none")

calendar[[6]] <- ggplot_calendar_heatmap(
  data.frame(
    date = as.Date(c(
      as.Date("2022-06-01"):as.Date("2022-08-31"),
      as.Date("2023-06-01"):as.Date("2023-08-31"))),
    value = c(
      apply(pp[,62,,,1,1] * pp[,62,,,2,1], 2, mean),
      apply(pp[,63,,,1,1] * pp[,63,,,2,1], 2, mean))),
  "date",
  "value"
) + xlab("") + ylab("") +
  scale_fill_gradient2(
    midpoint = 0.05,
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    limits = c(0,.5),
    name = "ERS"
  ) +
  theme_minimal() +
  theme(
    legend.key.size = unit(.5, 'cm'),
    legend.key.height = unit(.5, 'cm'),
    legend.key.width = unit(.5, 'cm'),
    legend.title = element_text(size=6),
    legend.text = element_text(size=6),
    plot.margin = unit(c(0,0,0,-0.7), 'lines'),
    axis.text.y = element_blank(), axis.ticks = element_blank(),
    strip.background = element_blank())

library(ggpubr)
library(cowplot)
calendar1 <- align_plots(calendar[[1]],
                         calendar[[3]],
                         calendar[[5]],
                         align = "hv",
                         axis="rlbt")
calendar2 <- align_plots(calendar[[2]],
                         calendar[[4]],
                         calendar[[6]],
                         align = "hv",
                         axis="rlbt")
calendar <- annotate_figure(
  ggarrange(
    ggarrange(plotlist = calendar1,
              nrow = 3, ncol=1, common.legend = FALSE),
    ggarrange(plotlist = calendar2,
              nrow = 3, ncol=1, common.legend = FALSE),
    nrow = 1, ncol = 2, align = "hv"
  ),
  left =  text_grob("Day of week", rot = 90),
  bottom = text_grob("Month", vjust = -1)
  )

ggplot2::ggsave("inst/img/MAIN_calendar_all.pdf", calendar,
                width = 11.69, height = 8.27 / 1.4)




# map daily
for (l in 1:10) {
  map <- mapSpain(apply(pp[,63,l+77,,1,1], 1:2, mean),
                Zobs = Ix[64,l+228,],
                grid = grid,
                stations = stations,
                coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
                ref = .5,
                zlim = c(0, 1),
                title = bquote(p[.(paste0("64,",l+228))](s)),
                legend.name = "",
                save = FALSE,
                contour = c("none", rep("auto", 8), "none")[l],
                breaks = c(.1, .2, .5, .8, .9) / .5,
                label = c(.1, .2, .5, .8, .9),
                smoothness = 5,
                threshold = 1e+05 / 2,
                threshold2 = 1e+05,
                dist = grid$dist)

  ggplot2::ggsave(paste0("inst/img/MAIN_Tx_pt64l", l+228, ".png"),
                  map, width = 8.27 / 2, height = 11.69 / 4)

  map <- mapSpain(apply(pp[,63,l+77,,2,1], 1:2, mean),
                  Zobs = In[64,l+228,],
                  grid = grid,
                  stations = stations,
                  coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
                  ref = .5,
                  zlim = c(0, 1),
                  title = bquote(p[.(paste0("64,",l+228))](s)),
                  legend.name = "",
                  save = FALSE,
                  contour = c("none", rep("auto", 8), "none")[l],
                  breaks = c(.1, .2, .5, .8, .9) / .5,
                  label = c(.1, .2, .5, .8, .9),
                  smoothness = 5,
                  threshold = 1e+05 / 2,
                  threshold2 = 1e+05,
                  dist = grid$dist)

  ggplot2::ggsave(paste0("inst/img/MAIN_Tn_pt64l", l+228, ".png"),
                  map, width = 8.27 / 2, height = 11.69 / 4)

  map <- mapSpain(apply(pp[,63,l+77,,1,1] * pp[,63,l+77,,2,1], 1:2, mean),
                  Zobs = Ix[64,l+228,] * In[64,l+228,],
                  grid = grid,
                  stations = stations,
                  coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
                  ref = .5,
                  zlim = c(0, 1),
                  title = bquote(p[.(paste0("64,",l+228))](s)),
                  legend.name = "",
                  save = FALSE,
                  contour = c("none", rep("auto", 8), "none")[l],
                  breaks = c(.1, .2, .5, .8, .9) / .5,
                  label = c(.1, .2, .5, .8, .9),
                  smoothness = 5,
                  threshold = 1e+05 / 2,
                  threshold2 = 1e+05,
                  dist = grid$dist)

  ggplot2::ggsave(paste0("inst/img/MAIN_Tx_Tn_pt64l", l+228, ".png"),
                  map, width = 8.27 / 2, height = 11.69 / 4)
}




## N max > N min
t <- 55:64
Z1 <- apply(pp[,t-1,,,1,2], c(1,4), mean)
Z2 <- apply(pp[,t-1,,,2,2], c(1,4), mean)
mapSpain(matrix(apply(Z1 > Z2, 2, mean) + apply(Z1 == Z2, 2, mean) / 2, nrow = 1),
         grid = grid,
         stations = stations,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = .5,
         zlim = c(0, 1),
         title = bquote(paste("Pr(", bar(N)["55:64,JJA"](s) > underline(N)["55:64,JJA"](s), ") + ",
                              frac(1, 2), "Pr(", bar(N)["55:64,JJA"](s) == underline(N)["55:64,JJA"](s), ")")),
         legend.name = "",
         save = TRUE,
         picture.name = "inst/img/MAIN_map_diff_Tx_Tn_2023_2014.png",
         contour = "auto",
         breaks = c(.1, .2, .5, .8, .9) / .5,
         label = c(.1, .2, .5, .8, .9),
         smoothness = 5,
         threshold = 1e+05 / 2,
         threshold2 = 1e+05,
         dist = grid$dist)

t <- 55:64 - 10
Z1 <- apply(pp[,t-1,,,1,2], c(1,4), mean)
Z2 <- apply(pp[,t-1,,,2,2], c(1,4), mean)
mapSpain(matrix(apply(Z1 > Z2, 2, mean) + apply(Z1 == Z2, 2, mean) / 2, nrow = 1),
         grid = grid,
         stations = stations,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = .5,
         zlim = c(0, 1),
         title = bquote(paste("Pr(", bar(N)["45:54,JJA"](s) > underline(N)["45:54,JJA"](s), ") + ",
                              frac(1, 2), "Pr(", bar(N)["45:54,JJA"](s) == underline(N)["45:54,JJA"](s), ")")),
         legend.name = "",
         save = TRUE,
         picture.name = "inst/img/MAIN_map_diff_Tx_Tn_2013_2004.png",
         contour = "auto",
         breaks = c(.1, .2, .5, .8, .9) / .5,
         label = c(.1, .2, .5, .8, .9),
         smoothness = 5,
         threshold = 1e+05 / 2,
         threshold2 = 1e+05,
         dist = grid$dist)

t <- 55:64 - 20
Z1 <- apply(pp[,t-1,,,1,2], c(1,4), mean)
Z2 <- apply(pp[,t-1,,,2,2], c(1,4), mean)
mapSpain(matrix(apply(Z1 > Z2, 2, mean) + apply(Z1 == Z2, 2, mean) / 2, nrow = 1),
         grid = grid,
         stations = stations,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = .5,
         zlim = c(0, 1),
         title = bquote(paste("Pr(", bar(N)["35:44,JJA"](s) > underline(N)["35:44,JJA"](s), ") + ",
                              frac(1, 2), "Pr(", bar(N)["35:44,JJA"](s) == underline(N)["35:44,JJA"](s), ")")),
         legend.name = "",
         save = TRUE,
         picture.name = "inst/img/MAIN_map_diff_Tx_Tn_2003_1994.png",
         contour = "auto",
         breaks = c(.1, .2, .5, .8, .9) / .5,
         label = c(.1, .2, .5, .8, .9),
         smoothness = 5,
         threshold = 1e+05 / 2,
         threshold2 = 1e+05,
         dist = grid$dist)


Z <- rep(0.5, 64)
for (t in 2:64) {
  Z1 <- apply(pp[,t-1,,,1,2], c(1,3), mean)
  Z2 <- apply(pp[,t-1,,,2,2], c(1,3), mean)
  Z[t] <- mean(apply(Z1 > Z2, 2, mean) + apply(Z1 == Z2, 2, mean) / 2)
}

ggJ(matrix(Z[-1], nrow = 1), rep(0.5, 63),
    "t (year)",
    bquote(paste("Pr(", bar(N)["t,JJA"](s) > underline(N)["t,JJA"](s), ") + ",
                 frac(1, 2), "Pr(", bar(N)["t,JJA"](s) == underline(N)["t,JJA"](s), ")")),
    c(0, 1),
    c(1, 21, 41, 61),
    c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
    picture.name = "inst/img/MAIN_ts_diff_Tx_Tn.pdf")




# More joint records now?
set.seed(223344)
Z <- array(rbinom(100 * 63 * 92 * 790, size = 1, prob = pp[,,,,1,1] * pp[,,,,2,1]), dim = c(100, 63, 92, 790))

t1 <- 55:64
Z1 <- apply(Z[,t1-1,,], c(1,4), mean)
for (tt in 1:4) {
  t2 <- t1 - tt * 10
  Z2 <- apply(Z[,t2-1,,], c(1,4), mean)
  mapSpain(matrix(apply(Z1 > Z2, 2, mean) + apply(Z1 == Z2, 2, mean) / 2, nrow = 1),
           grid = grid,
           stations = stations,
           coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
           ref = .5,
           zlim = c(0, 1),
           title = bquote(paste("Pr(", underline(bar(N))[.(paste0(t1[1],":",t1[10],",JJA"))](s) > underline(bar(N))[.(paste0(t2[1],":",t2[10],",JJA"))](s), ") + ",
                    frac(1, 2), "Pr(", underline(bar(N))[.(paste0(t1[1],":",t1[10],",JJA"))](s) == underline(bar(N))[.(paste0(t2[1],":",t2[10],",JJA"))](s), ")")),
           legend.name = "",
           save = TRUE,
           picture.name = paste0("inst/img/MAIN_map_diff_joint_last_", tt,".png"),
           contour = "auto",
           breaks = c(.1, .2, .5, .8, .9) / .5,
           label = c(.1, .2, .5, .8, .9),
           smoothness = 5,
           threshold = 1e+05 / 2,
           threshold2 = 1e+05,
           dist = grid$dist)
}


Zaux <- rep(0.5, 64)
for (t in 1:54) {
  t1 <- 55:64
  t2 <- 1:10 + t
  Z1 <- apply(Z[,t1-1,,], c(1,4), mean)
  Z2 <- apply(Z[,t2-1,,], c(1,4), mean)
  Zaux[t+10] <- mean(apply(Z1 > Z2, 2, mean) + apply(Z1 == Z2, 2, mean) / 2)
}

ggJ(matrix(Zaux[-1], nrow = 1), rep(0.5, 63),
    "t (year)",
    "Last decade > decade starting on year t",
    c(0, 1),
    c(1, 21, 41, 61),
    c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
    picture.name = "inst/img/MAIN_ts_diff_joint_last.pdf")


















ggplot(data = background) +
  geom_sf(fill = "antiquewhite") +
  xlab("Longitude") + ylab("Latitude") + ggtitle("Log-standard deviation (tx)") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.x=element_text(size = 6),
        axis.text.y=element_text(size = 6, angle = 90),
        axis.title=element_text(size = 10, face = "bold"),
        legend.justification = "left") +
  geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[, 1], y = st_coordinates(grid)[, 2], fill = apply(pp[,54:63,,,1,1], 4, mean)/ mean(1/55:64))) +
  geom_point(data = stations, aes(x = st_coordinates(stations)[, 1], y = st_coordinates(stations)[, 2],
                                  fill = apply(Ix[55:64,JJA,], 3, mean)/ mean(1/55:64)), color = "black", pch=21, size=5) +
  scale_fill_gradient2(midpoint = 1, low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                       space = "Lab", limits = c(0,5), name = expression(log(hat(s)[x]))) +
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])

ind <- pp[3,63,84,,1] * pp[3,63,84,,2]# / mean(1/55:64)
#78:87

ind <- apply(pp[,54:63-2,,,2,1], 4, mean)/ mean(1/53:62)
indq05 <- apply(apply(pp[,54:63-2,,,2,2], c(1,4), mean), 2, quantile, prob = 0.05)/ mean(1/53:62)
indq95 <- apply(apply(pp[,54:63-2,,,2,2], c(1,4), mean), 2, quantile, prob = 0.95)/ mean(1/53:62)
plot(indq05, ylim=c(min(indq05), max(indq95)), type="l")
lines(indq95, col="red")
abline(h=1)

ggplot(data = background) +
  geom_sf(fill = "antiquewhite") +
  xlab("Longitude") + ylab("Latitude") + ggtitle("Log-standard deviation (tx)") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.x=element_text(size = 6),
        axis.text.y=element_text(size = 6, angle = 90),
        axis.title=element_text(size = 10, face = "bold"),
        legend.justification = "left") +
  geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[, 1], y = st_coordinates(grid)[, 2], fill = ind)) + # apply(pp[,61,75,,1], 2, mean)
  geom_point(data = stations, aes(x = st_coordinates(stations)[, 1], y = st_coordinates(stations)[, 2],
                                  fill = apply(Ix[55:64,JJA,], 3, mean) / mean(1/55:64)), color = "black", pch=21, size=5) +
  scale_fill_gradient2(midpoint = 1, low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                       space = "Lab", limits = c(0,5), name = expression(log(hat(s)[x]))) +
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) +
  geom_point(data = grid[indq05 < 1 & indq95 > 1, ], aes(x = st_coordinates(grid)[indq05 < 1 & indq95 > 1, 1], y = st_coordinates(grid)[indq05 < 1 & indq95 > 1, 2]), shape = 4)

borrar <- mapSpain(Z = apply(pp[,54:63,,,1], c(1,4), mean) / mean(1/55:64),
         Zobs = apply(Ix[55:64,JJA,], 3, mean) / mean(1/55:64),
         grid = grid,
         stations = stations,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 6),
         picture.name = "inst/img/MAIN_Rxt5564JJAs.png",
         title = expression(paste(bar(R)["55:64,JJA"](s), " (mean)")),
         legend.name = "",
         save = TRUE)

i <- 64
plot(apply(Ix[i,JJA,], 1, mean), type = "l", ylim = c(0, 1))
lines(apply(In[i,JJA,], 1, mean), col = "red")
lines(apply(Ix[i,JJA,] * In[i,JJA,], 1, mean), col = "blue")

78:87

plot(1:64 * apply(Ix[,JJA,], 1, mean), type = "l", ylim = c(0, 8))
lines(1:64 * apply(In[,JJA,], 1, mean), col = "red")
lines(1:64 * apply(Ix[,JJA,] * In[,JJA,], 1, mean), col = "blue")


plot(
apply(Ix[64,JJA,], 1, mean),
apply(pp[,63,,,1,1], 2, mean)
)
abline(a=0,b=1)








# map number of records
mapSpain(rowMeans(pred[[2]][,1,]) / 365 + 1,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 1:62),
         zlim = c(4.55, 6),
         picture.name = "inst/img/MAIN_Nt62s.png",
         title = expression(paste(bar(N)["62"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 1:4 / 20,
         label = paste0(10 * 1:4 / 2, "%"),
         smoothness = 5,
         threshold = 1.5 * 1e+05,
         threshold2 = 1.5 * 1e+05,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,1,], 1, quantile, prob = 0.05) / 365 + 1,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 1:62),
         zlim = c(4.55, 6),
         picture.name = "inst/img/MAIN_Nt62s-q05.png",
         title = expression(paste(bar(N)["62"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:3 / 20,
         label =  c("CRM", paste0(10 * 1:3 / 2, "%")),
         smoothness = 5,
         threshold = 1.5 * 1e+05,
         threshold2 = 1.5 * 1e+05,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,1,], 1, quantile, prob = 0.95) / 365 + 1,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 1:62),
         zlim = c(4.55, 6),
         picture.name = "inst/img/MAIN_Nt62s-q95.png",
         title = expression(paste(bar(N)["62"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 1:5 / 20,
         label = paste0(10 * 1:5 / 2, "%"),
         smoothness = 5,
         threshold = 1.5 * 1e+05,
         threshold2 = 1.5 * 1e+05,
         dist = grid$dist)

summary(rowMeans(pred[[2]][,1,]) / 365 + 1, digits = 3)

interest <- colMeans((pred[[2]][,1,] / 365 + 1) > sum(1 / 1:62))
mean(interest)
quantile(interest, prob = c(0.05, 0.95))

mean(apply(pred[[2]][,1,] / 365 + 1, 1, quantile, prob = 0.05) > sum(1 / 1:62))
mean(apply(pred[[2]][,1,] / 365 + 1, 1, quantile, prob = 0.95) < sum(1 / 1:62))



# map number of records 21st century
mapSpain(rowMeans(pred[[2]][,2,]) / 90,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lDJFs.png",
         title = expression(paste(bar(N)["41:62,DJF"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:7 / 10,
         label = paste0(10 * 3:7, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,2,], 1, quantile, prob = 0.05) / 90,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lDJFs-q05.png",
         title = expression(paste(bar(N)["41:62,DJF"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:4 / 10,
         label = c("CRM", paste0(10 * 1:4, "%")),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,2,], 1, quantile, prob = 0.95) / 90,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lDJFs-q95.png",
         title = expression(paste(bar(N)["41:62,DJF"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 6:10 / 10,
         label = paste0(10 * 6:10, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,3,]) / 92,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lMAMs.png",
         title = expression(paste(bar(N)["41:62,MAM"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 2:9 / 10,
         label = paste0(10 * 2:9, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,3,], 1, quantile, prob = 0.05) / 92,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lMAMs-q05.png",
         title = expression(paste(bar(N)["41:62,MAM"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:6 / 10,
         label = c("CRM", paste0(10 * 1:6, "%")),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,3,], 1, quantile, prob = 0.95) / 92,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lMAMs-q95.png",
         title = expression(paste(bar(N)["41:62,MAM"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 5:12 / 10,
         label = paste0(10 * 5:12, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,4,]) / 92,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lJJAs.png",
         title = expression(paste(bar(N)["41:62,JJA"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 5:14 / 10,
         label = paste0(10 * 5:14, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,4,], 1, quantile, prob = 0.05) / 92,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lJJAs-q05.png",
         title = expression(paste(bar(N)["41:62,JJA"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 1:11 / 10,
         label = paste0(10 * 1:11, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,4,], 1, quantile, prob = 0.95) / 92,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lJJAs-q95.png",
         title = expression(paste(bar(N)["41:62,JJA"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 8:18 / 10,
         label = paste0(10 * 8:18, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,5,]) / 91,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lSONs.png",
         title = expression(paste(bar(N)["41:62,SON"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:8 / 10,
         label = paste0(10 * 3:8, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,5,], 1, quantile, prob = 0.05) / 91,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lSONs-q05.png",
         title = expression(paste(bar(N)["41:62,SON"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:5 / 10,
         label = c("CRM", paste0(10 * 1:5, "%")),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,5,], 1, quantile, prob = 0.95) / 91,
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = sum(1 / 41:62),
         zlim = c(0, 1.3),
         picture.name = "inst/img/MAIN_Nt41-62lSONs-q95.png",
         title = expression(paste(bar(N)["41:62,SON"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 6:11 / 10,
         label = paste0(10 * 6:11, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

## all year
interest <- colMeans((pred[[2]][,2,] + pred[[2]][,3,] + pred[[2]][,4,] + pred[[2]][,5,]) / 365)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## winter
interest <- colMeans(pred[[2]][,2,] / 90)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## spring
interest <- colMeans(pred[[2]][,3,] / 92)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## summer
interest <- colMeans(pred[[2]][,4,] / 92)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## autumn
interest <- colMeans(pred[[2]][,5,] / 91)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)

interest <- colMeans(((pred[[2]][,2,] + pred[[2]][,3,] + pred[[2]][,4,] + pred[[2]][,5,]) / 365) > sum(1 / 41:62))
mean(apply(((pred[[2]][,2,] + pred[[2]][,3,] + pred[[2]][,4,] + pred[[2]][,5,]) / 365), 1, quantile, prob = 0.025) > sum(1 / 41:62))



# maps of ratios in the last decade
mapSpain(rowMeans(pred[[2]][,6,]) / 90 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lDJF.png",
         title = expression(paste(R["53:62,DJF"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:5 / 5,
         label = 1 + 3:5 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,6,], 1, quantile, prob = 0.05) / 90 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lDJF-q05.png",
         title = expression(paste(R["53:62,DJF"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 1:2 / 5,
         label = 1 + 1:2 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,6,], 1, quantile, prob = 0.95) / 90 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lDJF-q95.png",
         title = expression(paste(R["53:62,DJF"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 6:9 / 5,
         label = 1 + 6:9 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,7,]) / 92 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lMAM.png",
         title = expression(paste(R["53:62,MAM"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 2:6 / 5,
         label = 1 + 2:6 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,7,], 1, quantile, prob = 0.05) / 92 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lMAM-q05.png",
         title = expression(paste(R["53:62,MAM"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + -1:4 / 5,
         label = 1 + -1:4 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,7,], 1, quantile, prob = 0.95) / 92 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lMAM-q95.png",
         title = expression(paste(R["53:62,MAM"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 4:9 / 5,
         label = 1 + 4:9 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,8,]) / 92 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lJJA.png",
         title = expression(paste(R["53:62,JJA"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:9 / 5,
         label = 1 + 3:9 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,8,], 1, quantile, prob = 0.05) / 92 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lJJA-q05.png",
         title = expression(paste(R["53:62,JJA"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:6 / 5,
         label = 1 + 0:6 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,8,], 1, quantile, prob = 0.95) / 92 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lJJA-q95.png",
         title = expression(paste(R["53:62,JJA"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 5:13 / 5,
         label = 1 + 5:13 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,9,]) / 91 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lSON.png",
         title = expression(paste(R["53:62,SON"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:6 / 5,
         label = 1 + 3:6 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,9,], 1, quantile, prob = 0.05) / 91 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lSON-q05.png",
         title = expression(paste(R["53:62,SON"](s), " (0.05 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:8 / 5,
         label = 1 + 0:8 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,9,], 1, quantile, prob = 0.95) / 91 / sum(1 / 53:62),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = 1,
         zlim = c(0, 3.8),
         picture.name = "inst/img/MAIN_Rt53-62lSON-q95.png",
         title = expression(paste(R["53:62,SON"](s), " (0.95 quantile)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 6:14 / 5,
         label = 1 + 6:14 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

## all year
interest <- colMeans((pred[[2]][,6,] + pred[[2]][,7,] + pred[[2]][,8,] + pred[[2]][,9,]) / 365 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## winter
interest <- colMeans(pred[[2]][,6,] / 90 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## spring
interest <- colMeans(pred[[2]][,7,] / 92 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## summer
interest <- colMeans(pred[[2]][,8,] / 92 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## autumn
interest <- colMeans(pred[[2]][,9,] / 91 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)

interest <- colMeans(((pred[[2]][,6,] + pred[[2]][,7,] + pred[[2]][,8,] + pred[[2]][,9,]) / 365) > sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
mean(apply(((pred[[2]][,6,] + pred[[2]][,7,] + pred[[2]][,8,] + pred[[2]][,9,]) / 365), 1, quantile, prob = 0.05) > sum(1 / 53:62))
mean(apply((pred[[2]][,6,] / 90), 1, quantile, prob = 0.05) > sum(1 / 53:62))
mean(apply((pred[[2]][,7,] / 92), 1, quantile, prob = 0.05) > sum(1 / 53:62))
mean(apply((pred[[2]][,8,] / 92), 1, quantile, prob = 0.05) > sum(1 / 53:62))
mean(apply((pred[[2]][,9,] / 91), 1, quantile, prob = 0.05) > sum(1 / 53:62))



# ERS
ggERS <- function(y, xchar, ychar, ylim, xnumbreaks, xlabbreaks,
                  title = NULL,
                  picture.name = "photo.png",
                  save = TRUE) {

  n <- ncol(y) + 1
  df <- data.frame(
    y = c(1, 2:n * colMeans(y)),
    t = 1:n,
    CI1 = c(1, 2:n * apply(y, 2, quantile, prob = 0.05)),
    CI2 = c(1, 2:n * apply(y, 2, quantile, prob = 0.95))
  )

  gg <- ggplot(data = df,
         mapping =  aes(x = t, y = y)) +
    geom_hline(yintercept=1,
               color = "gray") +
    theme_bw() +
    theme(legend.position="none") +
    ylab(ychar) +
    xlab(xchar) +
    ylim(ylim) +
    scale_x_continuous(breaks=xnumbreaks,
                       labels=xlabbreaks) +
    geom_ribbon(ggplot2::aes(
      ymin = CI1, ymax = CI2),
      alpha = 0.2) +
    geom_line(size = 0.2)

  if (!is.null(title)) {
    gg <- gg + ggtitle(title)
  }

  if (save) {
    ggplot2::ggsave(picture.name, gg, width = 8.27 / 2, height = 11.69 / 4)
  }

  gg
}

ggERS(apply(pp[,,,,1,2], 1:2, mean),
      "t (year)",
      expression(t %*% bar(ERS)[t](D)),
      c(0, 7.5),
      c(1, 21, 41, 61),
      c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
      picture.name = "inst/img/MAIN_ERS_Tx.pdf")

ggERS(apply(pp[,,,,2,2], 1:2, mean),
      "t (year)",
      expression(t %*% underline(ERS)[t](D)),
      c(0, 7.5),
      c(1, 21, 41, 61),
      c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
      picture.name = "inst/img/MAIN_ERS_Tn.pdf")

ggERS(apply(pp[,,,,1,2] * pp[,,,,2,2], 1:2, mean),
      "t (year)",
      expression(t %*% underline(bar(ERS))[t](D)),
      c(0, 7.5),
      c(1, 21, 41, 61),
      c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
      picture.name = "inst/img/MAIN_ERS_Tx_Tn.pdf")


# map daily
for (l in 1:8) {
  p <- mapSpain(apply(pp[,61,l+70,,1,1], 1:2, mean),
                Zobs = Ix[62,l+221,],
                grid = grid,
                stations = stations,
                coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
                ref = .5,
                zlim = c(0, 1),
                title = bquote(p[.(paste0("62,",l+221))](s)),
                legend.name = "",
                save = FALSE,
                contour = c("none", rep("auto", 7))[l],
                breaks = c(.1, .2, .5, .8, .9) / .5,
                label = c(.1, .2, .5, .8, .9),
                smoothness = 5,
                threshold = 1e+05 / 2,
                threshold2 = 1e+05,
                dist = grid$dist)

  ggplot2::ggsave(paste0("inst/img/MAIN_pt62l", l+221, ".png"),
                  p, width = 8.27 / 2, height = 11.69 / 4)
}


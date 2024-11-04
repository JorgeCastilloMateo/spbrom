#######################################################
### Section 4 - RESULTS SELECTED MODEL (PARAMETERS) ###
#######################################################

# This function returns the mean and quantiles by columns in X
my.summary <- function(X, prob = c(0.05, 0.95), round = 2) {

  if (is.matrix(X)) {
    x <- colMeans(X)
    y <- apply(X, 2, quantile, prob = prob)
    return(t(round(rbind("mean" = x, y), round)))
  } else {
    x <- mean(X)
    y <- quantile(X, prob = prob)
    return(round(c("mean" = x, y), round))
  }
}

# This function edits the output from my.summary to latex notation for a table
in.tex <- function(x) {
  internt.fun <- function(x) {
    if (as.numeric(x[2]) * as.numeric(x[3]) > 0) {
      paste0("\\texttt{", x[4], "} & $\\mathbf{", x[1], "}$ & $(", x[2], ",", x[3], ")$ \n")
    } else {
      paste0("\\texttt{", x[4], "} & $", x[1], "$ & $(", x[2], ",", x[3], ")$ \n")
    }

  }
  if (is.matrix(x)) {
    x <- data.frame(x, rownames(x))
    rownames(x) <- NULL
    return(apply(format(x, nsmall = 2), 1, internt.fun))
  } else {
    x <- format(x, nsmall = 2)
    x[4] <- "x"
    return(internt.fun(x))
  }
}

# Re-assign
model <- M5

for (i in 1:2) {
  # re-scale beta1
  model[,i]$params$beta1[,-1] <-
    sweep(model[,i]$params$beta1[,-1],  2, attr(model[,i]$x, "scaled:scale")[-1], FUN = '/')
  model[,i]$params$beta1[,1] <- model[,i]$params$beta1[,1] -
    rowSums(sweep(model[,i]$params$beta1[,-1],  2, attr(model[,i]$x, "scaled:center")[-1], FUN = '*'))
  # re-scale beta2
  model[,i]$params$beta2[,-1] <-
    sweep(model[,i]$params$beta2[,-1],  2, attr(model[,i]$x, "scaled:scale")[-1], FUN = '/')
  model[,i]$params$beta2[,1] <- model[,i]$params$beta2[,1] -
    rowSums(sweep(model[,i]$params$beta2[,-1],  2, attr(model[,i]$x, "scaled:center")[-1], FUN = '*'))
}



# MODEL PARAMETERS
## beta's
beta <- rbind(model[,1]$params$beta1, model[,2]$params$beta1)
cat(in.tex(my.summary(beta)))

beta <- rbind(model[,1]$params$beta2, model[,2]$params$beta2)
cat(in.tex(my.summary(beta)))

## decay's
decay <- rbind(model[,1]$params$decay, model[,2]$params$decay)
cat(in.tex(my.summary(3 / decay)))

## coregionalization a's (only for M4)
a11 <- c(model[,1]$params$a[,1], model[,2]$params$a[,1])
a22 <- c(model[,1]$params$a[,2], model[,2]$params$a[,2])
a21 <- c(model[,1]$params$a[,3], model[,2]$params$a[,3])
cat(in.tex(my.summary(a11)))
cat(in.tex(my.summary(a22)))
cat(in.tex(my.summary(a21)))
cat(in.tex(my.summary(a21 / sqrt(a21^2 + a22^2))))
cat(in.tex(my.summary(a11^2 / (a11^2 + 1))))
cat(in.tex(my.summary((a21^2 + a22^2) / (a21^2 + a22^2 + 1))))

## hyperparameter's (only for M5)
hpA <- rbind(model[,1]$params$hpA, model[,2]$params$hpA)
cat(in.tex(my.summary(hpA[,c(1,4,7)])))
cat(in.tex(my.summary(1 / sqrt(hpA[,c(2,5,8)]))))



# BOXPLOTS (for A(s), only for M5)
# Fix stations names
stations$NAME1 <-
  c("BADAJOZ", "MADRID (RETIRO)", "MALAGA", "NAVACERRADA", "SALAMANCA",
    "SAN SEBASTIAN", "TORTOSA", "VALENCIA", "ZARAGOZA", "BARCELONA (FABRA)",
    "ALBACETE", "BURGOS", "CIUDAD REAL", "CORUNA", "MURCIA", "SEVILLA",
    "SORIA", "BILBAO", "SANTIAGO", "PONFERRADA", "LEON", "LOGRONO", "ZAMORA",
    "REUS", "BARCELONA (AEROPUERTO)", "MADRID (TORREJON)", "VITORIA",
    "ALMERIA", "GIJON", "CACERES", "SANTANDER", "CASTELLON", "HUELVA",
    "LLEIDA", "MADRID (BARAJAS)", "MADRID (CUATROVIENTOS)", "MADRID (GETAFE)",
    "MORON", "VALLADOLID", "DAROCA")

ind <- order(stations$DIST)
a11 <- rbind(model[,1]$params$a[,ind], model[,2]$params$a[,ind])
a22 <- rbind(model[,1]$params$a[,ind + 40], model[,2]$params$a[,ind + 40])
a21 <- rbind(model[,1]$params$a[,ind + 80], model[,2]$params$a[,ind + 80])

## a11
png("inst/img/SUPP_boxplot_a11.png", width = 480 * 2, height = 480 * 2)
par(mar=c(10,3,4,1), cex = 2, cex.main = 1.75)
boxplot(a11, xaxt = "n", outline = FALSE)
title(expression(a["11"](s)))
axis(1, at = 1:40, labels = FALSE)
text(1:40, par("usr")[3] - diff(par("usr")[3:4]) * .05,
     srt = 60, adj = 1, xpd = TRUE,
     labels = stations$NAME1[ind], cex = 0.65)
dev.off()

## a22
png("inst/img/SUPP_boxplot_a22.png", width = 480 * 2, height = 480 * 2)
par(mar=c(10,3,4,1), cex = 2, cex.main = 1.75)
boxplot(a22, xaxt = "n", outline = FALSE)
title(expression(a["22"](s)))
axis(1, at = 1:40, labels = FALSE)
text(1:40, par("usr")[3] - diff(par("usr")[3:4]) * .05,
     srt = 60, adj = 1, xpd = TRUE,
     labels = stations$NAME1[ind], cex = 0.65)
dev.off()

## a21
png("inst/img/SUPP_boxplot_a21.png", width = 480 * 2, height = 480 * 2)
par(mar=c(10,3,4,1), cex = 2, cex.main = 1.75)
boxplot(a21, xaxt = "n", outline = FALSE)
title(expression(a["21"](s)))
axis(1, at = 1:40, labels = FALSE)
text(1:40, par("usr")[3] - diff(par("usr")[3:4]) * .05,
     srt = 60, adj = 1, xpd = TRUE,
     labels = stations$NAME1[ind], cex = 0.65)
dev.off()

## a (correlation)
png("inst/img/SUPP_boxplot_a.png", width = 480 * 2, height = 480 * 2)
par(mar=c(10,3,4,1), cex = 2, cex.main = 1.75)
boxplot(a21 / sqrt(a21^2 + a22^2), xaxt = "n", outline = FALSE)
title(expression(a(s)))
axis(1, at = 1:40, labels = FALSE)
text(1:40, par("usr")[3] - diff(par("usr")[3:4]) * .05,
     srt = 60, adj = 1, xpd = TRUE,
     labels = stations$NAME1[ind], cex = 0.65)
dev.off()

## a tx (prop spatial dep)
png("inst/img/SUPP_boxplot_a_tx.png", width = 480 * 2, height = 480 * 2)
par(mar=c(10,3,4,1), cex = 2, cex.main = 1.75)
boxplot(a11^2 / (a11^2 + 1), xaxt = "n", outline = FALSE,
        ylim = c(0.3, 1))
title(expression(bar(a)(s)))
axis(1, at = 1:40, labels = FALSE)
text(1:40, par("usr")[3] - diff(par("usr")[3:4]) * .05,
     srt = 60, adj = 1, xpd = TRUE,
     labels = stations$NAME1[ind], cex = 0.65)
dev.off()

## a tn (prop spatial dep)
png("inst/img/SUPP_boxplot_a_tn.png", width = 480 * 2, height = 480 * 2)
par(mar=c(10,3,4,1), cex = 2, cex.main = 1.75)
boxplot((a21^2 + a22^2) / (a21^2 + a22^2 + 1), xaxt = "n", outline = FALSE,
        ylim = c(0.3, 1))
title(expression(underline(a)(s)))
axis(1, at = 1:40, labels = FALSE)
text(1:40, par("usr")[3] - diff(par("usr")[3:4]) * .05,
     srt = 60, adj = 1, xpd = TRUE,
     labels = stations$NAME1[ind], cex = 0.65)
dev.off()



# MAPS and BLOCK AVERAGES (for A(s), only for M5)
a11 <- rbind(model[,1]$params$a[,1:40], model[,2]$params$a[,1:40])
a22 <- rbind(model[,1]$params$a[,1:40 + 40], model[,2]$params$a[,1:40 + 40])
a21 <- rbind(model[,1]$params$a[,1:40 + 80], model[,2]$params$a[,1:40 + 80])
hpA <- rbind(model[,1]$params$hpA, model[,2]$params$hpA)

## distance and correlation matrices (common decay)
SS <- 40
newSS <- nrow(grid)
d <- spbrom:::dist1(rbind(sf::st_coordinates(grid) / 1000, coords))
R11       <- exp(- hpA[1,"decaya11"] * d[1:newSS, 1:newSS])
R12       <- exp(- hpA[1,"decaya11"] * d[1:newSS, newSS + 1:SS])
R22       <- exp(- hpA[1,"decaya11"] * d[newSS + 1:SS, newSS + 1:SS])
R22inv    <- solve(R22)
R12R22inv <- R12 %*% R22inv
R         <- R11 - R12R22inv %*% t(R12)

library(Rcpp)
library(RcppArmadillo)
sourceCpp("inst/scripts/MVN.cpp")

set.seed(986484)
MVN <- mvrnormArma(3 * 1000, Sigma = R)

for (b in 1:1000) {
  media <- cbind(1, log1p(grid$dist)) %*% hpA[b,1:2] +
    R12R22inv %*% matrix(log(a11[b,]) - cbind(1, log1p(stations$DIST)) %*% hpA[b,1:2], ncol = 1)
  MVN[b,] <- exp(media + MVN[b,] / sqrt(hpA[b,"preca11"]))

  media <- cbind(1, log1p(grid$dist)) %*% hpA[b,5:6] +
    R12R22inv %*% matrix(log(a22[b,]) - cbind(1, log1p(stations$DIST)) %*% hpA[b,5:6], ncol = 1)
  MVN[b + 1000,] <- exp(media + MVN[b + 1000,] / sqrt(hpA[b,"preca22"]))

  media <- cbind(1, log1p(grid$dist)) %*% hpA[b,9:10] +
    R12R22inv %*% matrix(a21[b,] - cbind(1, log1p(stations$DIST)) %*% hpA[b,9:10], ncol = 1)
  MVN[b + 2000,] <- media + MVN[b + 2000,] / sqrt(hpA[b,"preca21"])
}

a11s0 <- MVN[1:1000,]
a22s0 <- MVN[1:1000 + 1000,]
a21s0 <- MVN[1:1000 + 2000,]

save(a11s0, a22s0, a21s0, file = "data/As0.RData")


mapSpain <- function(Z,
                     coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
                     ref = .5,
                     zlim = NULL,
                     picture.name = "photo.png",
                     title = expression(p["tl"](s)),
                     legend.name = "",
                     save = TRUE,
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

  #define grid of Spain
  spain <- ne_countries(scale = "large", country = "Spain", returnclass = "sf")
  spain <- st_transform(spain, 2062)

  spain_coords <- Polygons(
    list(Polygon(st_coordinates(spain)[st_coordinates(spain)[,"L2"] == 3,1:2])),
    ID = "spain")
  spain_coords <- SpatialPolygons(list(spain_coords))
  spain_coords <- as(spain_coords, "sf")
  st_crs(spain_coords) <- st_crs(spain)

  grid <- st_make_grid(spain, cellsize = 25 * 1000, what = "centers")
  grid <- st_intersection(grid, spain_coords)

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
    geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[,1], y = st_coordinates(grid)[,2], fill = Z)) +
    scale_fill_gradient2(midpoint = ref, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), space = "Lab", limits = zlim, name = legend.name)

  if (contour == "auto") {
    map <- map +
      geom_textcontour(data = grid,
                       mapping = ggplot2::aes(x = sf::st_coordinates(grid)[,1], y = sf::st_coordinates(grid)[,2], z = Z),
                       breaks = ref * breaks)
  } else if (contour == "smooth") {
    raster <- st_rasterize(st_sf(geometry = grid[dist > 2], value = Z[dist > 2]))
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

mapSpain(colMeans(a11s0),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = mean(a11s0),
         zlim = NULL,
         picture.name = "inst/img/SUPP_map_a11.png",
         title = expression(a["11"](s)),
         legend.name = "",
         save = TRUE)

mapSpain(colMeans(a22s0),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = mean(a22s0),
         zlim = NULL,
         picture.name = "inst/img/SUPP_map_a22.png",
         title = expression(a["22"](s)),
         legend.name = "",
         save = TRUE)

mapSpain(colMeans(a21s0),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = mean(a21s0),
         zlim = NULL,
         picture.name = "inst/img/SUPP_map_a21.png",
         title = expression(a["21"](s)),
         legend.name = "",
         save = TRUE)

mapSpain(colMeans(a21s0 / sqrt(a21s0^2 + a22s0^2)),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = mean(a21s0 / sqrt(a21s0^2 + a22s0^2)),
         zlim = NULL,
         picture.name = "inst/img/SUPP_map_a.png",
         title = expression(a(s)),
         legend.name = "",
         save = TRUE)

mapSpain(colMeans(a11s0^2 / (a11s0^2 + 1)),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = mean(a11s0^2 / (a11s0^2 + 1)),
         zlim = NULL,
         picture.name = "inst/img/SUPP_map_a_tx.png",
         title = expression(bar(a)(s)),
         legend.name = "",
         save = TRUE)

mapSpain(colMeans((a21s0^2 + a22s0^2) / (a21s0^2 + a22s0^2 + 1)),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = mean((a21s0^2 + a22s0^2) / (a21s0^2 + a22s0^2 + 1)),
         zlim = NULL,
         picture.name = "inst/img/SUPP_map_a_tn.png",
         title = expression(underline(a)(s)),
         legend.name = "",
         save = TRUE)

mapSpain(colMeans(a11s0^2 / (a11s0^2 + 1) - (a21s0^2 + a22s0^2) / (a21s0^2 + a22s0^2 + 1)),
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)),
         ref = mean(a11s0^2 / (a11s0^2 + 1) - (a21s0^2 + a22s0^2) / (a21s0^2 + a22s0^2 + 1)),
         zlim = c(0, 0.3),
         picture.name = "inst/img/SUPP_map_a_tx-tn.png",
         title = expression(bar(a)(s) - underline(a)(s)),
         legend.name = "",
         save = TRUE)

cat(in.tex(my.summary(rowMeans(a11s0))))
cat(in.tex(my.summary(rowMeans(a22s0))))
cat(in.tex(my.summary(rowMeans(a21s0))))
cat(in.tex(my.summary(rowMeans(a21s0 / sqrt(a21s0^2 + a22s0^2)))))
cat(in.tex(my.summary(rowMeans(a11s0^2 / (a11s0^2 + 1)))))
cat(in.tex(my.summary(rowMeans((a21s0^2 + a22s0^2) / (a21s0^2 + a22s0^2 + 1)))))



# SPATIAL CORRELATION (an-isotropic)
stations$NAME1 <-
  c("BADAJOZ", "MADRID (RETIRO)", "MALAGA", "NAVACERRADA", "SALAMANCA",
    "SAN SEBASTIAN", "TORTOSA", "VALENCIA", "ZARAGOZA", "BARCELONA (FABRA)",
    "ALBACETE", "BURGOS", "CIUDAD REAL", "CORUNA", "MURCIA", "SEVILLA",
    "SORIA", "BILBAO", "SANTIAGO", "PONFERRADA", "LEON", "LOGRONO", "ZAMORA",
    "REUS", "BARCELONA (AIRPORT)", "MADRID (TORREJON)", "VITORIA", "ALMERIA",
    "GIJON", "CACERES", "SANTANDER", "CASTELLON", "HUELVA", "LLEIDA",
    "MADRID (BARAJAS)", "MADRID (4VIENTOS)", "MADRID (GETAFE)", "MORON",
    "VALLADOLID", "DAROCA")

decayM3  <- c(M3[,1]$params$decay    , M3[,2]$params$decay    )
decayM5  <- c(M5[,1]$params$decay[,1], M5[,2]$params$decay[,1])
decayM5x <- c(M5[,1]$params$decay[,2], M5[,2]$params$decay[,2])

coords      <- sf::st_coordinates(stations) / 1000
coords.grid <- sf::st_coordinates(grid)     / 1000

for (i in 1:40) {
  print(i)
  coords_s_grid <- rbind(coords[i,], coords.grid)
  logS_s_grid   <- matrix(c(log(stations$SDx[i]), grid$logSDx), ncol = 1)

  r <- colMeans(exp(- decayM3 %*% spbrom:::dist1(coords_s_grid)[1, 2:791, drop = FALSE]))

  map <- ggplot(data = background) +
    geom_sf(fill = "antiquewhite") +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("Isotropic correlation with", subtitle = stations$NAME1[i]) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size = 6),
          axis.text.y=element_text(size = 6, angle = 90),
          axis.title=element_text(size = 10, face = "bold"),
          legend.position="none") +
    geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[, 1], y = st_coordinates(grid)[, 2], fill = r)) +
    scale_fill_gradient2(midpoint = 0.6, low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                         space = "Lab", limits = c(0, 1), name = "Correlation") +
    geom_point(aes(x = sf::st_coordinates(stations)[i,1], y = sf::st_coordinates(stations)[i,2]), shape = 4) +
    coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])

  ggplot2::ggsave(paste0("inst/img/SUPP_map_r", i, "_isotropic.png"), map, width = 8.27 / 2, height = 11.69 / 4)
}

for (i in 1:40) {
  print(i)
  coords_s_grid <- rbind(coords[i,], coords.grid)
  logS_s_grid   <- matrix(c(log(stations$SDx[i]), grid$logSDx), ncol = 1)

  r <- colMeans(exp(- decayM5 %*% spbrom:::dist1(coords_s_grid)[1, 2:791, drop = FALSE]) *
                exp(- decayM5x %*% spbrom:::dist1(logS_s_grid)[1, 2:791, drop = FALSE]))

  map <- ggplot(data = background) +
    geom_sf(fill = "antiquewhite") +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("Anisotropic correlation with", subtitle = stations$NAME1[i]) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size = 6),
          axis.text.y=element_text(size = 6, angle = 90),
          axis.title=element_text(size = 10, face = "bold"),
          legend.position="none") +
    geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[, 1], y = st_coordinates(grid)[, 2], fill = r)) +
    scale_fill_gradient2(midpoint = 0.6, low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                         space = "Lab", limits = c(0, 1), name = "Correlation") +
    geom_point(aes(x = sf::st_coordinates(stations)[i,1], y = sf::st_coordinates(stations)[i,2]), shape = 4) +
    coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2])

  ggplot2::ggsave(paste0("inst/img/SUPP_map_r", i, "_anisotropic.png"), map, width = 8.27 / 2, height = 11.69 / 4)
}


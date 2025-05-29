############################
### Section 2 - log(SDx) ###
############################

# Clear workspace
rm(list = setdiff(ls(),
  c("background", "data", "gg", "grid", "grid2",
    "In", "Ix", "limits", "stations")))

# OBTAIN log(SDx) (using kriging)
# Define stations as sp object
coords <- data.frame(st_coordinates(stations))
colnames(coords) <- c("x", "y")
stations2 <- data.frame(
  x = coords[,1], y = coords[,2],
  lon = stations$LON, lat = stations$LAT,
  elev = stations$HGHT, dist = stations$DIST,
  SDx = stations$SDx,
  row.names = stations$name)
coordinates(stations2) <- c("x", "y")
stations2@proj4string <- CRS("EPSG:2062")

# Krige the lodSDx surface
v.emp <- variogram(log(SDx) ~ log1p(elev) + log1p(dist), data = stations2,
                   cutoff = 800 * 1000, width = 50 * 1000)
v.fit <- fit.variogram(v.emp, vgm(model = "Exp"), fit.method = 6)

plot(v.emp, v.fit)

kriging <- krige(log(SDx) ~ log1p(elev) + log1p(dist), stations2, grid2, v.fit)

# Insert logSDx
grid$logSDx <- kriging$var1.pred

# Plot of logSDx
pp <- ggplot(data = background) +
  geom_sf(fill = "antiquewhite") +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle(expression("Log-standard deviation" ~ (T[x]))) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text.x=element_text(size = 6),
        axis.text.y=element_text(size = 6, angle = 90),
        axis.title=element_text(size = 10, face = "bold"),
        legend.justification = "left") +
  geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[, 1], y = st_coordinates(grid)[, 2], fill = logSDx)) +
  scale_fill_gradient2(midpoint = 1.5, low = scales::muted("blue"), mid = "white", high = scales::muted("red"),
                       space = "Lab", limits = c(0.8, 1.7), name = expression(widehat(log(s)))) +
  coord_sf(xlim = st_coordinates(limits)[, 1], ylim = st_coordinates(limits)[, 2]) +
  geom_point(data = stations,
             aes(x = st_coordinates(stations)[, 1], y = st_coordinates(stations)[, 2], fill = log(SDx)),
             color = "black", pch=21, size=2)

ggplot2::ggsave("inst/img/SUPP_map_logSDx.png", gg(pp), width = 8.27 / 2, height = 11.69 / 4)

# Cross-validation log(SDx)
K <- 10
SS <- 40
set.seed(63940)
folds <- matrix(nrow = SS / K, ncol = K)
remain <- 1:SS
for (k in 1:K) {
  folds[,k] <- sort(sample(remain, SS / K))
  remain    <- remain[!(remain %in% folds[,k])]
}

stations$logSDx.cv <- NA
for (k in 1:K) {
  v.emp <- variogram(log(SDx) ~ log1p(elev) + log1p(dist),
                     data = stations2[-folds[,k],],
                     cutoff = 800 * 1000, width = 50 * 1000)
  v.fit <- fit.variogram(v.emp, vgm(model = "Exp"), fit.method = 6)
  stations$logSDx.cv[folds[,k]] <-
    krige(log(SDx) ~ log1p(elev) + log1p(dist),
          stations2[-folds[,k],], stations2[folds[,k],], v.fit)$var1.pred
}
cor(stations$logSDx.cv, log(stations$SDx))^2
plot(x = stations$logSDx.cv,
     y = log(stations$SDx) - stations$logSDx.cv)



coords <- sf::st_coordinates(stations$geometry) / 1000


# Save data
save(coords, data, grid, In, Ix, stations, file = "data/main.RData")


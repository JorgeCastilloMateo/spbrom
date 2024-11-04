############################################
### Section 4 - 10-FOLD CROSS-VALIDATION ###
############################################

# Clear workspace
rm(list = setdiff(ls(),
  c("coords", "data", "grid", "In", "Ix", "stations")))

# Indexes
TT <- 64
LL <- 92
SS <- 40
tt <- 2:TT
K <- 10

# Choose folds at random
set.seed(63940)
folds <- matrix(nrow = SS / K, ncol = K)
remain <- 1:SS
for (k in 1:K) {
  folds[,k] <- sort(sample(remain, SS / K))
  remain    <- remain[!(remain %in% folds[,k])]
}

KFCVmetric <- array(dim = c(3 * 6, 12, 2),
  dimnames = list(
    model  = c(paste0("M", 0:5, "(max)"),
      paste0("M", 0:5, "(min)"),
      paste0("M", 0:5, "(joint)")),
    metric = c("AUC", "AUC1", "AUC2", "J", "J1", "J2",
               "AUPRC", "AUPRC1", "AUPRC2", "BS", "BS1", "BS2"),
    metric_type = c("CV", "pCV"))
)

METRICS <- array(dim = c(K, 12, 3, 2))



### M0
# Metrics
pp <- matrix(rep(1 / tt, LL * SS), 2, (TT - 1) * LL * SS, byrow = TRUE)
for (ii in 1:2) {
  KFCVmetric[ 1,,ii] <- spbrom::metrics(data$tx,           pp,      year = data$year, metric = dimnames(KFCVmetric)$metric_type[ii], D1 = 2:31, D2 = 32:64)
  KFCVmetric[ 7,,ii] <- spbrom::metrics(data$tn,           pp,      year = data$year, metric = dimnames(KFCVmetric)$metric_type[ii], D1 = 2:31, D2 = 32:64)
  KFCVmetric[13,,ii] <- spbrom::metrics(data$tx * data$tn, pp * pp, year = data$year, metric = dimnames(KFCVmetric)$metric_type[ii], D1 = 2:31, D2 = 32:64)
}

# Save metric
saveRDS(KFCVmetric, file = "data/KFCVmetric.rds")



### M1 (max)
KFCV_M1_tx <- list()
data$lag <- data$lag.tx
# Chain 1 (10 folds in parallel): < 18 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds"))
KFCV_M1_tx[[1]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter); spbrom::rom(
    formula = tx ~ trend * (lag + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    sp.A = FALSE,
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): < 18 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds"))
KFCV_M1_tx[[2]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter + 10); spbrom::rom(
    formula = tx ~ trend * (lag + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    sp.A = FALSE,
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence
MCMCconvergence(KFCV_M1_tx)

# Save data
saveRDS(KFCV_M1_tx, file = "data/KFCV_M1_tx.rds")

### M1 (min)
KFCV_M1_tn <- list()
data$lag <- data$lag.tn
# Chain 1 (10 folds in parallel): < 18 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds"))
KFCV_M1_tn[[1]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter + 100); spbrom::rom(
    formula = tn ~ trend * (lag + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    sp.A = FALSE,
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): < 18 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds"))
KFCV_M1_tn[[2]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter + 1000); spbrom::rom(
    formula = tn ~ trend * (lag + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    sp.A = FALSE,
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence
MCMCconvergence(KFCV_M1_tn)

# Save data
saveRDS(KFCV_M1_tn, file = "data/KFCV_M1_tn.rds")

# Metrics
p1 <- p2 <- array(dim = c(500, 23184, 2))
set.seed(12345)
for (k in 1:K) {
  print(k)
  data$lag <- data$lag.tx
  p1[,,1] <- spbrom::predict.rom(
    object  = KFCV_M1_tx[[1]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  data$lag <- data$lag.tn
  p1[,,2] <- spbrom::predict.rom(
    object  = KFCV_M1_tn[[1]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  data$lag <- data$lag.tx
  p2[,,1] <- spbrom::predict.rom(
    object  = KFCV_M1_tx[[2]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  data$lag <- data$lag.tn
  p2[,,2] <- spbrom::predict.rom(
    object  = KFCV_M1_tn[[2]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  for (ii in 1:2) {
    METRICS[k,,1,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,2,ii] <- spbrom::metrics(
      data$tn[data$site %in% folds[,k]],
      rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,3,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]] * data$tn[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]) * rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
  }
}
(KFCVmetric[ 2,,] <- apply(METRICS[,,1,], 2:3, mean))
(KFCVmetric[ 8,,] <- apply(METRICS[,,2,], 2:3, mean))
(KFCVmetric[14,,] <- apply(METRICS[,,3,], 2:3, mean))

# Save metric
saveRDS(KFCVmetric, file = "data/KFCVmetric.rds")



### M2
KFCV_M2 <- list()
# Chain 1 (10 folds in parallel): ~ 1.5 days
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds"))
KFCV_M2[[1]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    sp.A = FALSE,
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): ~ 1.5 days
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds"))
KFCV_M2[[2]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter + 10); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    sp.A = FALSE,
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence
MCMCconvergence(KFCV_M2)

# Save data
saveRDS(KFCV_M2, file = "data/KFCV_M2.rds")

# Metrics
set.seed(12345)
for (k in 1:K) {
  print(k)
  p1 <- spbrom::predict.brom(
    object  = KFCV_M2[[1]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  p2 <- spbrom::predict.brom(
    object  = KFCV_M2[[2]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  for (ii in 1:2) {
    METRICS[k,,1,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,2,ii] <- spbrom::metrics(
      data$tn[data$site %in% folds[,k]],
      rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,3,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]] * data$tn[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]) * rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
  }
}
(KFCVmetric[ 3,,] <- apply(METRICS[,,1,], 2:3, mean))
(KFCVmetric[ 9,,] <- apply(METRICS[,,2,], 2:3, mean))
(KFCVmetric[15,,] <- apply(METRICS[,,3,], 2:3, mean))

# Save metric
saveRDS(KFCVmetric, file = "data/KFCVmetric.rds")



### M3
KFCV_M3 <- list()
# Chain 1 (10 folds in parallel): ~ 1 day
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds", "stations"))
KFCV_M3[[1]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    sp.A = TRUE,
    matrix.mean.A = cbind(1, log1p(stations$DIST[!(1:SS %in% folds[, iter])])),
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, decay = c(2, 100),
                 sigma = c(0.1, 0.1), decay_A = 3 / 900),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): ~ 1 day
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds", "stations"))
KFCV_M3[[2]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter + 10); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    sp.A = TRUE,
    matrix.mean.A = cbind(1, log1p(stations$DIST[!(1:SS %in% folds[, iter])])),
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, decay = c(2, 100),
                 sigma = c(0.1, 0.1), decay_A = 3 / 900),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence
MCMCconvergence(KFCV_M3)

# Save data
saveRDS(KFCV_M3, file = "data/KFCV_M3.rds")

# Metrics
set.seed(12345)
for (k in 1:K) {
  print(k)
  p1 <- spbrom::predict.brom(
    object  = KFCV_M3[[1]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newmatrix.mean.A = cbind(1, log1p(stations$DIST[folds[, k]])),
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  p2 <- spbrom::predict.brom(
    object  = KFCV_M3[[2]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newmatrix.mean.A = cbind(1, log1p(stations$DIST[folds[, k]])),
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  for (ii in 1:2) {
    METRICS[k,,1,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,2,ii] <- spbrom::metrics(
      data$tn[data$site %in% folds[,k]],
      rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,3,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]] * data$tn[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]) * rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
  }
}
(KFCVmetric[ 4,,] <- apply(METRICS[,,1,], 2:3, mean))
(KFCVmetric[10,,] <- apply(METRICS[,,2,], 2:3, mean))
(KFCVmetric[16,,] <- apply(METRICS[,,3,], 2:3, mean))

# Save metric
saveRDS(KFCVmetric, file = "data/KFCVmetric.rds")



### M4
KFCV_M4 <- list()
# Chain 1 (10 folds in parallel): ~ 1.7 days
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds", "stations"))
KFCV_M4[[1]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    extra.coords = as.matrix(log(stations$SDx[!(1:SS %in% folds[, iter])])),
    sp.A = FALSE,
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 2, 100, 0.1)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): ~ 1.7 days
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds", "stations"))
KFCV_M4[[2]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter + 10); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    extra.coords = as.matrix(log(stations$SDx[!(1:SS %in% folds[, iter])])),
    sp.A = FALSE,
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 2, 100, 0.1)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence
MCMCconvergence(KFCV_M4)

# Save data
saveRDS(KFCV_M4, file = "data/KFCV_M4.rds")

# Metrics
set.seed(12345)
for (k in 1:K) {
  print(k)
  p1 <- spbrom::predict.brom(
    object  = KFCV_M4[[1]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newextra.coords = as.matrix(stations$logSDx.cv[folds[,k]]),
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  p2 <- spbrom::predict.brom(
    object  = KFCV_M4[[2]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newextra.coords = as.matrix(stations$logSDx.cv[folds[,k]]),
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  for (ii in 1:2) {
    METRICS[k,,1,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,2,ii] <- spbrom::metrics(
      data$tn[data$site %in% folds[,k]],
      rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,3,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]] * data$tn[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]) * rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
  }
}
(KFCVmetric[ 5,,] <- apply(METRICS[,,1,], 2:3, mean))
(KFCVmetric[11,,] <- apply(METRICS[,,2,], 2:3, mean))
(KFCVmetric[17,,] <- apply(METRICS[,,3,], 2:3, mean))

# Save metric
saveRDS(KFCVmetric, file = "data/KFCVmetric.rds")



### M5
KFCV_M5 <- list()
# Chain 1 (10 folds in parallel): ~ 1.2 days
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds", "stations"))
KFCV_M5[[1]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    extra.coords = as.matrix(log(stations$SDx[!(1:SS %in% folds[, iter])])),
    sp.A = TRUE,
    matrix.mean.A = cbind(1, log1p(stations$DIST[!(1:SS %in% folds[, iter])])),
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, decay = c(2, 2, 100, 0.1),
                 sigma = c(0.1, 0.1), decay_A = 3 / 900),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): ~ 1.2 days
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("data", "coords", "SS", "folds", "stations"))
KFCV_M5[[2]] <- parallel::parSapply(cl = cl, X = 1:K,
  FUN = function(iter) {set.seed(iter + 10); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    extra.coords = as.matrix(log(stations$SDx[!(1:SS %in% folds[, iter])])),
    sp.A = TRUE,
    matrix.mean.A = cbind(1, log1p(stations$DIST[!(1:SS %in% folds[, iter])])),
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    decay.prior = "gamma",
    prior = list(beta = 0.01, decay = c(2, 2, 100, 0.1),
                 sigma = c(0.1, 0.1), decay_A = 3 / 900),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence
MCMCconvergence(KFCV_M5)

# Save data
saveRDS(KFCV_M5, file = "data/KFCV_M5.rds")

# Metrics
set.seed(12345)
for (k in 1:K) {
  print(k)
  p1 <- spbrom::predict.brom(
    object  = KFCV_M5[[1]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newextra.coords = as.matrix(stations$logSDx.cv[folds[,k]]),
    newmatrix.mean.A = cbind(1, log1p(stations$DIST[folds[, k]])),
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  p2 <- spbrom::predict.brom(
    object  = KFCV_M5[[2]][,k],
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newextra.coords = as.matrix(stations$logSDx.cv[folds[,k]]),
    newmatrix.mean.A = cbind(1, log1p(stations$DIST[folds[, k]])),
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day,
    type = "KFCV")
  for (ii in 1:2) {
    METRICS[k,,1,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,2,ii] <- spbrom::metrics(
      data$tn[data$site %in% folds[,k]],
      rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
    METRICS[k,,3,ii] <- spbrom::metrics(
      data$tx[data$site %in% folds[,k]] * data$tn[data$site %in% folds[,k]],
      rbind(p1[,,1], p2[,,1]) * rbind(p1[,,2], p2[,,2]),
      data$year[data$site %in% folds[,k]],
      metric = dimnames(KFCVmetric)$metric_type[ii],
      D1 = 2:31, D2 = 32:64)
  }
}
(KFCVmetric[ 6,,] <- apply(METRICS[,,1,], 2:3, mean))
(KFCVmetric[12,,] <- apply(METRICS[,,2,], 2:3, mean))
(KFCVmetric[18,,] <- apply(METRICS[,,3,], 2:3, mean))

# Save metric
saveRDS(KFCVmetric, file = "data/KFCVmetric.rds")


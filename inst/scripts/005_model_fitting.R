##################################################
### Section 4 - MODEL FITTING (M1, M2, M3, M4) ###
##################################################

# M1 (max) - 2 chains in parallel < 20 hours
data$lag <- data$lag.tx
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
M1_tx <- parallel::parSapply(cl = cl, X = 1:2,
  FUN = function(iter) {set.seed(23 * iter); spbrom::rom(
    formula = tx ~ trend * (lag + log1p(dist)),
    data = data,
    coords = coords,
    sp.A = FALSE,
    site = data$site,
    year = data$year,
    day  = data$day,
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
MCMCconvergence(M1_tx)

# Save data
saveRDS(M1_tx, file = "data/M1_tx.rds")



# M1 (min) - 2 chains in parallel < 20 hours
data$lag <- data$lag.tn
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
M1_tn <- parallel::parSapply(cl = cl, X = 1:2,
  FUN = function(iter) {set.seed(23 * iter + 1000); spbrom::rom(
    formula = tn ~ trend * (lag + log1p(dist)),
    data = data,
    coords = coords,
    sp.A = FALSE,
    site = data$site,
    year = data$year,
    day  = data$day,
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
MCMCconvergence(M1_tn)

# Save data
saveRDS(M1_tn, file = "data/M1_tn.rds")



# M2 - 2 chains in parallel ~ 30 hours
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
M2 <- parallel::parSapply(cl = cl, X = 1:2,
  FUN = function(iter) {set.seed(23 * iter + 1000); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data,
    coords = coords,
    sp.A = FALSE,
    site = data$site,
    year = data$year,
    day  = data$day,
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 100)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
MCMCconvergence(M2)

# Save data
saveRDS(M2, file = "data/M2.rds")



# M3 - 2 chains in parallel ~ 30 hours
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
M3 <- parallel::parSapply(cl = cl, X = 1:2,
  FUN = function(iter) {set.seed(23 * iter + 1000); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data,
    coords = coords,
    sp.A = TRUE,
    matrix.mean.A = cbind(1, log1p(stations$DIST)),
    site = data$site,
    year = data$year,
    day  = data$day,
    decay.prior = "gamma",
    prior = list(beta = 0.01, decay = c(2, 100),
                 sigma = c(0.1, 0.1), decay_A = 3 / 900),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
MCMCconvergence(M3)

# Save data
saveRDS(M3, file = "data/M3.rds")



# M4 - 2 chains in parallel ~ 30 hours
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords", "stations"))
M4 <- parallel::parSapply(cl = cl, X = 1:2,
  FUN = function(iter) {set.seed(23 * iter + 1000); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data,
    coords = coords,
    extra.coords = as.matrix(log(stations$SDx)),
    sp.A = FALSE,
    site = data$site,
    year = data$year,
    day  = data$day,
    decay.prior = "gamma",
    prior = list(beta = 0.01, a = 1 / 5^2, decay = c(2, 2, 100, 0.1)),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
MCMCconvergence(M4)

# Save data
saveRDS(M4, file = "data/M4.rds")



# M5 - 2 chains in parallel ~ 30 hours
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords", "stations"))
M5 <- parallel::parSapply(cl = cl, X = 1:2,
  FUN = function(iter) {set.seed(23 * iter + 1000); spbrom::brom(
    formula = cbind(tx, tn) ~ trend * (lag.tx + lag.tn + log1p(dist)),
    data = data,
    coords = coords,
    extra.coords = as.matrix(log(stations$SDx)),
    sp.A = TRUE,
    matrix.mean.A = cbind(1, log1p(stations$DIST)),
    site = data$site,
    year = data$year,
    day  = data$day,
    decay.prior = "gamma",
    prior = list(beta = 0.01, decay = c(2, 2, 100, 0.1),
                 sigma = c(0.1, 0.1), decay_A = 3 / 900),
    n.report = 500001, n.burnin = 200000, n.sims = 400000, n.thin = 800)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
MCMCconvergence(M5)

# Save data
saveRDS(M5, file = "data/M5.rds")

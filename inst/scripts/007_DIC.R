#######################
###       DIC       ###
#######################

### DIC
TT <- 64
LL <- 92
SS <- 40
tt <- 2:TT

DICmetric <- matrix(nrow = 3 * 6, ncol = 3,
  dimnames = list(
    model  = c(paste0("M", 0:5, "(max)"),
               paste0("M", 0:5, "(min)"),
               paste0("M", 0:5, "(joint)")),
    metric = c("DIC", "D", "pD"))
)



# M0
pp <- matrix(rep(1 / tt, LL * SS), 2, (TT - 1) * LL * SS, byrow = TRUE)
DICmetric[ 1,] <- spbrom::metrics(data$tx, pp, metric = "DIC")
DICmetric[ 7,] <- spbrom::metrics(data$tn, pp, metric = "DIC")
DICmetric[13,] <- spbrom::metrics(data$tx * data$tn, pp * pp, metric = "DIC")

# save
saveRDS(DICmetric, file = "data/DICmetric.rds")



# M1 (max)
pp <- array(dim = c(1000, 231840, 2))
data$lag <- data$lag.tx
set.seed(12345)
pp[1:500,,1] <-
  spbrom::predict.rom(
    object  = M1_tx[, 1],
    newdata = data,
    type = "sample")
pp[501:1000,,1] <-
  spbrom::predict.rom(
    object  = M1_tx[, 2],
    newdata = data,
    type = "sample")
DICmetric[2,] <- spbrom::metrics(data$tx, pp[,,1], metric = "DIC")

# M1 (min)
data$lag <- data$lag.tn
set.seed(54321)
pp[1:500,,2] <-
  spbrom::predict.rom(
    object  = M1_tn[, 1],
    newdata = data,
    type = "sample")
pp[501:1000,,2] <-
  spbrom::predict.rom(
    object  = M1_tn[, 2],
    newdata = data,
    type = "sample")
DICmetric[8,] <- spbrom::metrics(data$tn, pp[,,2], metric = "DIC")

# M1 (joint)
DICmetric[14,] <- spbrom::metrics(data$tx * data$tn, pp[,,1] * pp[,,2], metric = "DIC")

# save
saveRDS(DICmetric, file = "data/DICmetric.rds")



# M2
pp <- array(dim = c(1000, 231840, 2))
set.seed(12345)
pp[1:500,,] <-
  spbrom::predict.brom(
    object  = M2[, 1],
    newdata = data,
    type = "sample")
pp[501:1000,,] <-
  spbrom::predict.brom(
    object  = M2[, 2],
    newdata = data,
    type = "sample")
DICmetric[ 3,] <- spbrom::metrics(data$tx, pp[,,1], metric = "DIC")
DICmetric[ 9,] <- spbrom::metrics(data$tn, pp[,,2], metric = "DIC")
DICmetric[15,] <- spbrom::metrics(data$tx * data$tn, pp[,,1] * pp[,,2], metric = "DIC")

# save
saveRDS(DICmetric, file = "data/DICmetric.rds")



# M3
pp <- array(dim = c(1000, 231840, 2))
set.seed(12345)
pp[1:500,,] <-
  spbrom::predict.brom(
    object  = M3[, 1],
    newdata = data,
    type = "sample")
pp[501:1000,,] <-
  spbrom::predict.brom(
    object  = M3[, 2],
    newdata = data,
    type = "sample")
DICmetric[ 4,] <- spbrom::metrics(data$tx, pp[,,1], metric = "DIC")
DICmetric[10,] <- spbrom::metrics(data$tn, pp[,,2], metric = "DIC")
DICmetric[16,] <- spbrom::metrics(data$tx * data$tn, pp[,,1] * pp[,,2], metric = "DIC")

# save
saveRDS(DICmetric, file = "data/DICmetric.rds")



# M4
pp <- array(dim = c(1000, 231840, 2))
set.seed(12345)
pp[1:500,,] <-
  spbrom::predict.brom(
    object  = M4[, 1],
    newdata = data,
    type = "sample")
pp[501:1000,,] <-
  spbrom::predict.brom(
    object  = M4[, 2],
    newdata = data,
    type = "sample")
DICmetric[ 5,] <- spbrom::metrics(data$tx, pp[,,1], metric = "DIC")
DICmetric[11,] <- spbrom::metrics(data$tn, pp[,,2], metric = "DIC")
DICmetric[17,] <- spbrom::metrics(data$tx * data$tn, pp[,,1] * pp[,,2], metric = "DIC")

# save
saveRDS(DICmetric, file = "data/DICmetric.rds")



# M5
pp <- array(dim = c(1000, 231840, 2))
set.seed(12345)
pp[1:500,,] <-
  spbrom::predict.brom(
    object  = M5[, 1],
    newdata = data,
    type = "sample")
pp[501:1000,,] <-
  spbrom::predict.brom(
    object  = M5[, 2],
    newdata = data,
    type = "sample")
DICmetric[ 6,] <- spbrom::metrics(data$tx, pp[,,1], metric = "DIC")
DICmetric[12,] <- spbrom::metrics(data$tn, pp[,,2], metric = "DIC")
DICmetric[18,] <- spbrom::metrics(data$tx * data$tn, pp[,,1] * pp[,,2], metric = "DIC")

# save
saveRDS(DICmetric, file = "data/DICmetric.rds")

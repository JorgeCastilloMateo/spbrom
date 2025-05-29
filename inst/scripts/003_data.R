###################################
### Section 2 - DATA PROCESSING ###
###################################

# Clear workspace
rm(list = setdiff(ls(),
  c("background", "gg", "grid", "grid2", "limits", "stations")))

# Load the R-package
library("spbrom")

# Read temperature data
library(readr)
tx <- read_csv("data/Tx_mat.csv")
tn <- read_csv("data/Tn_mat.csv")

# Indexes
TT <- 64
LL <- 365
SS <- 40
tt <- 2:TT
JJA <- 152:243

# temp data
tx365 <- array(dim = c(TT, LL, SS))
tn365 <- array(dim = c(TT, LL, SS))

for (i in 1:SS) {
  for (t in 1:TT) {
    tx365[t,,i] <- unlist(tx[1:LL + (t - 1) * LL, i + 1])
    tn365[t,,i] <- unlist(tn[1:LL + (t - 1) * LL, i + 1])
  }
}

tx365 <- tx365 / 10
tn365 <- tn365 / 10

# Compute sd JJA 1991-2020
stations$SDx <- apply(tx365[32:61, JJA, ], 3, sd, na.rm = TRUE)
stations$SDn <- apply(tn365[32:61, JJA, ], 3, sd, na.rm = TRUE)

# Compute indicators
Ix <- array(dim = c(TT, LL, SS))
In <- array(dim = c(TT, LL, SS))
for (i in 1:SS) {
  Ix[,,i] <- apply(tx365[,,i], 2, I.weak.record)
  In[,,i] <- apply(tn365[,,i], 2, I.weak.record)
}

# Main effects / indexes
data <- data.frame(
  tx      = c(Ix[-1, JJA, ]),
  tn      = c(In[-1, JJA, ]),
  year    = c(array(tt, dim = c(TT - 1, 92, SS))),
  day     = c(aperm(array(JJA, dim = c(92, TT - 1, SS)), perm = c(2, 1, 3))),
  yearday = NA,
  site    = c(aperm(array(1:SS, dim = c(SS, TT - 1, 92)), perm = c(2, 3, 1))),
  trend   = c(array(qnorm(1 / tt), dim = c(TT - 1, 92, SS))),
  sine    = c(aperm(array(sin(2 * pi * JJA / 365), dim = c(92, TT - 1, SS)), perm = c(2, 1, 3))),
  cosi    = c(aperm(array(cos(2 * pi * JJA / 365), dim = c(92, TT - 1, SS)), perm = c(2, 1, 3))),
  lon     = c(aperm(array(stations$LON, dim = c(SS, TT - 1, 92)), perm = c(2, 3, 1))),
  lat     = c(aperm(array(stations$LAT, dim = c(SS, TT - 1, 92)), perm = c(2, 3, 1))),
  dist    = c(aperm(array(stations$DIST, dim = c(SS, TT - 1, 92)), perm = c(2, 3, 1))),
  elev    = c(aperm(array(stations$HGHT, dim = c(SS, TT - 1, 92)), perm = c(2, 3, 1))),
  lag.tx  = c(Ix[-1, JJA - 1, ]),
  lag.tn  = c(In[-1, JJA - 1, ]))

data$yearday <- interaction(data$year, data$day)


# Save data
save(data, grid, In, Ix, stations, file = "data/main.RData")


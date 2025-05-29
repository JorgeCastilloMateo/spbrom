#######################
### Section 2 - EDA ###
#######################

### Script provided by:
###   Zeus Gracia-Tabuenca

rm(list = setdiff(ls(),
  c("background", "gg", "grid", "grid2", "limits", "stations")))

# Save graphical parameters
op <- par()

# Load libraries
if(!is.element("ggplot2", row.names(installed.packages()))) install.packages("ggplot2")
library("ggplot2")
if(!is.element("RColorBrewer", row.names(installed.packages()))) install.packages("RColorBrewer")
library("RColorBrewer")

# Set data directory
data_dir <- "data"
if(!dir.exists(data_dir)) stop("data_dir not found")

# Observatories of interest
stations$Zona <-
  c("MS","MS","A","MS","MN","C","L","L","E","L",
    "L","MN","MS","G","L","A","MN","C","G","MN",
    "MN","E","MN","L","L",NA,"E","A","C","MS",
    "C","L","A","E",NA,NA,NA,"A","MN","E")
idx <- which(!is.na(stations$Zona))
stations36 <- stations[idx,]
stations36$Zona <- factor(stations36$Zona)
levels(stations36$Zona) <- c("A","C","E","G","L","NP","SP")

# Data directory
outdir <- "inst/img/"
if(!dir.exists(outdir)) dir.create(outdir)

###############################################################################
# Figure 2
###############################################################################

# Trend - probability of record trends

# Read Tx data
tx <- read.csv(file = file.path(data_dir,"Tx_mat.csv"))
tn <- read.csv(file = file.path(data_dir,"Tn_mat.csv"))
# Reshape data
LL <- 365
TT <- nrow(tx)/LL
SS <- nrow(stations)
tx3d <- array(data = as.matrix(tx[,-1]), dim = c(LL,TT,SS))
if(sum(is.na(tx3d))>0) tx3d[is.na(tx3d)] <- -9999 # Remove NAs
tn3d <- array(data = as.matrix(tn[,-1]), dim = c(LL,TT,SS))
if(sum(is.na(tn3d))>0) tn3d[is.na(tn3d)] <- -9999 # Remove NAs

# Set seasons
summer_idx <- c(152:243)

# Compute upper records (indicators) for summer data only
upp.rcrd <- function(x) c(1,as.numeric(diff(cummax(x))>0))
itx3d <- apply(X = tx3d[summer_idx,,idx], MARGIN = c(1,3), upp.rcrd)
itn3d <- apply(X = tn3d[summer_idx,,idx], MARGIN = c(1,3), upp.rcrd)
voldim <- dim(itx3d)
if(!identical(voldim,dim(itn3d))) stop("TX and TN volumes dimension do not match!!!")

# Compute yearly trend
trend_tx <- apply(X = itx3d, MARGIN = 1, function(x) mean(x, na.rm = T))
trend_tn <- apply(X = itn3d, MARGIN = 1, function(x) mean(x, na.rm = T))
# Plot inputs
trnd_df <- data.frame(Itx=trend_tx,
                      Itn=trend_tn,
                      t=1:voldim[1])
# Plot
g1a <- ggplot(data = trnd_df,
             mapping =  aes(x=t)) +
  geom_hline(yintercept=1, color = "gray") +
  geom_line(aes(y=Itx*t, color="Tx")) +
  geom_line(aes(y=Itn*t, color="Tn")) +
  geom_smooth(aes(y=Itx*t, color="Tx"), se = F, linetype = "dashed") +
  geom_smooth(aes(y=Itn*t, color="Tn"), se = F, linetype = "dashed") +
  ylim(c(0, NA)) +
  ylab(expression(t%*%hat(p)[t])) +
  xlab("t (year)") +
  scale_x_continuous(breaks = c(1, 21, 41, 61),
                     labels = c("1 (1960)", "21 (1980)", "41 (2000)", " 61 (2020)")) +
  scale_color_manual(name = "",
                     breaks=c("Tx","Tn"),
                     values=c("Tx"="red", "Tn"="blue")) +
  theme_bw()
show(g1a)
ggsave(filename = "MAIN_trend_pt_txtn.pdf",
       plot = g1a,
       device = "pdf",
       path = outdir,
       width = 5,
       height = 3)

# Concurrence
# Load 3d volume functions

#' @title Lag one for margin two in 3D volume
#'
#' @description This function applies a lag-one into the second dimension of a
#' 3D volume object.
#'
#' @usage lag3d(vol3d)
#'
#' @param vol3d Input 3D volume.
#'
#' @details Default for lagged value is zero.
#'
#' @return Another 3D volume with same dimensions as input data.
#'
#' @author Zeus Gracia-Tabuenca
#'
#' @examples
#' \donttest{
#' Y365 <- readRDS("Data/Y365.rds")
#' I365 <- apply(X = Y365,
#'               MARGIN = c(2,3),
#'               FUN =  function(x) c(1,as.numeric(diff(cummax(x))>0)))
#' I365.lag1 <- lag3d(vol3D = I365)
#' }
#'
#' @export

lag3d <- function(vol3d){
  # Transpose
  vol3d_dim <- dim(vol3d)
  aux <- apply(vol3d, 3, function(x) t(x))
  # Apply lag1
  aux2 <- rbind(rep(0,ncol(aux)), aux[-nrow(aux),])
  # Transpose again
  aux3 <- array(data = apply(X = array(data = aux2, dim = vol3d_dim[c(2,1,3)]),
                             MARGIN = 3,
                             FUN = function(x) t(x)),
                dim = vol3d_dim)
  return(aux3)
}

#' @title Odds ratio calculation between two volumes
#'
#' @description This function calculates the odds ratio between two volumes
#' containing binary data.
#'
#' @usage or3d(vol1, vol2)
#'
#' @param vol1 Input 3D volume.
#' @param vol2 Input 3D volume.
#'
#' @details The OR is computed along the 2nd and 3rd dimensions.
#'
#' @return a data.frame summarizing the OR along the 1st dimension of both
#' volumes.
#'
#' @author Zeus Gracia-Tabuenca
#'
#' @examples
#' \donttest{
#' aux_df <- or3d(vol1, vol2)
#' }
#'
#' @export

or3d <- function(vol1, vol2){
  # Check dimension
  if(!identical(dim(vol1),dim(vol2))) stop("Dimensions differ!!")
  # Compute odds ratios by year
  n11 <- apply(vol1+vol2, 1, function(x) sum(x==2, na.rm = T))+0.5
  n00 <- apply(vol1+vol2, 1, function(x) sum(x==0, na.rm = T))+0.5
  n01 <- apply(vol1-vol2, 1, function(x) sum(x==-1, na.rm = T))+0.5
  n10 <- apply(vol1-vol2, 1, function(x) sum(x==1, na.rm = T))+0.5
  # Compute OR
  return(data.frame(OR=(n11*n00)/(n10*n01), t=1:dim(vol1)[1]))
}

# Compute lagged records
itx3d <- apply(X = tx3d[,,idx], MARGIN = c(1,3), upp.rcrd)
itn3d <- apply(X = tn3d[,,idx], MARGIN = c(1,3), upp.rcrd)
itx3d_lag1 <- lag3d(itx3d)       # Generate lag1 for TX
itn3d_lag1 <- lag3d(itn3d)       # Generate lag1 for Tn
# Extract only JJA
itx3d <- itx3d[,summer_idx,]
itn3d <- itn3d[,summer_idx,]
itx3d_lag1 <- itx3d_lag1[,summer_idx,]
itn3d_lag1 <- itn3d_lag1[,summer_idx,]

# Tx/Tn concurrence
aux_df <- or3d(itx3d,itn3d)
or_df <- data.frame(t = aux_df$t, or_xn.. = aux_df$OR)
# Tx/Tx(-1)
or_df$or_x.x. <- or3d(itx3d,itx3d_lag1)$OR
# Tn/Tn(-1)
or_df$or_.n.n <- or3d(itn3d,itn3d_lag1)$OR
# Tn/Tx(-1)
or_df$or_.nx. <- or3d(itn3d,itx3d_lag1)$OR
# Tx/Tn(-1)
or_df$or_x..n <- or3d(itx3d,itn3d_lag1)$OR

# Plot
colorcito <- brewer.pal(n = 8, name = "Set1")
g1b <- ggplot(data = or_df[-1,],
             mapping =  aes(x=t)) +
  geom_smooth(aes(y=log(or_xn..), color="xn|.."), se = F) +
  geom_smooth(aes(y=log(or_x.x.), color="x.|x'."), se = F) +
  geom_smooth(aes(y=log(or_.n.n), color=".n|.n'"), se = F) +
  geom_smooth(aes(y=log(or_.nx.), color=".n|x'."), se = F) +
  geom_smooth(aes(y=log(or_x..n), color="x.|.n'"), se = F) +
  #labs(title = "Tx,Tn concurrence and persistence") +
  ylab(expression(LOR[t])) +
  xlab("t (year)") +
  scale_x_continuous(breaks = c(1, 21, 41, 61),
                     labels = c("1 (1960)", "21 (1980)", "41 (2000)", " 61 (2020)")) +
  scale_color_manual(name = "",
                     breaks=c("xn|..",
                              "x.|x'.",
                              ".n|.n'",
                              ".n|x'.",
                              "x.|.n'"),
                     values=c("xn|.."=colorcito[1],
                              "x.|x'."=colorcito[2],
                              ".n|.n'"=colorcito[3],
                              ".n|x'."=colorcito[4],
                              "x.|.n'"=colorcito[5])) +
  theme_bw()
show(g1b)
# Save
ggsave(filename = "MAIN_lor_txtn_JJA.pdf",
       plot = g1b,
       device = "pdf",
       path = outdir,
       width = 5,
       height = 3)

###############################################################################
# Figure 3
###############################################################################

# Log-distance to the coast
# Modify Barcelona-Airport label
stations36$STANAME[10] <- "FABRA OBSERVATORY                       "
stations36$STANAME[25] <- "BCN/AEROPUERTO                          "
stations36$abb <- substr(stations36$STANAME,1,6)
stations36$b0.n <- stations36$b0.x <- as.numeric(NA)
stations36$b.trend.n <- stations36$b.trend.x <- as.numeric(NA)
# Local model
loc_frm <- as.formula(paste0("y ~ trend + Ix.lag1 + In.lag1"))
coef_n <- 2

# Extract coefficients
for(ss in 1:nrow(stations36)){
  print(stations36$STANAME[ss])

  # Update observatory records
  bin_df <- data.frame(Ix = c(t(itx3d[,,ss])),
                       In = c(t(itn3d[,,ss])),
                       t = rep(x = 1:voldim[1], each = voldim[2]),
                       l = rep(x = 1:voldim[2], times = voldim[1]))

  # Create yearly trend
  bin_df$trend <- qnorm(1/bin_df$t)
  # Remove infinite values
  bin_df <- bin_df[which(!is.infinite(bin_df$trend)),]

  # Apply GLM
  for(xx in 1:2){

    # Define dependent variable
    bin_df$y <- bin_df$Ix
    if(xx==2) bin_df$y <- bin_df$In

    # Fit model
    bin_df$Ix.lag1 <- dplyr::lag(bin_df$Ix,1,0)
    bin_df$In.lag1 <- dplyr::lag(bin_df$In,1,0)
    fit1 <- glm(formula = loc_frm, data = bin_df, family = binomial(link = "probit"))

    # Save coefficients
    if(xx==1){
      stations36$b0.x[ss] <- summary(fit1)$coef[1,1]
      stations36$b.trend.x[ss] <- summary(fit1)$coef[2,1]
    } else{
      stations36$b0.n[ss] <- summary(fit1)$coef[1,1]
      stations36$b.trend.n[ss] <- summary(fit1)$coef[2,1]
    }

  }# for xx in Ix or In
}# for each observatory ss

# Explore effects
stations36$log.dist <- log(stations36$DIST)
stations36$log.one.dist <- log1p(stations36$DIST)
# Check linear correlation effects
cor.x <- round(cor(stations36$b0.x, stations36$log.one.dist),2)
cor.n <- round(cor(stations36$b0.n, stations36$log.one.dist),2)
# Plot
xpos <- min(stations36$log.one.dist)
ypos <- min(min(stations36$b0.x),min(stations36$b0.n))
g1c <- ggplot(data = stations36, mapping = aes(x = log.one.dist)) +
  geom_point(mapping = aes(y=b0.x, color = "Tx"))  +
  geom_point(mapping = aes(y=b0.n, color = "Tn")) +
  geom_smooth(mapping = aes(y=b0.x, color = "Tx"), method = lm, se = F) +
  geom_smooth(mapping = aes(y=b0.n, color = "Tn"), method = lm, se = F) +
  theme_bw() +
  annotate("text", x=xpos+0.6, y=ypos+0.2, label= paste0("r = ",cor.x), size = 5, color="red") +
  annotate("text", x=xpos+0.6, y=ypos+0.05, label= paste0("r = ",cor.n), size = 5, color="blue") +
  ylab("intercept") +
  xlab("log(1+dist)") +
  scale_color_manual(name = "", breaks = c("Tx","Tn"),
                     values = c("Tx"="red","Tn"="blue"))
show(g1c)
# Save
ggsave(filename = "MAIN_b0_log1dist.pdf",
       plot = g1c,
       device = "pdf",
       path = outdir,
       width = 5,
       height = 3)

# Repeat for trend coefficients
cor.x <- round(cor(stations36$b.trend.x, stations36$log.one.dist),2)
cor.n <- round(cor(stations36$b.trend.n, stations36$log.one.dist),2)
# Plot
g1d <- ggplot(data = stations36, mapping = aes(x = log.one.dist)) +
  geom_point(mapping = aes(y = b.trend.x , color = "Tx"))  +
  geom_point(mapping = aes(y = b.trend.n, color = "Tn")) +
  geom_smooth(mapping = aes(y = b.trend.x, color = "Tx"), method = lm, se = F) +
  geom_smooth(mapping = aes(y = b.trend.n, color = "Tn"), method = lm, se = F) +
  theme_bw() +
  ylab("long-term trend") +
  xlab("log(1+dist)") +
  annotate("text", x=0.8, y=0.05, label= paste0("r = ",cor.x), size = 5, color="red") +
  annotate("text", x=0.83, y=-0.05, label= paste0("r =  ",cor.n), size = 5, color="blue") +
  scale_color_manual(name = "", breaks = c("Tx","Tn"),
                     values = c("Tx"="red","Tn"="blue"))
show(g1d)
# Save
ggsave(filename = "MAIN_b1_log1dist.pdf",
       plot = g1d,
       device = "pdf",
       path = outdir,
       width = 5,
       height = 3)

###############################################################################
# Figure 4
###############################################################################

# Calculate local model residuals
# Local model
loc_frm <- as.formula(paste0("y ~ trend*(Ix.lag1 + In.lag1)"))

# Create empty lists to store residuals
resid_x <- list()
resid_n <- list()

# Extract coefficients
for(ss in 1:voldim[3]){

  print(stations36$STANAME[ss])

  # Update observatory records
  bin_df <- data.frame(Ix = c(t(itx3d[,,ss])),
                       In = c(t(itn3d[,,ss])),
                       t = rep(x = 1:voldim[1], each = voldim[2]),
                       l = rep(x = 1:voldim[2], times = voldim[1]))

  # Calculate covariates
  bin_df$trend <- qnorm(1/bin_df$t)
  bin_df$Ix.lag1 <- dplyr::lag(bin_df$Ix,1,0)
  bin_df$In.lag1 <- dplyr::lag(bin_df$In,1,0)
  # Remove infinite values
  bin_df <- bin_df[which(!is.infinite(bin_df$trend)),]

  # Apply GLM
  bin_df$y <- bin_df$Ix
  resid_x[[ss]] <- resid(glm(formula = loc_frm, data = bin_df, family = binomial(link = "probit")))
  bin_df$y <- bin_df$In
  resid_n[[ss]] <- resid(glm(formula = loc_frm, data = bin_df, family = binomial(link = "probit")))

}# for each observatory ss

# Convert the residuals lists to matrices
resid_x <- do.call(cbind, resid_x)
resid_n <- do.call(cbind, resid_n)

###############################################################################
# Left plot. Relationship between residual variance

# Data.frame
plot_df <- data.frame(var.x = apply(resid_x,2,var),
                      var.n = apply(resid_n,2,var))
# Plot
plot_lim <- c(min(c(plot_df$var.x,plot_df$var.n)),
              max(c(plot_df$var.x,plot_df$var.n)))
g3a <- ggplot(data = plot_df, mapping = aes(x = var.x)) +
  geom_point(mapping = aes(y = var.n)) +
  geom_smooth(mapping = aes(y = var.n), color = "black", method = "lm", se = F) +
  geom_abline(slope = 1, color = "gray")  +
  ylim(plot_lim) +
  xlim(plot_lim) +
  xlab("Variance of Tx deviance residuals") +
  ylab("Variance of Tn deviance residuals") +
  theme_bw()
show(g3a)
ggsave(filename = "MAIN_var_resx_resn.pdf",
       plot = g3a,
       device = "pdf",
       path = outdir,
       width = 5,
       height = 3)

###############################################################################
# Rigth plot. Correlation between Tx and Tn residuals

# Store station-wise Tx and Tn correlations by intervals of 10 years
i20 <- 6
plot_df <- data.frame(cor.xn = vector(mode = "numeric",
                                      length = voldim[3]*i20),
                      t = c(rep("2-11",voldim[3]),
                            rep("12-21",voldim[3]),
                            rep("22-31",voldim[3]),
                            rep("32-41",voldim[3]),
                            rep("42-51",voldim[3]),
                            rep("52-61",voldim[3])))
# Compute correlations
for(ii in 1:i20){
  i20_ini <- 1+((ii-1)*10*voldim[2])+voldim[2]
  i20_end <- (ii*10*voldim[2])+voldim[2]
  plot_df$cor.xn[((ii-1)*voldim[3] + 1):(ii*voldim[3])] <- diag(cor(resid_x[i20_ini:i20_end,], resid_n[i20_ini:i20_end,]))
}

# Plot
plot_df$t <- factor(plot_df$t,
                    levels = c("2-11","12-21","22-31","32-41","42-51","52-61"))
g3b <- ggplot(data = plot_df, mapping = aes(x = t, y = cor.xn)) +
  geom_boxplot() +
  ylab("Correlation between Tx and Tn\ndeviance residuals") +
  theme_bw()
show(g3b)
ggsave(filename = "MAIN_cor_resx_resn_6bp.pdf",
       plot = g3b,
       device = "pdf",
       path = outdir,
       width = 5,
       height = 3)

###############################################################################
# Figure 5
###############################################################################

###############################################################################
# Left plot. Between stations correlation residuals vs. geodesic distance

# Read pairwise matrix of Euclidean distances
sfdist <- sf::st_distance(stations)
geomat <- matrix(sfdist[idx,idx], nrow = length(idx))
utri <- upper.tri(geomat)
plot_df <- data.frame(geodist = geomat[utri]/1000)
# Calculate correlations
plot_df$cor.xx <- cor(resid_x)[utri]
plot_df$cor.nn <- cor(resid_n)[utri]

# Plot
g3c <- ggplot(data = plot_df, mapping = aes(x = geodist)) +
  geom_point(mapping = aes(y=cor.xx), color = "red", alpha = 0.25)  +
  geom_point(mapping = aes(y=cor.nn), color = "blue", alpha = 0.25)  +
  geom_smooth(mapping = aes(y=cor.xx), color = "red", se = F) +
  geom_smooth(mapping = aes(y=cor.nn), color = "blue", se = F) +
  ylab("Deviance residuals correlations") +
  xlab("Euclidean distance (km)") +
  theme_bw()
show(g3c)
ggsave(filename = "MAIN_cor_res_dist.pdf",
       plot = g3c,
       device = "pdf",
       path = outdir,
       width = 5,
       height = 3)

###############################################################################
# Right plot. Between stations correlation residuals vs. SD diff.

# Set NA's for SD calculation
if(sum(tx3d==-9999)) tx3d[tx3d==-9999] <- NA
ref_idx <- c(32:61) # Reference period
tx_sd <- apply(tx3d[summer_idx,ref_idx,idx],3, sd, na.rm = TRUE)
# Calculate pairwise absolute SD differences
tx_sd_diff <- outer(tx_sd, tx_sd, FUN = function(x, y) abs(x - y))
# Store in data.frame
plot_df$sd.x <- tx_sd_diff[utri]

# Plot
g3d <- ggplot(data = plot_df, mapping = aes(x = sd.x)) +
  geom_point(mapping = aes(y=cor.xx), color = "red", alpha = 0.25)  +
  geom_point(mapping = aes(y=cor.nn), color = "blue", alpha = 0.25)  +
  geom_smooth(mapping = aes(y=cor.xx), color = "red", se = F) +
  geom_smooth(mapping = aes(y=cor.nn), color = "blue", se = F) +
  ylab("Deviance residuals correlations") +
  xlab("Absolute differences in s.d. of Tx") +
  theme_bw()
show(g3d)
ggsave(filename = "MAIN_cor_res_SDdiff.pdf",
       plot = g3d,
       device = "pdf",
       path = outdir,
       width = 5,
       height = 3)

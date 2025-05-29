#######################################
### Section 2 - DOWNLOAD Tx Tn DATA ###
#######################################

### Script provided by:
###   Zeus Gracia-Tabuenca

### DATA ACCESSED FOR THIS WORK ON
### 12-MAY-2024

### Note: ECA makes periodic updates,
###       so the data used in this work
###       and those accessed at another
###       time may differ


# Clear workspace
rm(list = ls())

# Load library
if(!is.element("data.table", row.names(installed.packages()))) install.packages("data.table")
library("data.table")

# Read list of observatories of interest. Based on previous work.
load(file.path("data", "stations.rda"))

###############################################################################
# Download
###############################################################################

# Download Tmax dataset from ECA
# Increase timeout time
getOption("timeout")
options(timeout=600)

# Follow: https://stackoverflow.com/questions/3053833/using-r-to-download-zipped-data-file-extract-and-import-data
# Space for 650.5MB will be required (2024-05-08)
tmp_dir <- tempdir()
if(!file.exists(tmp_dir)) dir.create(path = tmp_dir, recursive = T)
tmp_file <- file.path(tmp_dir,"ECA_blend_tx.zip")
download.file(url = "https://knmi-ecad-assets-prd.s3.amazonaws.com/download/ECA_blend_tx.zip",
              destfile = tmp_file)

# Get list of zipped files
# https://stackoverflow.com/a/32871141
zip_ls <- unzip(tmp_file, list=TRUE)$Name
head(zip_ls)
# Get number ID
num_ls <- grep("TX_STAID", zip_ls)
num_id <- sapply(num_ls, function(x) unlist(strsplit(zip_ls[x],"TX_STAID")[[1]][2]))
num_id <- as.integer(sapply(1:length(num_id), function(x) unlist(strsplit(num_id[x],".txt")[[1]][1])))

# Find the ones already in the stations data.frame and copy
raw_dir <- file.path("data", "raw", "tx")
if(!dir.exists(raw_dir)) dir.create(path = raw_dir, recursive = T)
# Extract files of interest
# Info file
unzip(zipfile = tmp_file,
      files = zip_ls[num_ls[match(stations$STAID, num_id)]],
      exdir = raw_dir)

# Delete temporary directory
unlink(tmp_dir, recursive = T)

###############################################################################
# Same for Minimum temperature

tmp_dir <- tempdir()
if(!file.exists(tmp_dir)) dir.create(path = tmp_dir, recursive = T)
tmp_file <- file.path(tmp_dir,"ECA_blend_tn.zip")
download.file(url = "https://knmi-ecad-assets-prd.s3.amazonaws.com/download/ECA_blend_tn.zip",
              destfile = tmp_file)

# Get list of zipped files
# https://stackoverflow.com/a/32871141
zip_ls <- unzip(tmp_file, list=TRUE)$Name
head(zip_ls)
# Get number ID
num_ls <- grep("TN_STAID", zip_ls)
num_id <- sapply(num_ls, function(x) unlist(strsplit(zip_ls[x],"TN_STAID")[[1]][2]))
num_id <- as.integer(sapply(1:length(num_id), function(x) unlist(strsplit(num_id[x],".txt")[[1]][1])))

# Find the ones already in the stations data.frame and copy
raw_dir <- file.path("data", "raw", "tn")
if(!dir.exists(raw_dir)) dir.create(path = raw_dir, recursive = T)
# Extract files of interest
# Info file
unzip(zipfile = tmp_file,
      files = zip_ls[num_ls[match(stations$STAID, num_id)]],
      exdir = raw_dir)

# Delete temporary directory
unlink(tmp_dir, recursive = T)

###############################################################################
# Create Tx and Tn matrices
###############################################################################

#Extract all the files in the folder Iberia_TX from that correspond to the data from each station
tx_ls <- list.files(path = file.path("data", "raw", "tx"),
                    full.names = TRUE)

# Create empty data.frame
# Daily date range
target_date <- seq(from = as.Date("1960-01-01"),
                   to = as.Date("2023-12-31"),
                   by = "day")
tail(target_date)
target_date <- target_date[grep("-02-29", target_date,invert = T)]
# Output file
tx_mat <- as.data.frame(matrix(data = as.numeric(NA),
                               nrow = length(target_date),
                               ncol = nrow(stations)))
tx_mat <- cbind(target_date, tx_mat)
names(tx_mat) <- c("Date", stations$STAID)

# Get data
cat("/n")
for(ii in 1:length(tx_ls)){
  cat(paste0("..",ii))
  # Read file
  aux <- read.csv(file = tx_ls[ii], skip = 19)
  # Remove scores without good label data (0)
  nozero <- which(aux$Q_TX != 0 | aux$TX == -9999)
  if(length(nozero)>0) aux$TX[nozero] <- NA
  # Take the days of interest
  aux$DATE <- as.Date(strptime(aux$DATE,"%Y%m%d"))
  tx_mat[,ii+1] <- aux$TX[match(tx_mat$Date, aux$DATE)]
}
cat("\n")
#Save matrix
fwrite(x = tx_mat, file = file.path("data", "Tx_mat.csv"))

###############################################################################
# Same for Minimum temperature

#Extract all the files in the folder Iberia_TX from that correspond to the data from each station
tn_ls <- list.files(path = file.path("data", "raw", "tn"),
                    full.names = TRUE)

# Create empty data.frame
# Output file
tn_mat <- as.data.frame(matrix(data = as.numeric(NA),
                               nrow = length(target_date),
                               ncol = nrow(stations)))
tn_mat <- cbind(target_date, tn_mat)
names(tn_mat) <- c("Date", stations$STAID)

# Get data
cat("/n")
for(ii in 1:length(tn_ls)){
  cat(paste0("..",ii))
  # Read file
  aux <- read.csv(file = tn_ls[ii], skip = 19)
  # Remove scores without good label data (0)
  nozero <- which(aux$Q_TN != 0 | aux$TN == -9999)
  if(length(nozero)>0) aux$TN[nozero] <- NA
  # Take the days of interest
  aux$DATE <- as.Date(strptime(aux$DATE,"%Y%m%d"))
  tn_mat[,ii+1] <- aux$TN[match(tn_mat$Date, aux$DATE)]
}
cat("\n")
#Save matrix
fwrite(x = tn_mat, file = file.path("data", "Tn_mat.csv"))

# Remove raw directories
unlink(file.path("data", "raw"), recursive = T)

# Restart R
.rs.restartR()

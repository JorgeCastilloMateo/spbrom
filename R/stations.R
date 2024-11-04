#' @title Weather stations
#' 
#' @description A data set containing information about the weather stations in 
#'   peninsular Spain. The variables are the following:
#'
#' \itemize{
#'   \item STAID  : Station identifier
#'   \item STANAME: Station name
#'   \item CN     : Country code (ISO3116 country codes)
#'   \item LAT    : Latitude in degrees
#'   \item LON    : Longitude in degrees
#'   \item HGTH   : Station elevation in meters
#' }
#'   
#' @docType data
#' @keywords datasets
#' @name stations
#' @source \href{https://www.ecad.eu}{EUROPEAN CLIMATE ASSESSMENT & DATASET (ECA&D)}
#' @references 
#' Klein Tank AMG and Coauthors (2002). 
#' Daily Dataset of 20th-Century Surface Air Temperature and Precipitation 
#' Series for the European Climate Assessment.
#' \emph{International Journal of Climatology}, \strong{22}(12), 1441-1453.
#' \doi{10.1002/joc.773}.
#' 
#' @usage 
#' data("stations")
#' @format 
#' \code{stations} is a data frame with 40 rows and 6 variables.
#' 
"stations"

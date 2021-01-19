#' Gap-fills and z-normalizes the data
#' @description This function gap-fills and z-normalizes the metabolism and
#' associated driver data.
#'
#' @param ts The individual site_year
#' 
#' @param PQ the photosynthetic quotient used to convert aquatic metabolism estimates 
#' to gC m-1 d-1
#' 
#' @param block The size of the largest block of missing data that the function will fill-in. 
#' Parameter for fillMiss3.R
#' 
#' @param pmiss The maximum amount of the missing data that can be missing in the 
#' dataset for fill-in procedure to be performed. Parameter for fillMiss3.R
#' 
#' @return Returns gap-filled and z-normalized time series
#'
#' @export
#===============================================================================
#Function for gap filling and z-normalizing time series. Adapted from the
#SUI_z_norm_func and SFS_z_norm_func code
#Created 8/22/2019
#===============================================================================
gapfill_norm <- function(ts, PQ, block, pmiss){
  #-------------------------------------------------
  #Replacing negative GPP and positive ER values
  #-------------------------------------------------
    #If there is negative GPP, then replace with NA
      if(length(ts[ts[, "GPP"] < 0 & !is.na(ts[, "GPP"]), "GPP"]) !=0){
        ts[ts[, "GPP"] < 0 & !is.na(ts[, "GPP"]), "GPP"] <- NA
      } #End if statement
  
    #If there is positive ER, then replace with NA
      if(length(ts[ts[, "ER"] > 0 & !is.na(ts[, "ER"]), "ER"]) !=0){
        ts[ts[, "ER"] > 0 & !is.na(ts[, "ER"]), "ER"] <- NA
      } #End if statement
      
    #Calculating NEP 
      ts$NEP <- ts[, "GPP"] + ts[, "ER"]  
  
  #-------------------------------------------------
  #Gap-filling the data
  #-------------------------------------------------
    #GPP
      ts$GPP_filled <- fillMiss3(ts, "GPP", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = TRUE)

    #ER
      ts$ER_filled <- fillMiss3(ts, "ER", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)

    #NEP
      ts$NEP_filled <- fillMiss3(ts, "NEP", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)

    #Water temperature
      ts$Wtemp_filled <- fillMiss3(ts, "temp_water", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)

    #PAR
      ts$PAR_filled <- fillMiss3(ts, "PAR_sum", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)

    #Discharge
      ts$Disch_filled  <- fillMiss3(ts, "discharge", block = block, pmiss = pmiss,
        model = "trend", smooth = TRUE, log = "y", plot = FALSE)

  #-------------------------------------------------
  #Z-normalization of the data
  #-------------------------------------------------
    #Defining a function for z normalization
      znorm <- function(timeseries, var){
        ts.mean <- mean(timeseries[, var], na.rm = TRUE)
        ts.dev <- sd(timeseries[, var], na.rm = TRUE)
        (timeseries[, var] - ts.mean)/ts.dev
      }

    #Cycling through and applying the Z normalization function to each year
      ts$GPP_norm <- znorm(ts, "GPP_filled")
      ts$ER_norm <- znorm(ts, "ER_filled")
      ts$NEP_norm <- znorm(ts, "NEP_filled")
      ts$Wtemp_norm <- znorm(ts, "Wtemp_filled")
      ts$PAR_norm <- znorm(ts, "PAR_filled")

  #-------------------------------------------------
  #Converting GPP from (gO2 m-2 d-1) to (gC m-2 d-1)
  #-------------------------------------------------
    ts$GPP_C <- ts[, "GPP"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$ER_C <- ts[, "ER"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$NEP_C <- ts[, "NEP"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$GPP_C_filled <- ts[, "GPP_filled"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$ER_C_filled <- ts[, "ER_filled"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    ts$NEP_C_filled <- ts[, "NEP_filled"] * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107

  #-------------------------------------------------
  #Writing the final output
  #-------------------------------------------------
    final <- ts[, c("Date", "U_ID", "Year", "DOY", "GPP", "ER", "NEP", "GPP_C",
      "ER_C", "NEP_C", "K600", "DO_obs", "DO_sat", "temp_water", "discharge",
      "PAR_sum", "Stream_PAR_sum", "LAI_proc", "GPP_filled", "ER_filled",
      "NEP_filled", "GPP_C_filled", "ER_C_filled", "NEP_C_filled", "Wtemp_filled",
      "Disch_filled", "PAR_filled", "GPP_norm", "ER_norm", "NEP_norm","Wtemp_norm", 
      "PAR_norm")]

} #End gapfill_norm function


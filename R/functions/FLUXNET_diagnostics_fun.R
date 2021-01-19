#' Calculates model performance diagnostics for FLUXNET estimates
#' @description This funciton calculates several diagnostic measures for FLUXNET data.
#'
#' @param ts The file with standardized FLUXNET data
#' @return Returns a table of model diagnostics
#'
#' @export
#'
#===============================================================================
#Calculating model diagnostics for the FLUXNET data as part of the metabolism synthesis 
#Created 2/6/2020
#===============================================================================
#Function for calculating model diagnostics 
  FLUXNET_diagnostics_fun <- function(ts){
    #Percentage of days with negative GPP
      GPP_neg <- (length(ts[ts[, "GPP"] < 0 & !is.na(ts[, "GPP"]), "GPP"]) / nrow(ts)) * 100

    #Percentage of days with positive ER
      ER_pos <- (length(ts[ts[, "ER"] > 0 & !is.na(ts[, "ER"]), "ER"]) / nrow(ts)) * 100

    #How many days does this year have (excluding neg GPP and Positive ER)
      num_days <- length(ts[ts[, "GPP"] >= 0 & !is.na(ts[, "GPP"]) & ts[, "ER"] <= 0 & !is.na(ts[, "ER"]), "GPP"])

    #Get the final set of diagnostic data
      final <- data.frame(GPP_neg, ER_pos, num_days)

    return(final)
      
  } #End FLUXNET_diagnostics_fun

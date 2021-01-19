#' Calculates model performance diagnostics for metabolism estimates
#' @description This funciton calculates several diagnostic measures for modeled
#' metabolism data.
#'
#' @param ts The file with standardized metabolism data
#' @return Returns a table of model diagnostics
#'
#' @export
#'
#===============================================================================
#Calculating model diagnostics for the metabolism synthesis 
#Created 7/9/2019
#===============================================================================
#Function for calculating model diagnostics 
  diagnostics_fun <- function(ts){
    #Correlation between ER and K
      #er_k <- cor(Year[, "ER"], Year[, "K600"], use = "pairwise.complete.obs")
      ER_K <- summary(lm(ts[, "ER"] ~ ts[, "K600"],
        na.action = na.omit))$adj.r.squared

    #What is the maximum daily K600
      K600_max <- mean(ts[, "K600"], na.rm = TRUE)

    #Percentage of days with negative GPP
      GPP_neg <- (length(ts[ts[, "GPP"] < 0 & !is.na(ts[, "GPP"]), "GPP"]) / nrow(ts)) * 100

    #Percentage of days with positive ER
      ER_pos <- (length(ts[ts[, "ER"] > 0 & !is.na(ts[, "ER"]), "ER"]) / nrow(ts)) * 100

    #How many days does this year have (excluding neg GPP and Positive ER)
      num_days <- length(ts[ts[, "GPP"] >= 0 & !is.na(ts[, "GPP"]) & ts[, "ER"] <= 0 & !is.na(ts[, "ER"]), "GPP"])

    #Get the final set of diagnostic data
      final <- data.frame(ER_K, GPP_neg, ER_pos, num_days, K600_max)

    return(final)
      
  } #End diagnostics_fun

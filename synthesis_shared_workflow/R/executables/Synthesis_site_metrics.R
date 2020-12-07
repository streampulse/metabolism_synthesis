#===============================================================================
#Calculate a series of metrics for the lotic synthesis sites including
#mean annual GPP, 95th percentile of GPP, mean annual ER, 5th percentile of ER,
#mean annual cumulative PAR (both incoming and estimated), magnificent 7
#statistics for GPP, ER, water temperature, discharge, PAR, and LAI 
#Created 9/14/2020
#===============================================================================
Synthesis_site_metrics <- function(){
#Reading in my generic magnificent 7 functions
  source(here::here("R", "functions", "magnificent_seven.R"))
        
#-------------------------------------------------
#Function to aid applying the mag7_fun to multiple variables
#-------------------------------------------------
  synthesis_mag7 <- function(compiled_data, var, standardize){
    #Calculating magnificent 7 with error catch for columns with missing data
      if(all(is.na(compiled_data[, var])) == FALSE){
        calc_mag7 <- mag7_fun(compiled_data, var = var, standardize = standardize)
      } #End if statement
      
    #For variables that have no data
      if(all(is.na(compiled_data[, var])) == TRUE){
        calc_mag7 <- setNames(data.frame(t(rep(NA, 7))), c("lam1", "tau2", 
          "tau3", "tau4", "amplitude", "phase", "ar1"))
      } #End if statement

    return(calc_mag7)      
  
  } #End synthesis_mag7 function

#-------------------------------------------------
#Define a function to calculate all site metrics
#-------------------------------------------------
  synthesis_summary_calc <- function(Site_ID, gap_filled){
    #Get the site of interest
      SOI <- gap_filled[[Site_ID]]

    #Mean annual cumulative GPP
      if(all(is.na(SOI[, "GPP_C_filled"])) == FALSE){
        gpp_C_sum <- mean(aggregate(GPP_C_filled ~ Year, data = SOI, 
          FUN = sum, na.rm = TRUE)$GPP_C_filled)
      } #End if statement
  
      if(all(is.na(SOI[, "GPP_C_filled"])) == TRUE){gpp_C_sum <- NA}
      
    #P95 of GPP
      if(all(is.na(SOI[, "GPP_C"])) == FALSE){
        gpp_C_95 <- quantile(SOI[, "GPP_C"], 0.95, na.rm = TRUE)[[1]]
      } else {
        gpp_C_95 <- NA
      } #End if else statement      
      
    #Mean annual cumulative ER
      if(all(is.na(SOI[, "ER_C_filled"])) == FALSE){
        er_C_sum <- mean(aggregate(ER_C_filled ~ Year, data = SOI, 
          FUN = sum, na.rm = TRUE)$ER_C_filled)
      } #End if statement
      
      if(all(is.na(SOI[, "ER_C_filled"])) == TRUE){er_C_sum <- NA}   
      
    #P05 of ER
      if(all(is.na(SOI[, "ER_C"])) == FALSE){
        er_C_05 <- quantile(SOI[, "ER_C"], 0.05, na.rm = TRUE)[[1]]
      } else {
        er_C_05 <- NA
      } #End if else statement       
      
    #Mean annual cumulative incoming PAR (kmol m-2 y-1)
      if(all(is.na(SOI[, "PAR_sum"])) == FALSE){
        PAR_sum <- mean(aggregate(PAR_sum ~ Year, data = SOI, 
          FUN = sum, na.rm = TRUE)$PAR_sum) * 0.001
      } else {
        PAR_sum <- NA
      } #End if else statement
  
    #Mean annual cumulative PAR at the stream surface (kmol m-2 y-1)
      if(all(is.na(SOI[, "Stream_PAR_sum"])) == FALSE){
        Stream_PAR_sum <- mean(aggregate(Stream_PAR_sum ~ Year, data = SOI, 
          FUN = sum, na.rm = TRUE)$Stream_PAR_sum) * 0.001
      } else {
        Stream_PAR_sum <- NA
      } #End if else statement  
      
    #List of column names to add
      mag7_colnames <- c("_mean", "_cv", "_skew", "_kurt", "_amp", "_phase", 
        "_ar1")
  
    #GPP
      mag7_gpp <- synthesis_mag7(SOI, var = "GPP_C_filled", standardize = "yes")
        colnames(mag7_gpp) <- paste("gpp_C", mag7_colnames, sep = "")    
  
    #ER
      mag7_er <- synthesis_mag7(SOI, var = "ER_C_filled", standardize = "yes")
        colnames(mag7_er) <- paste("er_C", mag7_colnames, sep = "")  
      
    #Water temperature
      mag7_wtemp <- synthesis_mag7(SOI, var = "Wtemp_filled", standardize = "yes")
        colnames(mag7_wtemp) <- paste("Wtemp", mag7_colnames, sep = "") 
        
    #Discharge
      mag7_disch <- synthesis_mag7(SOI, var = "Disch_filled", standardize = "yes")
        colnames(mag7_disch) <- paste("Disch", mag7_colnames, sep = "")           
      
    #PAR
      mag7_PAR <- synthesis_mag7(SOI, var = "PAR_filled", standardize = "yes")
        colnames(mag7_PAR) <- paste("PAR", mag7_colnames, sep = "")  
        
    #LAI
      mag7_LAI <- synthesis_mag7(SOI, var = "LAI_proc", standardize = "yes")
        colnames(mag7_LAI) <- paste("LAI", mag7_colnames, sep = "") 
        
    #Bind them all together
      bound <- cbind(mag7_gpp, mag7_er, mag7_wtemp, mag7_disch, mag7_PAR, mag7_LAI) 
  
      final <- data.frame(Site_ID, gpp_C_sum, gpp_C_95, er_C_sum, er_C_05, PAR_sum, 
        Stream_PAR_sum, bound, stringsAsFactors = FALSE)
      
        colnames(final)[1:5] <- c("Site_ID", "ann_GPP_C", "upper_GPP_C", "ann_ER_C", "lower_ER_C")
        
    return(final)
        
  } #End site_wrapper function

#-------------------------------------------------
#Calculate the site metrics, bind with MODIS NPP and final formatting   
#-------------------------------------------------
  #Read in the gap-filled data
    gap_filled <- readRDS(here::here("output", "synthesis_gap_filled.rds"))
  
  #Calculate the metrics
    metrics_compiled <- plyr::ldply(
      names(gap_filled), 
      synthesis_summary_calc, 
      gap_filled = gap_filled
    )
  
  #Read in the annual MODIS NPP (g C m-2 y-1)
    MODIS_annual_NPP <- readRDS(here::here("output", "MODIS_annual_NPP.rds"))
    
  #Calculate the mean annual MODIS NPP (g C m-2 y-1)  
    calc_mean_NPP <- function(Site_ID, gap_filled, npp_data){
      #Get the site of interest
        years <- unique((gap_filled[[Site_ID]][, "Year"]))
        
      #Subset for the site and filtered years
        mod_sub <- npp_data[npp_data[, "Site_ID"] %in% Site_ID & npp_data[, "Year"] %in% years, ]
        
      #Get the final output
        final <- setNames(data.frame(
          Site_ID,
          mean(mod_sub[, "MOD_ann_NPP"], na.rm = TRUE)
        ), c("Site_ID", "MOD_ann_NPP"))
        
      return(final)
       
    } #End calc_mean_NPP function
    
    mean_annual_NPP <- plyr::ldply(
      names(gap_filled),
      calc_mean_NPP,
      gap_filled = gap_filled,
      npp_data = MODIS_annual_NPP
    )
   
  #Merge the NPP data with the other site metrics
    final <- merge(metrics_compiled, mean_annual_NPP, by = "Site_ID", all.x = TRUE)
    
  #Export the output
    saveRDS(final, here::here("output", "synthesis_site_metrics.rds"))  
    write.csv(final, here::here("output", "synthesis_site_metrics.csv"), row.names = FALSE, quote = FALSE)  
  
} #End Synthesis_site_metrics wrapper

Synthesis_site_metrics()  
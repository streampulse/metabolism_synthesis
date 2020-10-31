#===============================================================================
#Calculate summary metrics for FLUXNET sites
#Created 2/6/2020
#===============================================================================
library("here")

FLUXNET_site_metrics <- function(){
#Read in the compiled annual FLUXNET data
  annual_compiled <- readRDS(here("output", "FLUXNET_annual_compiled.rds"))

#Read in the filtered FLUXNET data
  # filtered <- readRDS("C:/research/postdoc research/metab_synthesis/output/FLUXNET_filtered.rds")
  filtered <- readRDS(here("output", "FLUXNET_filtered.rds"))

#-------------------------------------------------
#Define a function to calculate metrics for FLUXNET sites
#-------------------------------------------------
  FLUXNET_calc <- function(Site_ID){
    #Get annual estimates for the site of interest
      SOI_annual <- annual_compiled[[Site_ID]]
      
    #Calculate mean annual rates
      mean_annual <- colMeans(SOI_annual[, c("GPP", "ER", "Net")], na.rm = TRUE)
      
    #Get the daily data for the site of interest
      SOI_daily <- filtered[[Site_ID]]  
      
    #P95 of GPP
      gpp_95 <- quantile(SOI_daily[, "GPP"], 0.95, na.rm = TRUE)[[1]]
      
    #P05 of ER
      er_05 <- quantile(SOI_daily[, "ER"], 0.05, na.rm = TRUE)[[1]]
      
    #Bind them all together
      bound <- data.frame(c(
        data.frame(Site_ID, stringsAsFactors = FALSE), 
        mean_annual, 
        data.frame(gpp_95, stringsAsFactors = FALSE), 
        data.frame(er_05, stringsAsFactors = FALSE)
      ), stringsAsFactors = FALSE)
      
      final <- bound[, c("Site_ID", "GPP", "gpp_95", "ER", "er_05")]
      
      colnames(final)[1:5] <- c("Site_ID", "ann_GPP", "upper_GPP", "ann_ER", "lower_ER")     
      
    return(final)
      
  } #End FLUXNET_calc

#-------------------------------------------------
#Calculate the metrics  
#-------------------------------------------------
  #Calculate the summary metrics  
    metrics_compiled <- plyr::ldply(names(filtered), FLUXNET_calc)
    
  #Export the output
    saveRDS(metrics_compiled, here("output", "FLUXNET_site_metrics.rds"))
  
} #End FLUXNET_site_metrics wrapper

FLUXNET_site_metrics()

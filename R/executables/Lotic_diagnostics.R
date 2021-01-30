#===============================================================================
#Calculate metabolism diagnostics for each site-year of data
#Created 1/22/2020
#===============================================================================
lotic_diagnostics <- function(){
#-------------------------------------------------
#Define functions to use  
#-------------------------------------------------  
  #Read in a function for calculating diagnostics
    source(here::here("R", "functions", "diagnostics_fun.R"))
 
  #Function to see if ancillary data is available for each site-year
    yearly_ancillary <- function(Site_ID, ts){
      #See if various time series data are available for the site-year
        ifelse(all(is.na(ts[, "PAR_sum"])) == FALSE, NLDAS_data <- TRUE, NLDAS_data <- FALSE)
        ifelse(all(is.na(ts[, "LAI_proc"])) == FALSE, LAI_data <- TRUE, LAI_data <- FALSE)
        ifelse(all(is.na(ts[, "Stream_PAR_sum"])) == FALSE, Predicted_light <- TRUE, Predicted_light <- FALSE)  
        
      #Year of interest
        YOI <- unique(ts[, "Year"])
        
      #See if there is annual MODIS NPP data for the year
        mod_check <- MODIS_NPP[MODIS_NPP[, "Site_ID"] == Site_ID & MODIS_NPP[, "Year"] == YOI, ]
        ifelse(nrow(mod_check) != 0, MODIS_NPP <- TRUE, MODIS_NPP <- FALSE)
        
      #Return the final information
        final <- data.frame(Site_ID, YOI, NLDAS_data, LAI_data, Predicted_light, 
          MODIS_NPP, stringsAsFactors = FALSE)
        
          colnames(final)[2] <- "Year"
          
        return(final)
        
    } #End yearly_ancillary function

  #Define a wrapper function to calculate diagnostics for each site-year
    yearly_diagnostics <- function(Site_ID){
      #Get the timeseries for the site of interest
        SOI <- synthesis_standardized[[Site_ID]]    
      
      #Split by year
        year_split <- split(SOI, SOI[, "Year"])
   
      #Calculate diagnostics for each site-year of data
        sy_diag <- do.call(rbind, lapply(year_split, diagnostics_fun))   
        
      #Add site-year and site name
        sy_diag$Year <- as.numeric(row.names(sy_diag))
        sy_diag$Site_ID <- Site_ID
        
      #Check for availability of ancillary data for each site-year of data
        sy_ancillary <- do.call(rbind,lapply(year_split, FUN = yearly_ancillary, 
          Site_ID = Site_ID))       
        
      #Merge the tables together
        merged <- merge(sy_diag, sy_ancillary, by = c("Site_ID", "Year"))
        
      #Reorder for the final output
        final <- merged[, c("Site_ID", "Year", "ER_K", "K600_max", "GPP_neg", "ER_pos", "num_days", 
          "NLDAS_data", "LAI_data", "Predicted_light", "MODIS_NPP")]
        
      return(final)
        
    } #End yearly_diagnostics function
  
#-------------------------------------------------
#Calculate diagnostics and ancillary data availability  
#-------------------------------------------------    
  #Read in the standardized data
    if(!exists("synthesis_standardized")){synthesis_standardized <- readRDS(here::here("output", "lotic_standardized_metabolism.rds"))}
  
  #Read in the MODIS annual NPP
    MODIS_NPP <- readRDS(here::here("output", "MODIS_annual_NPP.rds"))   
    
  #Apply the function to calculate diagnostics for each site-year
    compiled_diagnostics <- do.call(rbind, lapply(names(synthesis_standardized), yearly_diagnostics))
  
  #Save the output
    saveRDS(compiled_diagnostics, here::here("output", "lotic_yearly_diagnostics.rds"))
  
} #End Synthesis_diagnostics wrapper function

lotic_diagnostics()
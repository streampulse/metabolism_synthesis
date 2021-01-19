#===============================================================================
#Calculate diagnostics for each site-year of FLUXNET data
#Created 2/6/2020
#===============================================================================
library("here")

FLUXNET_diagnostics <- function(){
#Read in a function for calculating diagnostics
  if(exists("FLUXNET_diagnostics_fun", mode = "function") == FALSE){
    source(here::here("R", "functions", "FLUXNET_diagnostics_fun.R"))
  } #End if statement

#Read in the standardized data
  if(!exists("fluxnet_standardized")){fluxnet_standardized <- readRDS(here::here("output", "fluxnet_standardized.rds"))}

#Define a wrapper function to calculate diagnostics for each site-year
  yearly_diagnostics <- function(Site_ID){
    #Get the timeseries for the site of interest
      SOI <- fluxnet_standardized[[Site_ID]]    
    
    #Split by year
      year_split <- split(SOI, SOI[, "Year"])
 
    #Calculate diagnostics for each site-year of data
      sy_diag <- do.call(rbind, lapply(year_split, FLUXNET_diagnostics_fun))   

    #Add site-year and site name
      sy_diag$Year <- as.numeric(row.names(sy_diag))
      sy_diag$Site_ID <- Site_ID
      
    #Reorder for the final output
      sy_diag[, c("Site_ID", "Year", "GPP_neg", "ER_pos", "num_days")]
      
  } #End yearly_diagnostics function
  
#Apply the function to calculate diagnostics for each site-year
  compiled_diagnostics <- do.call(rbind, lapply(names(fluxnet_standardized), yearly_diagnostics))

#Save the output
  saveRDS(compiled_diagnostics, here::here("output", "FLUXNET_yearly_diagnostics.rds"))
  
} #End FLUXNET_diagnostics wrapper function

FLUXNET_diagnostics()
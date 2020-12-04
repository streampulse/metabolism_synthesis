#===============================================================================
#Metabolism synthesis shared workflow
#Functions all use relative paths, the here package means you can place the
#file structure anywhere on your comp and they will work
#Updated 9/14/2020
#===============================================================================
library("here")

#List of packages you might need if you don't have them (everything
#except "here" is referenced with :: so shouldn't need to load them)
#c("here", "plyr", "stringr", "dplyr", "lmomco")
  
#Create a here file in the project location (first time only)
  if(file.exists(".here") == FALSE){here::set_here()}

#-------------------------------------------------
#1a-b. Filter and gap-fill the lotic metabolism dataset
#-------------------------------------------------  
  #Read in the function to filter the standardized dataset
    source(here("R", "functions", "filter_metab.R"))
      
  #Filter the standardized dataset
    synthesis_filtered <- filter_metab(
      diag = readRDS(here("output", "lotic_yearly_diagnostics.rds")), 
      filters = paste0("num_days >= ", 365 * 0.6), 
      metab_rds = readRDS(here("output", "lotic_standardized_metabolism.rds")),
      has_NLDAS = TRUE, 
      has_LAI = TRUE, 
      has_StreamLight = TRUE, 
      has_MODIS_NPP = TRUE
    )
      
  #Save the output
    saveRDS(synthesis_filtered, here("output", "synthesis_filtered.rds"))

  #Gap-fill the data
    #Read in a suite of gap-filling functions
      source(here("R", "functions", "fillMiss3.R"))
      source(here("R", "functions", "gapfill_norm.R"))
      source(here("R", "functions", "synthesis_gapfill.R"))
        
      gap_filled <- synthesis_gapfill(
        synthesis_filtered = readRDS(here("output", "synthesis_filtered.rds")), 
        PQ = 1.25,
        block = 150,
        pmiss = 50
      )    
        
    #Save the output
      saveRDS(gap_filled, here("output", "synthesis_gap_filled.rds"))
      
#-------------------------------------------------
#1c. Filter FLUXNET sites
#------------------------------------------------- 
  #Filter the standardized dataset
    FLUXNET_filtered <- filter_metab(
      diag = readRDS(here("output", "FLUXNET_yearly_diagnostics.rds")), 
      filters = paste0("num_days >= ", 365 * 0.6), 
      metab_rds = readRDS(here("output", "FLUXNET_standardized.rds"))
    )  
      
  #Remove FLUXNET tier-2 sites, that grant data producers collaboration rights
    FLUXNET_filtered_final <- FLUXNET_filtered[!names(FLUXNET_filtered) %in% c("RU-Sam", "RU-SkP", "RU-Tks", "RU-Vrk", "SE-St1", "ZA-Kru")]
        
  #Save the output
    saveRDS(FLUXNET_filtered_final, here("output", "FLUXNET_filtered.rds"))      

#-------------------------------------------------
#2. Calculate a series of data metrics for each site     
#-------------------------------------------------
  #2a. Calculate site metrics for synthesis sites
    source(here("R", "executables", "Synthesis_site_metrics.R"))
      
  #2b. Calculate site metrics for FLUXNET sites
    source(here("R", "executables", "FLUXNET_site_metrics.R"))

#===============================================================================
#Metabolism synthesis shared workflow
#Functions all use relative paths, the here package means you can place the
#file structure anywhere on your comp and they will work
#Updated 9/14/2020
#===============================================================================

#List of packages you might need if you don't have them (everything
#is referenced with :: so shouldn't need to load them)
#c("here", "plyr", "stringr", "dplyr", "lmomco")
  
# #Create a here file in the project location (first time only)
#   if(file.exists(".here") == FALSE){here::set_here()}

#===============================================================================
#II. Filtering, gap-filling, and calculating metrics
#===============================================================================
#-------------------------------------------------
#1a-b. Calculate diagnostics for both lotic and terrestrial sites      
#-------------------------------------------------      
  #Calculate diagnostics for each site-year of lotic data as well as ancillary
  #data availability*
    source(here::here("R", "executables", "Lotic_diagnostics.R")) 
        
  #Calculate diagnostics for each site-year of FLUXNET data*
    source(here::here("R", "executables", "FLUXNET_diagnostics.R"))    

#-------------------------------------------------
#2a-b. Filter and gap-fill the lotic metabolism dataset
#-------------------------------------------------  
  #Read in the function to filter the standardized dataset
    source(here::here("R", "functions", "filter_metab.R"))
      
  #Filter the standardized dataset
    synthesis_filtered <- filter_metab(
      diag = readRDS(here::here("output", "lotic_yearly_diagnostics.rds")), 
      filters = paste0(paste0("num_days >= ", 365 * 0.6), " & ER_K < 0.6", " & K600_max < 100"),
      metab_rds = readRDS(here::here("output", "lotic_standardized_metabolism.rds")),
      has_NLDAS = TRUE, 
      has_MODIS_NPP = TRUE
    )
      
  #Save the output
    saveRDS(synthesis_filtered, here::here("output", "synthesis_filtered.rds"))

  #Gap-fill the data
    #Read in a suite of gap-filling functions
      source(here::here("R", "functions", "fillMiss3.R"))
      source(here::here("R", "functions", "gapfill_norm.R"))
      source(here::here("R", "functions", "synthesis_gapfill.R"))
        
      gap_filled <- synthesis_gapfill(
        synthesis_filtered = readRDS(here::here("output", "synthesis_filtered.rds")), 
        PQ = 1.25,
        block = 150,
        pmiss = 50
      )    
        
    #Save the output
      saveRDS(gap_filled, here::here("output", "synthesis_gap_filled.rds"))
      
#-------------------------------------------------
#2c. Filter FLUXNET sites
#------------------------------------------------- 
  #Filter the standardized dataset
    FLUXNET_filtered <- filter_metab(
      diag = readRDS(here::here("output", "FLUXNET_yearly_diagnostics.rds")), 
      filters = paste0("num_days >= ", 365 * 0.6), 
      metab_rds = readRDS(here::here("output", "FLUXNET_standardized.rds"))
    )  
      
  #Remove FLUXNET tier-2 sites, that grant data producers collaboration rights
    FLUXNET_filtered_final <- FLUXNET_filtered[!names(FLUXNET_filtered) %in% c("RU-Sam", "RU-SkP", "RU-Tks", "RU-Vrk", "SE-St1", "ZA-Kru")]
        
  #Save the output
    saveRDS(FLUXNET_filtered_final, here::here("output", "FLUXNET_filtered.rds"))      

#-------------------------------------------------
#3. Calculate a series of data metrics for each site     
#-------------------------------------------------
  #3a. Calculate site metrics for synthesis sites
    source(here::here("R", "executables", "Synthesis_site_metrics.R"))
      
  #3b. Calculate site metrics for FLUXNET sites
    source(here::here("R", "executables", "FLUXNET_site_metrics.R"))
    

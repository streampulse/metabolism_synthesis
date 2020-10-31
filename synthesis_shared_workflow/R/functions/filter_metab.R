#' Filters metabolism data based on selection criteria
#' @description This funciton filters the combined metabolism data based on
#' the calculated diagnostics for each site-year.
#'
#' @param diag The site-year diagnostics calculated from diagnostics_fun
#' @param filters The filtering criteria. For example "num_days >= 180". Multiple
#' criteria can also be used, for example "num_days >= 180 & ER_K <= 0.5".
#' @param metab_rds The standardized metabolism file
#' @param has_NLDAS Optional parameter only for aquatic sites to check if NLDAS data is present
#' @param has_LAI Optional parameter only for aquatic sites to check if LAI data is present
#' @param has_StreamLight Optional parameter only for aquatic sites to check if StreamLight predictions are present
#' @param has_MODIS_NPP Optional parameter only for aquatic sites to check if annual MODIS NPP data is present
#' @return Returns a filtered set of metabolism data
#'
#' @export
#'
#===============================================================================
#Defining a function to filter the dataset
#===============================================================================
filter_metab <- function(diag, filters, metab_rds, has_NLDAS = NULL, has_LAI = NULL, 
  has_StreamLight = NULL, has_MODIS_NPP = NULL){
  #Catch if no filters are applied
    if(is.null(filters)){}
  
  #List of optional parameters
    optional_params <- list(has_NLDAS, has_LAI, has_StreamLight, has_MODIS_NPP)
      names(optional_params) <- c("NLDAS_data", "LAI_data", "Predicted_light", "MODIS_NPP")
      
  #Check if/which optional parameters were used
    optional_use <- do.call(rbind, lapply(optional_params, FUN = is.null))
  
  #Set filters based on which optional parameters are not null
    if(sum(optional_use) != length(optional_use)){
      message("Checking availability of ancillary data for aquatic sites")
      
      final_filters <- paste0(
        filters, 
        "&",
        paste0(rownames(optional_use)[which(optional_use == FALSE)] , " == TRUE", collapse = "&")  
      )
  
    } else {
      final_filters <- filters
    } #End if check for not null optional parameters  

  #Filter the site-years (Update this so it is not deprecated)
    sy_filter <- dplyr::filter_(diag, final_filters)

  #Split the filtered sites by Site_ID
    site_split <- split(sy_filter, sy_filter[, "Site_ID"])
    
  #Helper function to retrieve the standardized data that meets filtering requirements
    retrieve_metab <- function(Site_ID){
      #Years that meet the filtering requirement
        years <- site_split[[Site_ID]][, "Year"]
        
      #Get the standardized metabolism for the year of interest
        metab_SOI <- metab_rds[[Site_ID]]
        
      return(metab_SOI[metab_SOI[, "Year"] %in% years, ])
      
    } #End retrieve_metab function
   
  #Get the filtered metabolism data 
    filtered <- lapply(names(site_split), retrieve_metab)
      names(filtered) <- names(site_split)
      
    message(paste0("A total of ", length(filtered), " sites and ", nrow(sy_filter),
      " site-years meet the filtering criteria"))  
     
  return(filtered)

} #End filter_metab function

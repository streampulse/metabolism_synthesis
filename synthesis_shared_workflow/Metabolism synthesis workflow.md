---
title: "Metabolism synthesis workflow"
author: "Phil Savoy"
date: "10/16/2019"
output: html_document
---

# Description

Description of the workflow for the metabolism synthesis.

Last updated: 9/14/2020

***

Key:

* <i>relative paths</i>
* <span style="color:orange">**R functions**</span>
* <span style="color:SlateBlue">**R executable scripts**</span>
* **Package**

***

# Setup
Before proceeding, load the necessary packages.

```{r Setup envrionment, message = FALSE, results = "hide", eval = FALSE}
#Load the packages  
  package_list <- as.list(c("StreamLightUtils", "StreamLight", "plyr", "stringr", "dplyr", "mgcv", "here"))
  lapply(package_list, library, character.only = TRUE)
```

And create a here file in the project location (first time only)

```{r Create here file, message = FALSE, results= "hide", eval = FALSE}
#Create a here file in the project location (first time only)
  if(file.exists(".here") == FALSE){here::set_here()}
```

# {.tabset .tabset-pills}

---

## I. Acquiring and harmonizing data

---

ADD SOME DESCRIPTIVE TEXT HERE

### 1. Unpacking StreamPULSE data and generating basic site information for the CONUS

Download the following files:

* Appling et al. (2018) metabolism data from [see paper reference](https://www.nature.com/articles/sdata2018292). The full dataset is located [here](https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982) and just the metabolism estimates are available [here](https://www.sciencebase.gov/catalog/item/59eb9c0ae4b0026a55ffe389).
* Download a table of site information titled as "All basic site data" (all_basic_site_data.csv.zip) from https://data.streampulse.org/download_bulk and place in <i>.../data/StreamPULSE</i>
* Bulk download all site data ("All public model objects in RDS format") from https://data.streampulse.org/download_bulk and place in <i>.../data/StreamPULSE</i>


The first step is to prepare the project by getting the necessary metabolism data. 
* A copy of the metabolism estimates from Appling et al. 2018 were placed into <i>../data/Powell</i>. 

```{r Place Appling et al. 2018, results = "hide", eval = FALSE}
#Move a copy of metabolism estimates from Appling et al. 2018 to my local project
  file.copy("C:/research/Datasets/Appling 2018 metabolism/daily_predictions.tsv",
    here("data", "Powell"))  
```

Unzip the downloaded StreamPULSE site information. For ease and transparency, the acquisition date of the data pull is recorded in the workflow.

```{r Site info unzip, results = "hide", eval = FALSE}
#Date of StreamPULSE data acquisition
  acquire_date <- "8_13_2020"

#Unzip the StreamPULSE site information
  unzip(
    here("data", "StreamPULSE", paste0("all_basic_site_data.csv_", acquire_date, ".zip")),
    exdir = here("data", "StreamPULSE")
  )
```

Unzip the downloaded StreamPULSE metabolism estimates (all_sp_model_objects.zip) 

```{r StreamPULSE unzip, results = "hide", eval = FALSE}
#Unzip the StreamPULSE metabolism data
  unzip(
    here("data", "StreamPULSE", paste0("all_sp_model_objects_", acquire_date, ".zip")),
    exdir = here("data", "StreamPULSE")
  )
```


### 2. Generate site information for the CONUS

Basic site information was generated for all sites with modeled metabolism data in the continental United States using <span style="color:SlateBlue">**Lotic_site_basic.R**</span>. A few steps are done at this stage including:

Export a set of basic site information for sites in the CONUS

```{r Basic lotic site information, results = "hide", eval = FALSE}
#Get a formatted table of basic site information for all lotic sites
  source(here("R", "executables", "Lotic_site_basic.R")) 
```

**Output** <i>../output/lotic_site_basic.rds</i>

Calculate the start and end year of metabolism data for each site, as well as individual site-years for each site using <span style="color:SlateBlue">**Lotic_date_ranges.R**</span>.

```{r Calculate date ranges, results = "hide", eval = FALSE}
#Find the date ranges for each site
  source(here("R", "executables", "Lotic_date_ranges.R"))
```

* **Output** <i>../output/lotic_timeframes.rds</i>
* **Output** <i>../output/sp_site_years.rds</i>
* **Output** <i>../output/powell_site_years.rds</i>

Sites were then filtered for the CONUS so they could be associated with various standardized products. Where possible, each site was associated with a NHDPlus_v2 COMID through the use of the **nhdplusTools** package and the simple wrapper functions <span style="color:orange">**comid_from_point**</span>. Similarly, the vector processing unit was associated with each site using the simple wrapper function <span style="color:orange">**vpu_from_point.R**</span>.

Once sites were assigned a COMID, they could be associated with other datasets. Stream order information for each reach was pulled from the NHDplus_v2 flowlines. Additionally, a host of StreamCat metrics were then associated with each reach. All of these operations were done with <span style="color:SlateBlue">**Lotic_site_info.R**</span>.

```{r Advanced site information, results = "hide", eval = FALSE}
source(here("R", "executables", "Lotic_site_info.R"))
```

**Output** <i>../output/lotic_site_info.rds</i>

The final dataset contains the following columns:

* **Site_ID:** Unique site identifier
* **Name:** Site long name
* **Source:** Site data source
* **Lat:** Site Latitude
* **Lon:** Site Longitude
* **epsg_crs:** Site coordinate reference system as EPSG code
* **COMID:** NHDPlus_v2 reach COMID
* **VPU:** Vector processing unit
* **StreamOrde:** Stream order (from NHDPlus_v2 flowlines)
* **PctImp2011Cat:** Mean imperviousness of anthropogenic surfaces (NLCD 2011) within catchment
* **PctImp2011Ws:** Mean imperviousness of anthropogenic surfaces (NLCD 2011) within watershed
* **PctImp2011CatRp100:** Mean imperviousness of anthropogenic surfaces (NLCD 2011) within catchment and within a 100-m buffer of NHD stream lines
* **PctImp2011WsRp100:** Mean imperviousness of anthropogenic surfaces (NLCD 2011) within watershed and within a 100-m buffer of NHD stream lines
* **PermCat:** Mean permeability (cm/hour) of soils (STATSGO) within catchment
* **PermWs:** Mean permeability (cm/hour) of soils (STATSGO) within watershed
* **DamDensCat:** Density of georeferenced dams within catchment (dams/ square km) based on the National Inventory of Dams (https://catalog.data.gov/dataset/national-inventory-of-dams)
* **DamDensWs:** Density of georeferenced dams within watershed (dams/ square km) based on the National Inventory of Dams (https://catalog.data.gov/dataset/national-inventory-of-dams)
* **RdDensCat:** Density of roads (2010 Census Tiger Lines) within catchment (km/square km)
* **RdDensWs:** Density of roads (2010 Census Tiger Lines) within watershed (km/square km)
* **RdDensCatRp100:** Density of roads (2010 Census Tiger Lines) within catchment and within a 100-m buffer of NHD stream lines (km/square km)
* **RdDensWsRp100:** Density of roads (2010 Census Tiger Lines) within watershed and within a 100-m buffer of NHD stream lines (km/square km)
* **RunoffCat:** Mean runoff (mm) within catchment
* **RunoffWs:** Mean runoff (mm) within watershed
* **BFICat:** Baseflow is the component of streamflow that can be attributed to ground-water discharge into streams. The Baseflow Index (BFI) is the ratio of baseflow to total flow, expressed as a percentage, within catchment.
* **BFIWs:** Baseflow is the component of streamflow that can be attributed to ground-water discharge into streams. The Baseflow Index (BFI) is the ratio of baseflow to total flow, expressed as a percentage, within watershed.
* **PctDecid2011Cat:** % of catchment area classified as deciduous forest land cover (NLCD 2011 class 41)
* **PctDecid2011Ws:** % of watershed area classified as deciduous forest land cover (NLCD 2011 class 41)
* **PctDecid2011CatRp100:** % of catchment area classified as deciduous forest land cover (NLCD 2011 class 41) within a 100-m buffer of NHD streams
* **PctDecid2011WsRp100:** % of watershed area classified as deciduous forest land cover (NLCD 2011 class 41) within a 100-m buffer of NHD streams
* **PctConif2011Cat:** % of catchment area classified as evergreen forest land cover (NLCD 2011 class 42)
* **PctConif2011Ws:** % of watershed area classified as evergreen forest land cover (NLCD 2011 class 42)
* **PctConif2011CatRp100:** % of catchment area classified as evergreen forest land cover (NLCD 2011 class 42) within a 100-m buffer of NHD streams
* **PctConif2011WsRp100:** % of watershed area classified as evergreen forest land cover (NLCD 2011 class 42) within a 100-m buffer of NHD streams
* **WetIndexCat:** Mean Composite Topographic Index (CTI) [Wetness Index] within catchment
* **WetIndexWs:** Mean Composite Topographic Index (CTI) [Wetness Index] within watershed
* **prG_BMMI:** Predicted probability that a stream segment is in good biologial condition based on a random forest model of the NRSA benthic invertebrate multimetric index (BMMI)
* **CCONN:** Hydrologic connectivity component score calculated using catchment metrics
* **CHABT:** Habitat provision component score calculated using catchment metrics
* **ICI:** Index of catchment integrity
* **IWI:** Index of watershed integrity

### 3. Downloading and processing remote products (NLDAS and MODIS)    

### 3a. NLDAS data
North American Data Assimilation Systems (NLDAS) downwelling shortwave radiation (W m^-2^) data was downloaded via the data rods service (https://disc.gsfc.nasa.gov/information/tools?title=Hydrology%20Data%20Rods). The download and processing of NLDAS data was handled by the <span style="color:orange">**NLDAS_DL_bulk**</span> and <span style="color:orange">**NLDAS_proc**</span> functions from **StreamLightUtils**

```{r Downloading and processing NLDAS data, results = "hide", eval = FALSE}
source(here("R", "executables", "NLDAS_data_prep.R"))
```

**Output** <i>../output/NLDAS_processed.rds</i>

#### 3b. MODIS data

The MODIS leaf area index (LAI) (m^2^ m^-2^) product (MCD15A2H-006) (500m, 8-day) was acquired through [AppEEARS](https://lpdaacsvc.cr.usgs.gov/appeears/) by generating a set of sites to request data for from 01-01-2005 to 09-01-2020. 

```{r AppEEARS request sites, results = "hide", eval = FALSE}
#Get AppEEARS request sites
  source(here("R", "executables", "appeears_request_sites.R"))
```

After submitting and downloading the AppEEARS request the zipped files were unpacked and processed using the <span style="color:orange">**AppEEARS_unpack**</span> and <span style="color:orange">**AppEEARS_proc**</span> functions from **StreamLightUtils**. 

Not all sites had LAI data in the AppEEARS request, and for those sites data was then downloading using the **MODISTools** package for a 1km buffer around the site location. The timeseries for each pixel was processed separately and then a representative median timeseries of filtered LAI was returned for the site.

Timeseries of LAI were smoothed and gap-filled to extract the seasonal phenological signal following the approach of [Gu et al. (2009)](https://link.springer.com/chapter/10.1007/978-1-4419-0026-5_2) using the **phenofit** R package [found here](https://github.com/kongdd/phenofit) (Kong 2020)

```{r Unpacking and processing LAI data, results = "hide", eval = FALSE}
#After getting the AppEEARS request, download data for missing sites
  source(here("R", "executables", "LAI_get_missing.R"))
        
#Unpack & process AppEEARS data, combine with single pixel data 
  source(here("R", "executables", "LAI_data_prep.R"))
  LAI_data_prep(fit_method = "Gu")
```

**Output** <i>../output/MOD_LAI_processed.rds</i>

### 4. Estimates of light at the stream surface

Estimates of light at the stream surface were generated using the **StreamLight** package. Running the model requires the creating of a parameter file for each site as well as driver files. 

#### 4a. Generate a parameter file for each site.

```{r Generate parameter file, results = "hide", eval = FALSE}
#Generate light model parameter files
  source(here("R", "executables", "StreamLight_params.R"))
```

**Output** <i>../output/lotic_streamlight_params.rds</i>

Some specific details about parameters include:

* **Stream azimuth:** I calculated the circular mean of azimuths for each stream reach based on latitude and longitude of the site location using NHDPlus_v2 hydrography data.

* **Bankfull width** Stream widths were estimated based on data from McManamay & DeRolph (2019) (ADD REFERENCE)

* **Tree height** LiDAR-derived estimates of tree height (m) from Simard et al. (2011) were collected at each site using the <span style="color:orange">**extract_height**</span> function in the **StreamLightUtils** package. Locations that returned a height of 0 were replaced with a weighted mean of tree heights based on landcover information. ADD MORE DETAIL HERE, BUT THE APPROACH IS THE SAME AS DESCRIBED IN MY CURRENT FRESHWATER SCIENCE SUBMISSION.

#### 4b. Generate driver files

Driver files for **StreamLight** were created using NLDAS downwelling shortwave radiation and MODIS LAI data processed earlier. The driver files were made using the <span style="color:orange">**make_driver**</span> function in the **StreamLightUtils** package.

```{r Generate driver files, results = "hide", eval = FALSE}
#Generate driver files
  source(here("R", "executables", "StreamLight_drivers.R"))
```

**Output** Inidivual driver files were output in the form "Site_ID"_driver.rds  in <i>../output/driver files</i>

#### 4c. Run StreamLight to estimate light at the stream surface

The **StreamLight** package was used to estimate light reaching the stream surface at hourly timesteps.

```{r Run StreamLight model, results = "hide", eval = FALSE}
#Estimate light reaching the stream surface
  source(here("R", "executables", "StreamLight_estimates.R"))
```

**Output** <i>../output/modeled_light.rds</i>

Light estimates at the stream surface were then aggregated to daily sums (mol m^-2^ d^-1^) so they could be incorporated into the metabolism synthesis data.

```{r Daily light estimates, results = "hide", eval = FALSE}
#Generate daily results of modeled light
  source(here("R", "executables", "StreamLight_daily.R"))
```

**Output** <i>../output/modeled_light.rds</i>

### 5. Compiling and standardizing the datasets

Compile all of the StreamPULSE data so that the same workflow can be used on both the Powell center and StreamPULSE datasets.

```{r Compile StreamPULSE data, results = "hide", eval = FALSE}
#Compile StreamPULSE data
  source(here("R", "executables", "Compile_StreamPULSE.R"))
```

**Output** <i>../output/sp_daily_predictions.rds</i>

Stream metabolism, riparian LAI and GPP, incoming light, and modeled estimates of light at the stream surface were compiled into a single set of standardized files using <span style="color:SlateBlue">**Lotic_standardized_metabolism.R**</span>, which leverages the <span style="color:orange">**synthesis_format**</span> function. In addition to stream data, incoming PAR, predicted PAR, and processed LAI generated in earlier steps are also included. To ensure that all files are formatted in the same manner, if any datasets are not present then the corresponding columns were populated with NA values.

```{r Standardized files, results = "hide", eval = FALSE}
#Generate standardized files
  source(here("R", "executables", "Lotic_standardized_metabolism.R"))
```

**Output** <i>../output/lotic_standardized_metabolism.rds</i>

The final dataset contains the following columns:

* **Date:** Date in YYYY-MM-DD format
* **U_ID:** Unique date identifier (format as year + DOY)
* **Year:** Year
* **DOY:** Day of year (1:365 or 366)
* **GPP:** Stream GPP estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018)
* **ER:** Stream ER estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018)
* **K600:** Model estimate of K600, the mean reaeration rate coefficient, scaled to a Schmidt number of 600, for this date. Value is the median of the post-warmup MCMC distribution
* **DO.obs:** Mean dissolved oxygen concentration (mg O2 L^-1^) for the date (4am to 359am)
* **DO.sat:** Mean theoretical saturation concentration (mg O2 L^-1^) for the date (4am to 359am)
* **temp.water:** Mean water temperature (degreees C) for the date (4am to 359pm)
* **discharge:** Mean discharge (m^3^ s^-1^) for the date (4am to 359pm)
* **PAR_sum:** Daily sum of incoming above canopy PAR (mol m^-2^ d^-1^)
* **Stream_PAR_sum:** Daily sum of PAR estimated at the stream surface (mol m^-2^ d^-1^)
* **LAI_proc:** MODIS LAI data that has been processed and gap-filled (m^2^ m^-2^)

### 6. Processing terrestrial data

#### 6a. Compile daily FLUXNET data

FLUXNET 2015 daily data was compiled into a similar standardized format as the lotic sites using <span style="color:SlateBlue">**Lotic_standardized_metabolism.R**</span>. Important to note though is that some associated environmental variables are reported as daily means rather than sums (for example precip or SW), please see the individual column descriptions below. **I WOULD NOT USE THE ASSOCIATED ENVIRONMENTAL DATA I INCLUDED FOR FLUXNET DATA BECAUSE OF THESE DIFFERENCES.**

```{r Standardized FLUXNET daily files, results = "hide", eval = FALSE}
#Compile all of the daily FLUXNET data
  source(here("R", "executables", "FLUXNET_daily_compile.R"))
```

**Output** <i>../output/fluxnet_standardized.rds</i>

The final dataset contains the following columns:

* **Date:** Date in YYYY-MM-DD format
* **U_ID:** Unique date identifier (format as year + DOY)
* **Year:** Year
* **DOY:** Day of year (1:365 or 366)
* **GPP:** FLUXNET GPP estimates (g C m^-2^ d^-1^) "GPP_NT_VUT_REF"
* **ER:** FLUXNET ER estimates (g C m^-2^ d^-1^) "RECO_NT_VUT_REF
* **Net:** FLUXNET NEE (changed to NEP by * -1) "NEE_VUT_REF"
* **Temp:** Average air temperature from half-hourly data (degrees C)
* **Precip:** Average precipition from half-hourly data (mm)
* **VPD:** Average vapor pressure deficit from half-hourly data (hPa)
* **SW:** Average incoming shortwave radiation from half-hourly data (W m^-2)

#### 6b. Compile annual FLUXNET data

Rather than making my own calculations, I also compiled the annual productivity data that is already compiled and provided by FLUXNET using <span style="color:SlateBlue">**FLUXNET_annual_compile.R**</span>.

```{r Compile annual FLUXNET, results = "hide", eval = FALSE}
#Compile all of the annnual FLUXNET data*    
  source(here("R", "executables", "FLUXNET_annual_compile.R"))
```

**Output** <i>../output/fluxnet_standardized.rds</i>

The final dataset contains the following columns:

* **Year:** Year
* **GPP:** FLUXNET annual GPP (sum from daily estimates) (g C m^-2^ y^-1^) "GPP_NT_VUT_REF". Gross Primary Production, from Nighttime partitioning method, reference version selected from GPP versions using a model efficiency approach.
* **ER:** FLUXNET annual ER (sum from daily estimates) (g C m^-2^ d^-y^) "RECO_NT_VUT_REF". Ecosystem Respiration, from Nighttime partitioning method, reference selected from RECO versions using a model efficiency approach.
* **Net:** FLUXNET NEE (changed to NEP by * -1) "NEE_VUT_REF". I derived this from data of this description:Net Ecosystem Exchange, using Variable Ustar Threshold (VUT) for each year, reference selected on the basis of the model efficiency.
* **Temp:** Average air temperature from daily data (degrees C)
* **Precip:** Average precipition from daily data (mm)
* **VPD:** Average vapor pressure deficit from daily data (hPa)
* **SW:** Average incoming shortwave radiation from daily data (W m^-2)

#### 6c. Get annual MODIS NPP data

Annual MODIS NPP data (g C m^-22^ y^-1^) (MOD17A3HGF v006) was collected for each site-year of data through a combination of an AppEEARS request and downloading of single pixels similar to the LAI data.

```{r MODIS NPP data, results = "hide", eval = FALSE}
#Unpack MODIS NPP (MOD17A3HGF v006) data from AppEEARS, download data for
#sites without AppEEARS data, and combine into annual NPP for each site-year
  source(here("R", "executables", "NPP_data_prep.R"))
```

**Output** <i>../output/MODIS_annual_NPP.rds</i>

### 7. Calculate diagnostics for each site (both lotic and terrestrial)

#### 7a. Lotic diagnostics

For each site-year the following set of diagnostics were calculated with the help of the <span style="color:orange">**diagnostics_fun**</span> located in this project folder (<i>../R/functions)</i>.

```{r Metabolism diagnostics, results = "hide", eval = FALSE}
#Calculate diagnostics for each site-year of aquatic data as well as ancillary
#data availability*
  source(here("R", "executables", "Lotic_diagnostics.R")) 
```

**Output** <i>../output/yearly_diagnostics.rds</i>

The output contains the following the following columns:

* **Site_ID:** Unique site identifier
* **Year:** Year
* **ER_K** Correlation (R^2^) between ER and K600
* **GPP_neg** Percentage of days with negative stream GPP estimates
* **ER_pos** Percentage of days with positive stream ER estimates
* **num_days** The number of days with estimated stream metabolism (excluding days with negative GPP and positive ER)
* **NLDAS_data:** Logical indicator (TRUE/FALSE) for availability of NLDAS data for this site-year
* **LAI_data:** Logical indicator (TRUE/FALSE) for availability of LAI data for this site-year
* **Predicted_light:** Logical indicator (TRUE/FALSE) for availability of StreamLight predictions of light at the stream surface for this site-year
* **MODIS_NPP:** Logical indicator (TRUE/FALSE) for availability of MODIS annual NPP data for this site-year

---

## II. Filtering, gap-filling, and calculating metrics

---

### 1. Filtering and gap-filling the datasets   

#### 1a. Filter the lotic metabolism dataset

At this point all of the data for any site with metabolism data in the CONUS has been compiled. Analysis for the metabolism synthesis paper will be derived from a filtered subset of this dataset. At present the only filter applied is that site-years must have stream metabolism estimates (excluding days with negative GPP estimates or positive ER estimates) for at least 60% of the year. For the lotic metabolism data I also assumed that all NLDAS, LAI, light predictions, and MODIS NPP data should be available, but this can easily be changed by removing these arguments from <span style="color:orange">**filter_metab**</span>.

```{r Filter metabolism data, results = "hide", eval = FALSE}
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
```

**Output** <i>../output/synthesis_filtered.rds</i>

#### 1b. Gap-fill lotic sites sites

Each site-year of data was fed through several functions that handled gap-filling and z-normalization of data (ADD MORE DETAILS ABOUT THESE ROUTINES AT A LATER DATE). 

At this stage several notable steps occurred:

* All negative stream GPP estimates were replaced with NA
* All positve stream ER estimates were replaced with NA
* Stream NEP (GPP + ER) was calculated
* Additional columns for stream GPP, ER, and NEP expressed in units of (g C m^-2^ d^-1) were calculated based on a photosynthetic quotient (PQ) of 1.25.
* Gap-filling of stream metabolism, water temperature, daily incoming PAR, discharge
* Z-normalization of stream metabolism, water temperature, daily incoming PAR

```{r Gap-fill data, results = "hide", eval = FALSE}
#Read in a suite of gap-filling functions
  source(here("R", "functions", "gapfill_toolbox.R"))  
  
  gap_filled <- synthesis_gapfill(
    synthesis_filtered = readRDS(here("output", "synthesis_filtered.rds")), 
    PQ = 1.25
  )    
      
#Save the output
  saveRDS(gap_filled, here("output", "synthesis_gap_filled.rds"))
```

**Output** <i>../output/synthesis_gap_filled.rds</i>

#### 1c. Filter the FLUXNET dataset 

Note, the daily FLUXNET data does not need to be gap-filled to calculate annual sums, since that is taken care of already in the yearly FLUXNET products.

```{r Filter FLUXNET, results = "hide", eval = FALSE}
#Filter the standardized dataset
  FLUXNET_filtered <- filter_metab(
    diag = readRDS(here("output", "FLUXNET_yearly_diagnostics.rds")), 
    filters = paste0("num_days >= ", 365 * 0.6), 
    metab_rds = readRDS(here("output", "FLUXNET_standardized.rds"))
  )  
      
#Save the output
  saveRDS(FLUXNET_filtered, here("output", "FLUXNET_filtered.rds"))
```

### 2. Calculate a series of data metrics for each site 

#### 2a. Metrics for the lotic synthesis sites

```{r Synthesis site metrics, results = "hide", eval = FALSE}
#Calculate site metrics for synthesis sites
  source(here("R", "executables", "Synthesis_site_metrics.R"))
```

**Output** <i>../output/synthesis_site_metrics.rds</i>

The output contains the following the following columns:

* **Site_ID:** Unique site identifier
* **ann_GPP_C:** Mean annual cumulative stream GPP (g C m^-2^ y^-1^). This was calculated by first calculating annual sums of GPP (g C m^-2^ y^-1^) for each site year, and then taking the mean annual rate for each site.
* **upper_GPP_C:** 95th percentile of daily rates of stream GPP (g C m^-2^ d^-1^).
* **ann_ER_C:** Mean annual cumulative stream ER (g C m^-2^ y^-1^). This was calculated by first calculating annual sums of ER (g C m^-2^ y^-1^) for each site year, and then taking the mean annual rate for each site.
* **lower_ER_C:** 5th percentile of daily rates of stream GPP (g C m^-2^ d^-1^). Since ER is negative, you can think of this as equivalent to the 95th percentile done for GPP. 
* **PAR_sum:**  Mean annual cumulative incoming PAR (kmol m^-2^ y^-1^). This was calculated by first calculating annual sums of PAR (kmol m^-2^ y^-1^) for each site year, and then taking the mean annual rate for each site.
* **Stream_PAR_sum:**  Mean annual cumulative predicted PAR at the stream surface (kmol m^-2^ y^-1^). This was calculated by first calculating annual sums of predicted PAR at the stream surface (kmol m^-2^ y^-1^) for each site year, and then taking the mean annual rate for each site.
* **gpp_C_mean:** Mean daily GPP (g C m^-2^ d^-1^)
* **gpp_C_cv:** CV of daily GPP
* **gpp_C_skew:** Skewness of daily GPP
* **gpp_C_kurt:** Kurtosis of daily GPP
* **gpp_C_amp:** Amplitude of daily GPP
* **gpp_C_phase:** Phase of daily GPP (day of year)
* **gpp_C_ar1:** Autoregressive lag-one correlation coefficient of daily GPP
* **er_C_mean:** Mean daily ER (g C m^-2^ d^-1^)
* **er_C_cv:** CV of daily ER
* **er_C_skew:** Skewness of daily ER
* **er_C_kurt:** Kurtosis of daily ER
* **er_C_amp:** Amplitude of daily ER
* **er_C_phase:** Phase of daily ER (day of year)
* **er_C_ar1:** Autoregressive lag-one correlation coefficient of daily ER
* **Wtemp_mean:** Mean daily water temperature
* **Wtemp_cv:** CV of daily GPP
* **Wtemp_skew:** Skewness of daily water temperature
* **Wtemp_kurt:** Kurtosis of daily water temperature
* **Wtemp_amp:** Amplitude of daily water temperature
* **Wtemp_phase:** Phase of daily water temperature (day of year)
* **Wtemp_ar1:** Autoregressive lag-one correlation coefficient of daily water temperature
* **Disch_mean:** Mean daily discharge
* **Disch_cv:** CV of daily discharge
* **Disch_skew:** Skewness of daily discharge
* **Disch_kurt:** Kurtosis of daily discharge
* **Disch_amp:** Amplitude of daily discharge
* **Disch_phase:** Phaste of daily discharge (day of year)
* **Disch_ar1:** Autoregressive lag-one correlation coefficient of daily discharge
* **PAR_mean:** Mean daily PAR
* **PAR_cv:** CV of daily PAR
* **PAR_skew:** Skewness of daily PAR
* **PAR_kurt:** Kurtosis of daily PAR
* **PAR_amp:** Amplitude of daily PAR
* **PAR_phase:** Phase of daily PAR (day of year)
* **PAR_ar1:** Autoregressive lag-one correlation coefficient of daily PAR
* **LAI_mean:** Mean daily LAI (m^2^ m^-2^)
* **LAI_cv:** CV of LAI
* **LAI_skew:** Skewness of LAI
* **LAI_kurt:** Kurtosis of LAI
* **LAI_amp:** Amplitude of LAI
* **LAI_phase:** Phase of LAI  (day of year)
* **LAI_ar1:** Autoregressive lag-one correlation coefficient of daily LAI
* **MOD_ann_NPP:** Mean annual MODIS NPP (g C m^-2^ y^-1^) for the concurrent period of record for stream metabolism data at a site. Annual sums of NPP (g C m^-2^ d^-y^) were available for each site year and then the mean was taken to get a mean annual rate for each site.

#### 2b. Metrics for the FLUXNET sites

```{r FLUXNET site metrics, results = "hide", eval = FALSE}
#Calculate site metrics for FLUXNET sites
  source(here("R", "executables", "FLUXNET_site_metrics.R"))
```

**Output** <i>../output/synthesis_site_metrics.rds</i>

The output contains the following the following columns:

* **Site_ID:** Unique site identifier
* **ann_GPP_C:** Mean annual cumulative GPP (g C m^-2^ y^-1^). This was calculated from the annual sums of GPP (g C m^-2^ y^-1^) for each site year provided by FLUXNET, and then taking the mean annual rate for each site.
* **upper_GPP_C:** 95th percentile of daily rates of GPP (g C m^-2^ d^-1^).
* **ann_ER_C:** Mean annual cumulative ER (g C m^-2^ y^-1^). This was calculated from the annual sums of ER (g C m^-2^ y^-1^) for each site year provided by FLUXNET, and then taking the mean annual rate for each site
* **lower_ER_C:** 5th percentile of daily rates of GPP (g C m^-2^ d^-1^). Since ER is negative, you can think of this as equivalent to the 95th percentile done for GPP.  
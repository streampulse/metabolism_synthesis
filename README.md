<font size="6">Summary of files in this repository </font> 

This repository contains data and code used in [Bernhardt et al. (2022)](https://doi.org/10.1073/pnas.2121976119). The preparation and analysis of this dataset consisted of a four part workflow:

* **Part I: Standardized dataset creation:** Part I is the preparation of standardized datasets of stream metabolism from the [StreamPULSE data portal](https://data.streampulse.org/) and terrestrial ecosystem fluxes from the [FLUXNET2015 dataset](https://fluxnet.org/data/fluxnet2015-dataset/). Everything generated in Part I is included in `output_data` so that other users could use this larger dataset if they choose to work with a different subset of data than we used.
* **Part II. Filtering, gap-filling, and calculating metrics:** In Part II the standardized dataset of lotic and terrestrial metabolism data was filtered down to the subset of sites used in Bernhardt et al. (2022). After filtering, additional descriptive metrics were calculated to use in analysis.
* **Part III. Plotting and export of stats dataset:** Datasets from part II are minimally subsetted and recast for plotting convenience. Figures 1-6 are generated. An analysis-ready dataset is compiled.
* **Part IV. Structural equation modeling:** Code used for structural equation model (SEM) analysis on annual metabolism dataset. We used an observed variables model to estimates the effect of light (PAR reaching the stream surface) and hydrologic variability (skewness of daily discharge) on annual river GPP. bles model to estimates the effect of light (PAR reaching the stream surface) and hydrologic variability (skewness of daily discharge) on annual river GPP.  


<font size="4">The repository contains: </font> 
-  A documentation of each parts of the workflow, located in the documentation folder of this repository. The original R markdown files are also included in the root directory of the repository.
- An associated package called **BernhardtMetabolism**, which contains data and code used to produce Part II of the workflow
- Final data outputs

<font size="4">Key: </font> 
<font size="3"><i><span style="color:#00cc99">Exported datasets</span></i></font>
**<span style="color:#009688">"Column name"</span>**

## 1. Workflow documentation

- **Part-I-Standardized-dataset-creation.html :** Documentation of Part I of the workflow. This document only contains a description of the data, methods, and and steps used to complete Part I of the workflow.
- **Part-II-Filtering-and-gap-filling.html:** Documentation of Part II of the workflow. This docoument, along with the associated R markdown file (Part II-Filtering and gap filling.Rmd) goes through the entire set of steps used to complete Part II of the workflow. In conjunction with the **BernhardtMetabolism** package, this document contains all of the data and code to fully reproduce Part II of the workflow.

## 2. Exported datasets

A final set of outputs were then exported to the output_data directory of this repository.

### 2.1 Standardized, unfiltered, lotic metabolism

This is the full, unfiltered, standardized dataset of stream metabolism that was compiled.

<font size="3"><i><span style="color:#00cc99">lotic_site_info_full.rds</span></i></font>

Format: A single data frame with the following columns:

* **<span style="color:#009688">"Site_ID":</span>** Unique site identifier
* **<span style="color:#009688">"Name":</span>** Site long name
* **<span style="color:#009688">"Source":</span>** Site data source
* **<span style="color:#009688">"Lat":</span>** Site Latitude
* **<span style="color:#009688">"Lon":</span>** Site Longitude
* **<span style="color:#009688">"epsg_crs":</span>** Site coordinate reference system as EPSG code
* **<span style="color:#009688">"COMID":</span>** NHDPlus_v2 reach COMID
* **<span style="color:#009688">"VPU":</span>** Vector processing unit
* **<span style="color:#009688">"StreamOrde":</span>** Stream order (from NHDPlus_v2 flowlines)
* **<span style="color:#009688">"Azimuth":</span>** Channel azimuth calculated as the circular mean of azimuths for each stream reach based on latitude and longitude of the site location using NHDPlus_v2 hydrography data.
* **<span style="color:#009688">"TH":</span>** Tree height (m). Estimates were derived from Tree heights were derived using 30m resolution global canopy height estimates from [Potapov et al. (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0034425720305381?via%3Dihub)
* **<span style="color:#009688">"Width":</span>** Channel width (m)
* **<span style="color:#009688">"Width_src":</span>** Source of channel width estimates. Values include: (NWIS field measurements, Regional geomorphic scaling coeff, StreamPULSE estimates)
* **<span style="color:#009688">"WS_area_km2":</span>** Watershed area (km ^-2^)
* **<span style="color:#009688">"WS_area_src":</span>** Source of watershed area estimates. Values include: (Appling2018_USGS2013, Appling2018_StreamStats, nwis_site_description, StreamStats, localuser_HBFLTER, localuser_UNHWQAL)

<font size="3"><i><span style="color:#00cc99">lotic_standardized_full.rds</span></i></font>

Format: A list of data frames, where each element of the list is a data frame of timeseries for a single site. The names of each list element correspond to the unique site identifier (Site_ID) for a site. Each data frame contains the following columns:

* **<span style="color:#009688">"Site_ID":</span>** Unique site identifier
* **<span style="color:#009688">"Date":</span>** Date in YYYY-MM-DD format
* **<span style="color:#009688">"U_ID":</span>** Unique date identifier (format as year + DOY)
* **<span style="color:#009688">"Year":</span>** Year
* **<span style="color:#009688">"DOY":</span>** Day of year (1 to 365 or 366)
* **<span style="color:#009688">"GPP_raw":</span>** Stream GPP estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018), raw data
* **<span style="color:#009688">"ER_raw":</span>** Stream ER estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018), raw data
* **<span style="color:#009688">"GPP":</span>** Stream GPP estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018), negative values replaced with NA
* **<span style="color:#009688">"ER":</span>** Stream ER estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018), positive values replaced with NA
* **<span style="color:#009688">"K600":</span>** Model estimate of K600, the mean reaeration rate coefficient, scaled to a Schmidt number of 600, for this date. Value is the median of the post warmup MCMC distribution
* **<span style="color:#009688">"DO.obs":</span>** Mean dissolved oxygen concentration (mg O2 L^-1^) for the date (4am to 3:59am)
* **<span style="color:#009688">"DO.sat":</span>** Mean theoretical saturation concentration (mg O2 L^-1^) for the date (4am to 3:59am)
* **<span style="color:#009688">"temp.water":</span>** Mean water temperature (degreees C) for the date (4am to 3:59pm)
* **<span style="color:#009688">"discharge":</span>** Mean discharge (m^3^ s^-1^) for the date (4am to 3:59pm)
* **<span style="color:#009688">"PAR_sum":</span>** Daily sum of incoming above canopy PAR (mol m^-2^ d^-1^)
* **<span style="color:#009688">"Stream_PAR_sum":</span>** Daily sum of PAR estimated at the stream surface (mol m^-2^ d^-1^)
* **<span style="color:#009688">"LAI_proc":</span>** MODIS LAI data that has been processed and gap-filled (m^2^ m^-2^)

### 2.2 Standardized and filtered data used in Bernhardt et al. (2022) 

#### 2.2.1 Standardized, filtered, gap-filled lotic metabolism

This is the dataset of stream metabolism used in Bernhardt et al. (2022). This data has been filtered based on several measures of data quality and availability and gaps were gap-filled. Finally, a suite of site metrics were calculated for use in analysis.

<font size="3"><i><span style="color:#00cc99">lotic_site_info_filtered.rds</span></i></font>

Format: A single data frame with the following columns:

* **<span style="color:#009688">"Site_ID":</span>** Unique site identifier
* **<span style="color:#009688">"Name":</span>** Site long name
* **<span style="color:#009688">"Source":</span>** Site data source
* **<span style="color:#009688">"Lat":</span>** Site Latitude
* **<span style="color:#009688">"Lon":</span>** Site Longitude
* **<span style="color:#009688">"epsg_crs":</span>** Site coordinate reference system as EPSG code
* **<span style="color:#009688">"COMID":</span>** NHDPlus_v2 reach COMID
* **<span style="color:#009688">"VPU":</span>** Vector processing unit
* **<span style="color:#009688">"StreamOrde":</span>** Stream order (from NHDPlus_v2 flowlines)
* **<span style="color:#009688">"Azimuth":</span>** Channel azimuth calculated as the circular mean of azimuths for each stream reach based on latitude and longitude of the site location using NHDPlus_v2 hydrography data.
* **<span style="color:#009688">"TH":</span>** Tree height (m). Estimates were derived from Tree heights were derived using 30m resolution global canopy height estimates from [Potapov et al. (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0034425720305381?via%3Dihub)
* **<span style="color:#009688">"Width":</span>** Channel width (m)
* **<span style="color:#009688">"Width_src":</span>** Source of channel width estimates. Values include: (NWIS field measurements, Regional geomorphic scaling coeff, StreamPULSE estimates)
* **<span style="color:#009688">"WS_area_km2":</span>** Watershed area (km ^-2^)
* **<span style="color:#009688">"WS_area_src":</span>** Source of watershed area estimates. Values include: (Appling2018_USGS2013, Appling2018_StreamStats, nwis_site_description, StreamStats, localuser_HBFLTER, localuser_UNHWQAL)
* **<span style="color:#009688">"ann_GPP_C":</span>** Mean annual cumulative stream GPP (g C m^-2^ y^-1^). This was calculated by first calculating annual sums of GPP (g C m^-2^ y^-1^) for each site year, and then taking the mean annual rate for each site.
* **<span style="color:#009688">"upper_GPP_C":</span>** 95th percentile of daily rates of stream GPP (g C m^-2^ d^-1^).
* **<span style="color:#009688">"ann_ER_C":</span>** Mean annual cumulative stream ER (g C m^-2^ y^-1^). This was calculated by first calculating annual sums of ER (g C m^-2^ y^-1^) for each site year, and then taking the mean annual rate for each site.
* **<span style="color:#009688">"lower_ER_C":</span>** 5th percentile of daily rates of stream ER (g C m^-2^ d^-1^). Since ER is negative, you can think of this as equivalent to the 95th percentile done for GPP. 
* **<span style="color:#009688">"PAR_sum":</span>**  Mean annual cumulative incoming PAR (kmol m^-2^ y^-1^). This was calculated by first calculating annual sums of PAR (kmol m^-2^ y^-1^) for each site year, and then taking the mean annual rate for each site.
* **<span style="color:#009688">"Stream_PAR_sum":</span>**  Mean annual cumulative predicted PAR at the stream surface (kmol m^-2^ y^-1^). This was calculated by first calculating annual sums of predicted PAR at the stream surface (kmol m^-2^ y^-1^) for each site year, and then taking the mean annual rate for each site.
* **<span style="color:#009688">"gpp_C_mean":</span>** Mean daily GPP (g C m^-2^ d^-1^)
* **<span style="color:#009688">"gpp_C_cv":</span>** CV of daily GPP
* **<span style="color:#009688">"gpp_C_skew":</span>** Skewness of daily GPP
* **<span style="color:#009688">"gpp_C_kurt":</span>** Kurtosis of daily GPP
* **<span style="color:#009688">"gpp_C_amp":</span>** Amplitude of daily GPP
* **<span style="color:#009688">"gpp_C_phase":</span>** Phase of daily GPP (day of year)
* **<span style="color:#009688">"gpp_C_ar1":</span>** Autoregressive lag-one correlation coefficient of daily GPP
* **<span style="color:#009688">"er_C_mean":</span>** Mean daily ER (g C m^-2^ d^-1^)
* **<span style="color:#009688">"er_C_cv":</span>** CV of daily ER
* **<span style="color:#009688">"er_C_skew":</span>** Skewness of daily ER
* **<span style="color:#009688">"er_C_kurt":</span>** Kurtosis of daily ER
* **<span style="color:#009688">"er_C_amp":</span>** Amplitude of daily ER
* **<span style="color:#009688">"er_C_phase":</span>** Phase of daily ER (day of year)
* **<span style="color:#009688">"er_C_ar1":</span>** Autoregressive lag-one correlation coefficient of daily ER
* **<span style="color:#009688">"Wtemp_mean** Mean daily water temperature
* **<span style="color:#009688">"Wtemp_cv":</span>** CV of daily GPP
* **<span style="color:#009688">"Wtemp_skew":</span>** Skewness of daily water temperature
* **<span style="color:#009688">"Wtemp_kurt":</span>** Kurtosis of daily water temperature
* **<span style="color:#009688">"Wtemp_amp":</span>** Amplitude of daily water temperature
* **<span style="color:#009688">"Wtemp_phase":</span>** Phase of daily water temperature (day of year)
* **<span style="color:#009688">"Wtemp_ar1":</span>** Autoregressive lag-one correlation coefficient of daily water temperature
* **<span style="color:#009688">"Disch_mean":</span>** Mean daily discharge
* **<span style="color:#009688">"Disch_cv":</span>** CV of daily discharge
* **<span style="color:#009688">"Disch_skew":</span>** Skewness of daily discharge
* **<span style="color:#009688">"Disch_kurt":</span>** Kurtosis of daily discharge
* **<span style="color:#009688">"Disch_amp":</span>** Amplitude of daily discharge
* **<span style="color:#009688">"Disch_phase":</span>** Phaste of daily discharge (day of year)
* **<span style="color:#009688">"Disch_ar1":</span>** Autoregressive lag-one correlation coefficient of daily discharge
* **<span style="color:#009688">"PAR_mean":</span>** Mean daily PAR
* **<span style="color:#009688">"PAR_cv":</span>** CV of daily PAR
* **<span style="color:#009688">"PAR_skew":</span>** Skewness of daily PAR
* **<span style="color:#009688">"PAR_kurt":</span>** Kurtosis of daily PAR
* **<span style="color:#009688">"PAR_amp":</span>** Amplitude of daily PAR
* **<span style="color:#009688">"PAR_phase":</span>** Phase of daily PAR (day of year)
* **<span style="color:#009688">"PAR_ar1":</span>** Autoregressive lag-one correlation coefficient of daily PAR
* **<span style="color:#009688">"LAI_mean":</span>** Mean daily LAI (m^2^ m^-2^)
* **<span style="color:#009688">"LAI_cv":</span>** CV of LAI
* **<span style="color:#009688">"LAI_skew":</span>** Skewness of LAI
* **<span style="color:#009688">"LAI_kurt":</span>** Kurtosis of LAI
* **<span style="color:#009688">"LAI_amp":</span>** Amplitude of LAI
* **<span style="color:#009688">"LAI_phase":</span>** Phase of LAI  (day of year)
* **<span style="color:#009688">"LAI_ar1":</span>** Autoregressive lag-one correlation coefficient of daily LAI
* **<span style="color:#009688">"MOD_ann_NPP":</span>** Mean annual MODIS NPP (g C m^-2^ y^-1^) for the concurrent period of record for stream metabolism data at a site. Annual sums of NPP (g C m^-2^ d^-y^) were available for each site year and then the mean was taken to get a mean annual rate for each site.
* **<span style="color:#009688">"ndays":</span>** Total number of days with daily GPP (non gap-filled) for the site in the filtered dataset
* **<span style="color:#009688">"nyears":</span>** Total number of years for the site in the filtered dataset
* **<span style="color:#009688">"coverage":</span>** Total coverage of daily GPP (non gap-filled) for the site, calculated as ndays / all possible days for all site-years included in the filtered dataset. Ranges from 0-1.

<font size="3"><i><span style="color:#00cc99">lotic_gap_filled.rds</span></i></font>

Format: A list of data frames, where each element of the list is a data frame of timeseries for a single site. The names of each list element correspond to the unique site identifier (Site_ID) for a site. Each data frame contains the following columns:

* **<span style="color:#009688">"Site_ID":</span>** Unique site identifier
* **<span style="color:#009688">"Date":</span>** Date in YYYY-MM-DD format
* **<span style="color:#009688">"U_ID":</span>** Unique date identifier (format as year + DOY)
* **<span style="color:#009688">"Year":</span>** Year
* **<span style="color:#009688">"DOY":</span>** Day of year (1 to 365 or 366)
* **<span style="color:#009688">"GPP_raw":</span>** Stream GPP estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018), raw data
* **<span style="color:#009688">"ER_raw":</span>** Stream ER estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018), raw data
* **<span style="color:#009688">"GPP":</span>** Stream GPP estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018), negative values replaced with NA
* **<span style="color:#009688">"ER":</span>** Stream ER estimates (g O2 m^-2^ d^-1^) from Appling et al. (2018), positive values replaced with NA
* **<span style="color:#009688">"K600":</span>** Model estimate of K600, the mean reaeration rate coefficient, scaled to a Schmidt number of 600, for this date. Value is the median of the post warmup MCMC distribution
* **<span style="color:#009688">"DO.obs":</span>** Mean dissolved oxygen concentration (mg O2 L^-1^) for the date (4am to 3:59am)
* **<span style="color:#009688">"DO.sat":</span>** Mean theoretical saturation concentration (mg O2 L^-1^) for the date (4am to 3:59am)
* **<span style="color:#009688">"temp.water":</span>** Mean water temperature (degreees C) for the date (4am to 3:59pm)
* **<span style="color:#009688">"discharge":</span>** Mean discharge (m^3^ s^-1^) for the date (4am to 3:59pm)
* **<span style="color:#009688">"PAR_sum":</span>** Daily sum of incoming above canopy PAR (mol m^-2^ d^-1^)
* **<span style="color:#009688">"Stream_PAR_sum":</span>** Daily sum of PAR estimated at the stream surface (mol m^-2^ d^-1^)
* **<span style="color:#009688">"LAI_proc":</span>** MODIS LAI data that has been processed and gap-filled (m^2^ m^-2^)
* **<span style="color:#009688">"GPP_filled":</span>** Gap-filled stream GPP estimates (g O2 m^-2^ d^-1^)
* **<span style="color:#009688">"ER_filled":</span>** Gap-filled stream ER estimates (g O2 m^-2^ d^-1^)
* **<span style="color:#009688">"NEP_filled":</span>** Gap-filled stream NEP estimates (g O2 m^-2^ d^-1^)
* **<span style="color:#009688">"GPP_C_filled":</span>** Gap-filled stream GPP estimates, expressed in carbon (g C m^-2^ d^-1^)
* **<span style="color:#009688">"ER_C_filled":</span>** Gap-filled stream ER estimates, expressed in carbon (g C m^-2^ d^-1^)
* **<span style="color:#009688">"NEP_c_filled":</span>** Gap-filled stream NEP estimates, expressed in carbon (g C m^-2^ d^-1^)
* **<span style="color:#009688">"Wtemp_filled":</span>** Gap-filled mean water temperature from the "temp.water" column.
* **<span style="color:#009688">"Disch_filled":</span>** Gap-filled mean discharge from the "discharge" column
* **<span style="color:#009688">"PAR_filled":</span>** Gap-filled daily sum of incoming above canopy PAR, from the "PAR_sum column
* **<span style="color:#009688">"GPP_norm":</span>** Z-normalized GPP
* **<span style="color:#009688">"ER_norm":</span>** Z-normalized ER
* **<span style="color:#009688">"NEP_norm":</span>** Z-normalized NEP
* **<span style="color:#009688">"Wtemp_norm":</span>** Z-normalized water temperature
* **<span style="color:#009688">"PAR_norm":</span>** Z-normalized daily sum of incoming above canopy PAR

#### 2.2.2 Standardized and filtered FLUXNET 

This is the dataset of terrestrial ecosystem fluxes used in Bernhardt et al. (2022). This data has been filtered based on data availability and several basic site metrics were calculated for use in analysis.

<font size="3"><i><span style="color:#00cc99">fluxnet_site_info_filtered.rds</span></i></font>

Format: A single data frame with the following columns:

* **<span style="color:#009688">"Site_ID":</span>** Unique site identifier
* **<span style="color:#009688">"Name":</span>** Site long name
* **<span style="color:#009688">"Lat":</span>** Site Latitude
* **<span style="color:#009688">"Lon":</span>** Site Longitude
* **<span style="color:#009688">"ann_GPP":</span>** Mean annual cumulative GPP (g C m^-2^ y^-1^). This was calculated from the annual sums of GPP (g C m^-2^ y^-1^) for each site year provided by FLUXNET, and then taking the mean annual rate for each site.
* **<span style="color:#009688">"upper_GPP":</span>** 95th percentile of daily rates of GPP (g C m^-2^ d^-1^).
* **<span style="color:#009688">"ann_ER":</span>** Mean annual cumulative ER (g C m^-2^ y^-1^). This was calculated from the annual sums of ER (g C m^-2^ y^-1^) for each site year provided by FLUXNET, and then taking the mean annual rate for each site
* **<span style="color:#009688">"lower_ER":</span>** 5th percentile of daily rates of GPP (g C m^-2^ d^-1^). Since ER is negative, you can think of this as equivalent to the 95th percentile done for GPP. 
* **<span style="color:#009688">"ndays":</span>** Total number of days with daily GPP (non gap-filled) for the site in the filtered dataset
* **<span style="color:#009688">"nyears":</span>** Total number of years for the site in the filtered dataset
* **<span style="color:#009688">"coverage":</span>** Total coverage of daily GPP (non gap-filled) for the site, calculated as ndays / all possible days for all site-years included in the filtered dataset. Ranges from 0-1.

<font size="3"><i><span style="color:#00cc99">fluxnet_filtered_metabolism.rds</span></i></font>

Format: A list of data frames, where each element of the list is a data frame of timeseries for a single site. The names of each list element correspond to the unique site identifier (Site_ID) for a site. Each data frame contains the following columns:

* **<span style="color:#009688">"Date":</span>** Date in YYYY-MM-DD format
* **<span style="color:#009688">"U_ID":</span>** Unique date identifier (format as year + DOY)
* **<span style="color:#009688">"Year":</span>** Year
* **<span style="color:#009688">"DOY":</span>** Day of year (1 to 365 or 366)
* **<span style="color:#009688">"GPP_raw":</span>** FLUXNET annual GPP (sum from daily estimates) (g C m^-2^ y^-1^) "GPP_NT_VUT_REF", raw data. Gross Primary Production, from Nighttime partitioning method, reference version selected from GPP versions using a model efficiency approach.
* **<span style="color:#009688">"ER_raw":</span>** FLUXNET annual ER (sum from daily estimates) (g C m^-2^ d^-y^) "RECO_NT_VUT_REF", raw data. Ecosystem Respiration, from Nighttime partitioning method, reference selected from RECO versions using a model efficiency approach.
* **<span style="color:#009688">"GPP":</span>** FLUXNET annual GPP (sum from daily estimates) (g C m^-2^ y^-1^) "GPP_NT_VUT_REF", negative values replaced with NA. Gross Primary Production, from Nighttime partitioning method, reference version selected from GPP versions using a model efficiency approach. 
* **<span style="color:#009688">"ER":</span>** FLUXNET annual ER (sum from daily estimates) (g C m^-2^ d^-y^) "RECO_NT_VUT_REF", positive values replaced with NA. Ecosystem Respiration, from Nighttime partitioning method, reference selected from RECO versions using a model efficiency approach.
* **<span style="color:#009688">"Net":</span>** FLUXNET NEE (changed to NEP by * -1) "NEE_VUT_REF". I derived this from data of this description:Net Ecosystem Exchange, using Variable Ustar Threshold (VUT) for each year, reference selected on the basis of the model efficiency.
* **<span style="color:#009688">"Temp":</span>** Average air temperature from daily data (degrees C)
* **<span style="color:#009688">"Precip":</span>** Average precipition from daily data (mm)
* **<span style="color:#009688">"VPD":</span>** Average vapor pressure deficit from daily data (hPa)
* **<span style="color:#009688">"SW":</span>** Average incoming shortwave radiation from daily data (W m^-2)
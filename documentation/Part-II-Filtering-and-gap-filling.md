Part II: Filtering, gap-filling, and calculating metrics
================
Phil Savoy, Mike Vlah, Lauren Koenig
Last updated on 18 January, 2022

# Description

<font size="5">Summary of major workflow components</font>

This workflow is part of a series that describes the preparation of data
used in Bernhardt et al. (2022). The creation of this dataset consists
of several sub workflows as described below.

-   **Part I: Standardized dataset creation:** Part I is the preparation
    of standardized datasets of stream metabolism from the [StreamPULSE
    data portal](https://data.streampulse.org/) and terrestrial
    ecosystem fluxes from the [FLUXNET2015
    dataset](https://fluxnet.org/data/fluxnet2015-dataset/). Everything
    generated in Part I is included in `output_data` so that other users
    could use this larger dataset if they choose to work with a
    different subset of data than we used.

-   **Part II. Filtering, gap-filling, and calculating metrics:** In
    Part II the standardized dataset of lotic and terrestrial metabolism
    data was filtered down to the subset of sites used in Bernhardt et
    al. (2022). After filtering, additional descriptive metrics were
    calculated to use in analysis.

-   **Part III. Plotting and export of stats dataset:** Datasets from
    part II are minimally subsetted and recast for plotting convenience.
    Figures 1-6 are generated. An analysis-ready dataset is compiled.

<font size="5">In this document</font>

-   **Filtering and gap-filling data:** Data diagnostics were calculated
    for each site-year of data that described estimates of GPP and ER
    from the lotic and terrestrial sites. Based on these diagnostics,
    the datasets were then filtered to retain sites with sufficient
    information and model performance for later analysis. The filtered
    sites were then gap-filled to facilitate the calculation of
    descriptive metrics in the next step.

-   **Calculating site metrics:** A series of metrics were calculated
    for each site, for both metabolism and environmental data. For lotic
    sites, an additional series of metrics were calculated based on
    adapted ecologically relevant streamflow statistics described in
    [Archfield et
    al. (2013)](https://onlinelibrary.wiley.com/doi/abs/10.1002/rra.2710).

-   **Exporting the final dataset:** Several datasets are exported for
    use in the analysis described in Part III of this workflow
    documentation. This section details the exported datasets and
    provides column definitions for each piece of exported data.

<font size="5">Key</font>

Throughout this document certain terms are indicated with changes in the
color or emphasis of text. Please refer to the following key for their
meaning:

------------------------------------------------------------------------

-   <span style="color:orange">**R functions**</span>
-   **R Package**
-   **<span style="color:#009688">“Column name”</span>**

------------------------------------------------------------------------

This workflow utilizes a package called **BernhardtMetabolism**. This
package bundles together several files created in Part I of the
workflow, as well as all of the necessary code to reproduce the analysis
described in this document (Part II of the workflow). Before proceeding,
make sure that the package is installed.

``` r
install.packages('./BernhardtMetabolism_0.1.1.zip', repos = NULL)

#if the above fails, try this:
# remotes::install_local('./BernhardtMetabolism_0.1.1.zip')
```

# 

------------------------------------------------------------------------

## Filtering and gap-filling the dataset

------------------------------------------------------------------------

<font size="6">Filtering and gap-filling the dataset</font>

### 1. Calculate diagnostics for each site (both lotic and terrestrial)

To help filter down the full standardized datasets, a series of
diagnostic measures were calculated for each site-year of data. These
diagnostics include information about model fits and the availability of
various data sources. These diagnostics differ slightly between lotic
and terrestrial sites due to differences in how they were modeled, and
because additional data sources were used to further investigate
patterns in lotic metabolism.

#### 1.1 Calculate diagnostics for lotic sites

For each site-year the following set of diagnostics were calculated with
the help of the <span style="color:orange">**diagnostics\_fun**</span>
from the included **BernhardtMetabolism** package.

The output contains the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Year”:</span>** Year
-   **<span style="color:#009688">“ER\_K”:</span>** The correlation
    coefficient (R<sup>2</sup>) between ecosystem respiration and K600.
    This was calculated from the “ER” column and thus does not include
    positive ER values as part of the calculation.
-   **<span style="color:#009688">“GPP\_neg”:</span>** Percentage of
    days with negative stream GPP estimates based on the raw GPP data
    (“GPP\_raw”)
-   **<span style="color:#009688">“ER\_pos”:</span>** Percentage of days
    with positive stream ER estimates based on the raw ER data
    ("ER\_raw)
-   **<span style="color:#009688">“num\_days”:</span>** The number of
    days with estimated stream metabolism (excluding days with negative
    GPP and positive ER)
-   **<span style="color:#009688">“NLDAS\_data”:</span>** Logical
    indicator (TRUE or FALSE) for availability of NLDAS data for this
    site-year
-   **<span style="color:#009688">“LAI\_data”:</span>** Logical
    indicator (TRUE or FALSE) for availability of LAI data for this
    site-year
-   **<span style="color:#009688">“Predicted\_light”:</span>** Logical
    indicator (TRUE or FALSE) for availability of StreamLight
    predictions of light at the stream surface for this site-year
-   **<span style="color:#009688">“MODIS\_NPP”:</span>** Logical
    indicator (TRUE or FALSE) for availability of MODIS annual NPP data
    for this site-year

#### 1.2 Calculate diagnostics for terrestrial sites

For each site-year the following set of diagnostics were calculated with
the help of the <span style="color:orange">**diagnostics\_fun**</span>
from the included **BernhardtMetabolism** package.

The output contains the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Year”:</span>** Year
-   **<span style="color:#009688">“GPP\_neg”:</span>** Percentage of
    days with negative GPP estimates based on the raw GPP data
    (“GPP\_raw”)
-   **<span style="color:#009688">“ER\_pos”:</span>** Percentage of days
    with positive ER estimates based on the raw ER data (“ER\_raw”)
-   **<span style="color:#009688">“num\_days”:</span>** The number of
    days (excluding days with negative GPP and positive ER)

### 2. Filtering and gap-filling the datasets

#### 2.1 Lotic sites

##### 2.1.1 Filter the lotic metabolism dataset

At this point all of the data for any site with metabolism data in the
CONUS has been compiled. Analysis for the metabolism synthesis paper
will be derived from a filtered subset of this dataset based on the
previously calculated data diagnostic measures.

<font size="3">Filtering criteria:</font>

-   Site-years must have stream metabolism estimates (excluding days
    with negative GPP estimates or positive ER estimates) for at least
    60% of the year
-   R<sup>2</sup> between ER and K600 less than 0.6 for the site-year
-   The maximum K600 value is less than 100 for the site-year
-   The site-year has NLDAS incoming shortwave radiation data (used to
    calculate annual total incoming PAR)
-   The site-year has light estimates at the stream surface
-   The site-year has MODIS terrestrial NPP data available

``` r
#===============================================================================
#Filter the lotic metabolism data
#===============================================================================
#Filter the site-years
  lotic_filtered <- lotic_yearly_diagnostics %>%
    filter(num_days >= 365 * 0.6) %>%
    filter(ER_K < 0.6) %>%
    filter(K600_max < 100) %>%
    filter(NLDAS_data == TRUE) %>%
    filter(Predicted_light == TRUE) %>%
    filter(modis_npp == TRUE)

#Retrieve site-years of timeseries from the standardized metabolism dataset
  lotic_filtered_metabolism <- BernhardtMetabolism::filter_metab_data(
    filtered_sy = lotic_filtered,
    standardized_timeseries = lotic_standardized_full
  )
```

    ## A total of 222 sites and 865 site-years meet the filtering criteria

##### 2.1.2 Gap-fill lotic sites

To calculate annual sums, as well as several descriptive statistics, it
was necessary to first perform gap-filling. GPP, ER, NEP, water
temperature, daily incoming PAR, and discharge were all gap-filled.
Gap-filling was done using a structural timeseries. The code for
gap-filling is a slight modification of the <span
style="color:orange">**fillMiss**</span> function from the
[**waterData** package](https://github.com/USGS-R/waterData).

At this stage several notable steps occurred:

-   Stream NEP (GPP + ER) was calculated
-   Additional columns for stream GPP, ER, and NEP expressed in units of
    (g C m<sup>-2</sup> d<sup>-1</sup>) were calculated based on a
    photosynthetic quotient (PQ) of 1.25.
-   Gap-filling of stream metabolism, water temperature, daily incoming
    PAR, discharge
-   Z-normalization of stream metabolism, water temperature, daily
    incoming PAR

``` r
#===============================================================================
#Gap-fill the filtered lotic metabolism
#===============================================================================
#Gap-fill the data       
  lotic_gap_filled <- BernhardtMetabolism::synthesis_gapfill(
    synthesis_filtered = lotic_filtered_metabolism, 
    PQ = 1.25,
    block = 150,
    pmiss = 50
  )   
```

#### 2.2 Terrestrial sites

The terrestrial sites only need to be filtered, but not gap-filled
since: 1.) It is not necessary to perform gap-filling to derive annual
sums because the FLUXNET2015 dataset already provides annual sums, and
2.) because terrestrial sites are not the primary focus of this paper,
only a few derived metrics are calculated and none require gap-filled
data.

<font size="3">Filtering criteria:</font>

-   Site-years must have metabolism estimates (excluding days with
    negative GPP estimates or positive ER estimates) for at least 60% of
    the year
-   Additionally, tier-2 FLUXNET sites were removed based on their
    stricter data usage agreements (“RU-Sam”, “RU-SkP”, “RU-Tks”,
    “RU-Vrk”, “SE-St1”, “ZA-Kru”).

``` r
#===============================================================================
#Filter the FLUXNET metabolism data
#===============================================================================
#Filter the site-years, including removing tier-2 sites, that grant data producers collaboration rights
  fluxnet_filtered <- fluxnet_yearly_diagnostics %>%
    filter(num_days >= 365 * 0.6) %>%
    filter(!(Site_ID %in% c("RU-Sam", "RU-SkP", "RU-Tks",  "RU-Vrk", "SE-St1", "ZA-Kru")))
        
#Retrieve site-years of timeseries from the standardized metabolism dataset
  fluxnet_filtered_metabolism <- BernhardtMetabolism::filter_metab_data(
    filtered_sy = fluxnet_filtered,
    standardized_timeseries = fluxnet_standardized_full
  )  
```

    ## A total of 163 sites and 1143 site-years meet the filtering criteria

------------------------------------------------------------------------

## Calculating site metrics

------------------------------------------------------------------------

<font size="6">Calculating site metrics</font>

### 1. Metrics for lotic synthesis sites

We calculated a number of metrics to describe both metabolism and
environmental conditions for the final subset of filtered and gap-filled
lotic sites. A number of these metrics are based on a set of
ecologically relevant streamflow statistics as described by [Archfield
et
al. (2013)](https://onlinelibrary.wiley.com/doi/abs/10.1002/rra.2710).

<font size="4">“Magnificent 7” metrics from [Archfield et
al. (2013)](https://onlinelibrary.wiley.com/doi/abs/10.1002/rra.2710)</font>

-   Mean
-   Coefficient of variation (CV)
-   Skewness
-   Kurtosis
-   Autoregressive lag-one correleation coefficient (AR1)
-   Amplitude
-   Phase

The “magnificent 7” statistics help to capture variation in the
magnitude, frequency, duration, timing, and rate of change in a
timeseries. Therefore, we extended the initial application of this suite
of metrics to include GPP, ER, water temperature, discharge, PAR, and
LAI. The code used to calculate these statistics is based on code from
the [**EflowStats** package](https://github.com/USGS-R/EflowStats),
which we modified to be generalized across different types of data.

<font size="4">Final compiled table of site information</font>

Finally, the site metrics were joined with the basic site information to
provide a single table of site information and calculated metrics for
the filtered set of lotic sites.

The output contains the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Name”:</span>** Site long name
-   **<span style="color:#009688">“Source”:</span>** Site data source
-   **<span style="color:#009688">“Lat”:</span>** Site Latitude
-   **<span style="color:#009688">“Lon”:</span>** Site Longitude
-   **<span style="color:#009688">“epsg\_crs”:</span>** Site coordinate
    reference system as EPSG code
-   **<span style="color:#009688">“COMID”:</span>** NHDPlus\_v2 reach
    COMID
-   **<span style="color:#009688">“VPU”:</span>** Vector processing unit
-   **<span style="color:#009688">“StreamOrde”:</span>** Stream order
    (from NHDPlus\_v2 flowlines)
-   **<span style="color:#009688">“Azimuth”:</span>** Channel azimuth
    calculated as the circular mean of azimuths for each stream reach
    based on latitude and longitude of the site location using
    NHDPlus\_v2 hydrography data.
-   **<span style="color:#009688">“TH”:</span>** Tree height (m).
    Estimates were derived from Tree heights were derived using 30m
    resolution global canopy height estimates from [Potapov et
    al. (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0034425720305381?via%3Dihub)
-   **<span style="color:#009688">“Width”:</span>** Channel width (m)
-   **<span style="color:#009688">“Width\_src”:</span>** Source of
    channel width estimates. Values include: (NWIS field measurements,
    Regional geomorphic scaling coeff, StreamPULSE estimates)
-   **<span style="color:#009688">“WS\_area\_km2”:</span>** Watershed
    area (km <sup>-2</sup>)
-   **<span style="color:#009688">“WS\_area\_src”:</span>** Source of
    watershed area estimates. Values include: (Appling2018\_USGS2013,
    Appling2018\_StreamStats, nwis\_site\_description, StreamStats,
    localuser\_HBFLTER, localuser\_UNHWQAL)
-   **<span style="color:#009688">“ann\_GPP\_C”:</span>** Mean annual
    cumulative stream GPP (g C m<sup>-2</sup> y<sup>-1</sup>). This was
    calculated by first calculating annual sums of GPP (g C
    m<sup>-2</sup> y<sup>-1</sup>) for each site year, and then taking
    the mean annual rate for each site.
-   **<span style="color:#009688">“upper\_GPP\_C”:</span>** 95th
    percentile of daily rates of stream GPP (g C m<sup>-2</sup>
    d<sup>-1</sup>).
-   **<span style="color:#009688">“ann\_ER\_C”:</span>** Mean annual
    cumulative stream ER (g C m<sup>-2</sup> y<sup>-1</sup>). This was
    calculated by first calculating annual sums of ER (g C
    m<sup>-2</sup> y<sup>-1</sup>) for each site year, and then taking
    the mean annual rate for each site.
-   **<span style="color:#009688">“lower\_ER\_C”:</span>** 5th
    percentile of daily rates of stream ER (g C m<sup>-2</sup>
    d<sup>-1</sup>). Since ER is negative, you can think of this as
    equivalent to the 95th percentile done for GPP.
-   **<span style="color:#009688">“PAR\_sum”:</span>** Mean annual
    cumulative incoming PAR (kmol m<sup>-2</sup> y<sup>-1</sup>). This
    was calculated by first calculating annual sums of PAR (kmol
    m<sup>-2</sup> y<sup>-1</sup>) for each site year, and then taking
    the mean annual rate for each site.
-   **<span style="color:#009688">“Stream\_PAR\_sum”:</span>** Mean
    annual cumulative predicted PAR at the stream surface (kmol
    m<sup>-2</sup> y<sup>-1</sup>). This was calculated by first
    calculating annual sums of predicted PAR at the stream surface (kmol
    m<sup>-2</sup> y<sup>-1</sup>) for each site year, and then taking
    the mean annual rate for each site.
-   **<span style="color:#009688">“gpp\_C\_mean”:</span>** Mean daily
    GPP (g C m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“gpp\_C\_cv”:</span>** CV of daily GPP
-   **<span style="color:#009688">“gpp\_C\_skew”:</span>** Skewness of
    daily GPP
-   **<span style="color:#009688">“gpp\_C\_kurt”:</span>** Kurtosis of
    daily GPP
-   **<span style="color:#009688">“gpp\_C\_amp”:</span>** Amplitude of
    daily GPP
-   **<span style="color:#009688">“gpp\_C\_phase”:</span>** Phase of
    daily GPP (day of year)
-   **<span style="color:#009688">“gpp\_C\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily GPP
-   **<span style="color:#009688">“er\_C\_mean”:</span>** Mean daily ER
    (g C m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“er\_C\_cv”:</span>** CV of daily ER
-   **<span style="color:#009688">“er\_C\_skew”:</span>** Skewness of
    daily ER
-   **<span style="color:#009688">“er\_C\_kurt”:</span>** Kurtosis of
    daily ER
-   **<span style="color:#009688">“er\_C\_amp”:</span>** Amplitude of
    daily ER
-   **<span style="color:#009688">“er\_C\_phase”:</span>** Phase of
    daily ER (day of year)
-   **<span style="color:#009688">“er\_C\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily ER
-   **<span style="color:#009688">"Wtemp\_mean** Mean daily water
    temperature (degreees C)
-   **<span style="color:#009688">“Wtemp\_cv”:</span>** CV of daily GPP
-   **<span style="color:#009688">“Wtemp\_skew”:</span>** Skewness of
    daily water temperature
-   **<span style="color:#009688">“Wtemp\_kurt”:</span>** Kurtosis of
    daily water temperature
-   **<span style="color:#009688">“Wtemp\_amp”:</span>** Amplitude of
    daily water temperature
-   **<span style="color:#009688">“Wtemp\_phase”:</span>** Phase of
    daily water temperature (day of year)
-   **<span style="color:#009688">“Wtemp\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily water temperature
-   **<span style="color:#009688">“Disch\_mean”:</span>** Mean daily
    discharge (m<sup>3</sup> s<sup>-1</sup>)
-   **<span style="color:#009688">“Disch\_cv”:</span>** CV of daily
    discharge
-   **<span style="color:#009688">“Disch\_skew”:</span>** Skewness of
    daily discharge
-   **<span style="color:#009688">“Disch\_kurt”:</span>** Kurtosis of
    daily discharge
-   **<span style="color:#009688">“Disch\_amp”:</span>** Amplitude of
    daily discharge
-   **<span style="color:#009688">“Disch\_phase”:</span>** Phaste of
    daily discharge (day of year)
-   **<span style="color:#009688">“Disch\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily discharge
-   **<span style="color:#009688">“PAR\_mean”:</span>** Mean daily sum
    of incoming above canopy PAR (mol m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“PAR\_cv”:</span>** CV of daily PAR
-   **<span style="color:#009688">“PAR\_skew”:</span>** Skewness of
    daily PAR
-   **<span style="color:#009688">“PAR\_kurt”:</span>** Kurtosis of
    daily PAR
-   **<span style="color:#009688">“PAR\_amp”:</span>** Amplitude of
    daily PAR
-   **<span style="color:#009688">“PAR\_phase”:</span>** Phase of daily
    PAR (day of year)
-   **<span style="color:#009688">“PAR\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily PAR
-   **<span style="color:#009688">“LAI\_mean”:</span>** Mean daily LAI
    (m<sup>2</sup> m<sup>-2</sup>)
-   **<span style="color:#009688">“LAI\_cv”:</span>** CV of LAI
-   **<span style="color:#009688">“LAI\_skew”:</span>** Skewness of LAI
-   **<span style="color:#009688">“LAI\_kurt”:</span>** Kurtosis of LAI
-   **<span style="color:#009688">“LAI\_amp”:</span>** Amplitude of LAI
-   **<span style="color:#009688">“LAI\_phase”:</span>** Phase of LAI
    (day of year)
-   **<span style="color:#009688">“LAI\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily LAI
-   **<span style="color:#009688">“MOD\_ann\_NPP”:</span>** Mean annual
    MODIS NPP (g C m<sup>-2</sup> y<sup>-1</sup>) for the concurrent
    period of record for stream metabolism data at a site. Annual sums
    of NPP (g C m<sup>-2</sup> d<sup>-y</sup>) were available for each
    site year and then the mean was taken to get a mean annual rate for
    each site.
-   **<span style="color:#009688">“ndays”:</span>** Total number of days
    with daily GPP (non gap-filled) for the site in the filtered dataset
-   **<span style="color:#009688">“nyears”:</span>** Total number of
    years for the site in the filtered dataset
-   **<span style="color:#009688">“coverage”:</span>** Total coverage of
    daily GPP (non gap-filled) for the site, calculated as ndays / all
    possible days for all site-years included in the filtered dataset.
    Ranges from 0-1.

### 2. Metrics for terrestrial sites

For the terrestrial sites, only a small number of summary measures were
calculated to capture the magnitude of fluxes of GPP and ER. These were
then joined to some basic site information to provide a single table of
site information and calculated metrics for the filtered set of
terrestrial sites.

The output contains the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Name”:</span>** Site long name
-   **<span style="color:#009688">“Lat”:</span>** Site Latitude
-   **<span style="color:#009688">“Lon”:</span>** Site Longitude
-   **<span style="color:#009688">“ann\_GPP”:</span>** Mean annual
    cumulative GPP (g C m<sup>-2</sup> y<sup>-1</sup>). This was
    calculated from the annual sums of GPP (g C m<sup>-2</sup>
    y<sup>-1</sup>) for each site year provided by FLUXNET, and then
    taking the mean annual rate for each site.
-   **<span style="color:#009688">“upper\_GPP”:</span>** 95th percentile
    of daily rates of GPP (g C m<sup>-2</sup> d<sup>-1</sup>).
-   **<span style="color:#009688">“ann\_ER”:</span>** Mean annual
    cumulative ER (g C m<sup>-2</sup> y<sup>-1</sup>). This was
    calculated from the annual sums of ER (g C m<sup>-2</sup>
    y<sup>-1</sup>) for each site year provided by FLUXNET, and then
    taking the mean annual rate for each site
-   **<span style="color:#009688">“lower\_ER”:</span>** 5th percentile
    of daily rates of ER (g C m<sup>-2</sup> d<sup>-1</sup>). Since ER
    is negative, you can think of this as equivalent to the 95th
    percentile done for GPP.  
-   **<span style="color:#009688">“ndays”:</span>** Total number of days
    with daily GPP (non gap-filled) for the site in the filtered dataset
-   **<span style="color:#009688">“nyears”:</span>** Total number of
    years for the site in the filtered dataset
-   **<span style="color:#009688">“coverage”:</span>** Total coverage of
    daily GPP (non gap-filled) for the site, calculated as ndays / all
    possible days for all site-years included in the filtered dataset.
    Ranges from 0-1.

------------------------------------------------------------------------

## Exporting the final dataset

------------------------------------------------------------------------

<font size="6">Exporting the final dataset</font>

A final set of outputs were then exported to the output\_data directory
of this repository. All files were output as both a .rds and .csv to
facilitate use both within R and other environments. However, this repo
only contains the .rds files and the .csv files can be found as part of
the encapsulated data repository HERE.

### 1. Standardized, unfiltered, lotic metabolism

This is the full, unfiltered, standardized dataset of stream metabolism
that was compiled.

<font size="3"><i><span
style="color:#00cc99">lotic\_site\_info\_full(.rds/.csv)</span></i></font>

Format: A single data frame with the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Name”:</span>** Site long name
-   **<span style="color:#009688">“Source”:</span>** Site data source
-   **<span style="color:#009688">“Lat”:</span>** Site Latitude
-   **<span style="color:#009688">“Lon”:</span>** Site Longitude
-   **<span style="color:#009688">“epsg\_crs”:</span>** Site coordinate
    reference system as EPSG code
-   **<span style="color:#009688">“COMID”:</span>** NHDPlus\_v2 reach
    COMID
-   **<span style="color:#009688">“VPU”:</span>** Vector processing unit
-   **<span style="color:#009688">“StreamOrde”:</span>** Stream order
    (from NHDPlus\_v2 flowlines)
-   **<span style="color:#009688">“Azimuth”:</span>** Channel azimuth
    calculated as the circular mean of azimuths for each stream reach
    based on latitude and longitude of the site location using
    NHDPlus\_v2 hydrography data.
-   **<span style="color:#009688">“TH”:</span>** Tree height (m).
    Estimates were derived from Tree heights were derived using 30m
    resolution global canopy height estimates from [Potapov et
    al. (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0034425720305381?via%3Dihub)
-   **<span style="color:#009688">“Width”:</span>** Channel width (m)
-   **<span style="color:#009688">“Width\_src”:</span>** Source of
    channel width estimates. Values include: (NWIS field measurements,
    Regional geomorphic scaling coeff, StreamPULSE estimates)
-   **<span style="color:#009688">“WS\_area\_km2”:</span>** Watershed
    area (km<sup>2</sup>)
-   **<span style="color:#009688">“WS\_area\_src”:</span>** Source of
    watershed area estimates. Values include: (Appling2018\_USGS2013,
    Appling2018\_StreamStats, nwis\_site\_description, StreamStats,
    localuser\_HBFLTER, localuser\_UNHWQAL)

<font size="3"><i><span
style="color:#00cc99">lotic\_standardized\_full(.rds/.csv)</span></i></font>

Format: A list of data frames, where each element of the list is a data
frame of timeseries for a single site (note the .csv version is one
single data frame that contains all data for all sites). The names of
each list element correspond to the unique site identifier (Site\_ID)
for a site. Each data frame contains the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Date”:</span>** Date in YYYY-MM-DD
    format
-   **<span style="color:#009688">“U\_ID”:</span>** Unique date
    identifier (format as year + DOY)
-   **<span style="color:#009688">“Year”:</span>** Year
-   **<span style="color:#009688">“DOY”:</span>** Day of year (1 to 365
    or 366)
-   **<span style="color:#009688">“GPP\_raw”:</span>** Stream GPP
    estimates (g O2 m<sup>-2</sup> d<sup>-1</sup>) from Appling et
    al. (2018), raw data
-   **<span style="color:#009688">“ER\_raw”:</span>** Stream ER
    estimates (g O2 m<sup>-2</sup> d<sup>-1</sup>) from Appling et
    al. (2018), raw data
-   **<span style="color:#009688">“GPP”:</span>** Stream GPP estimates
    (g O2 m<sup>-2</sup> d<sup>-1</sup>) from Appling et al. (2018),
    negative values replaced with NA
-   **<span style="color:#009688">“ER”:</span>** Stream ER estimates (g
    O2 m<sup>-2</sup> d<sup>-1</sup>) from Appling et al. (2018),
    positive values replaced with NA
-   **<span style="color:#009688">“K600”:</span>** Model estimate of
    K600, the mean reaeration rate coefficient, scaled to a Schmidt
    number of 600, for this date. Value is the median of the post warmup
    MCMC distribution
-   **<span style="color:#009688">“DO.obs”:</span>** Mean dissolved
    oxygen concentration (mg O2 L<sup>-1</sup>) for the date (4am to
    3:59am)
-   **<span style="color:#009688">“DO.sat”:</span>** Mean theoretical
    saturation concentration (mg O2 L<sup>-1</sup>) for the date (4am to
    3:59am)
-   **<span style="color:#009688">“temp.water”:</span>** Mean water
    temperature (degrees C) for the date (4am to 3:59pm)
-   **<span style="color:#009688">“discharge”:</span>** Mean discharge
    (m<sup>3</sup> s<sup>-1</sup>) for the date (4am to 3:59pm)
-   **<span style="color:#009688">“PAR\_sum”:</span>** Daily sum of
    incoming above canopy PAR (mol m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“Stream\_PAR\_sum”:</span>** Daily sum
    of PAR estimated at the stream surface (mol m<sup>-2</sup>
    d<sup>-1</sup>)
-   **<span style="color:#009688">“LAI\_proc”:</span>** MODIS LAI data
    that has been processed and gap-filled (m<sup>2</sup>
    m<sup>-2</sup>)

### 2. Standardized and filtered data used in Bernhardt et al. (2022)

#### 2.1 Standardized, filtered, gap-filled lotic metabolism

This is the dataset of stream metabolism used in Bernhardt et
al. (2022). This data has been filtered based on several measures of
data quality and availability and gaps were gap-filled. Finally, a suite
of site metrics were calculated for use in analysis.

<font size="3"><i><span
style="color:#00cc99">lotic\_site\_info\_filtered(.rds/.csv)</span></i></font>

Format: A single data frame with the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Name”:</span>** Site long name
-   **<span style="color:#009688">“Source”:</span>** Site data source
-   **<span style="color:#009688">“Lat”:</span>** Site Latitude
-   **<span style="color:#009688">“Lon”:</span>** Site Longitude
-   **<span style="color:#009688">“epsg\_crs”:</span>** Site coordinate
    reference system as EPSG code
-   **<span style="color:#009688">“COMID”:</span>** NHDPlus\_v2 reach
    COMID
-   **<span style="color:#009688">“VPU”:</span>** Vector processing unit
-   **<span style="color:#009688">“StreamOrde”:</span>** Stream order
    (from NHDPlus\_v2 flowlines)
-   **<span style="color:#009688">“Azimuth”:</span>** Channel azimuth
    calculated as the circular mean of azimuths for each stream reach
    based on latitude and longitude of the site location using
    NHDPlus\_v2 hydrography data.
-   **<span style="color:#009688">“TH”:</span>** Tree height (m).
    Estimates were derived from Tree heights were derived using 30m
    resolution global canopy height estimates from [Potapov et
    al. (2021)](https://www.sciencedirect.com/science/article/abs/pii/S0034425720305381?via%3Dihub)
-   **<span style="color:#009688">“Width”:</span>** Channel width (m)
-   **<span style="color:#009688">“Width\_src”:</span>** Source of
    channel width estimates. Values include: (NWIS field measurements,
    Regional geomorphic scaling coeff, StreamPULSE estimates)
-   **<span style="color:#009688">“WS\_area\_km2”:</span>** Watershed
    area (km<sup>2</sup>)
-   **<span style="color:#009688">“WS\_area\_src”:</span>** Source of
    watershed area estimates. Values include: (Appling2018\_USGS2013,
    Appling2018\_StreamStats, nwis\_site\_description, StreamStats,
    localuser\_HBFLTER, localuser\_UNHWQAL)
-   **<span style="color:#009688">“ann\_GPP\_C”:</span>** Mean annual
    cumulative stream GPP (g C m<sup>-2</sup> y<sup>-1</sup>). This was
    calculated by first calculating annual sums of GPP (g C
    m<sup>-2</sup> y<sup>-1</sup>) for each site year, and then taking
    the mean annual rate for each site.
-   **<span style="color:#009688">“upper\_GPP\_C”:</span>** 95th
    percentile of daily rates of stream GPP (g C m<sup>-2</sup>
    d<sup>-1</sup>).
-   **<span style="color:#009688">“ann\_ER\_C”:</span>** Mean annual
    cumulative stream ER (g C m<sup>-2</sup> y<sup>-1</sup>). This was
    calculated by first calculating annual sums of ER (g C
    m<sup>-2</sup> y<sup>-1</sup>) for each site year, and then taking
    the mean annual rate for each site.
-   **<span style="color:#009688">“lower\_ER\_C”:</span>** 5th
    percentile of daily rates of stream ER (g C m<sup>-2</sup>
    d<sup>-1</sup>). Since ER is negative, you can think of this as
    equivalent to the 95th percentile done for GPP.
-   **<span style="color:#009688">“PAR\_sum”:</span>** Mean annual
    cumulative incoming PAR (kmol m<sup>-2</sup> y<sup>-1</sup>). This
    was calculated by first calculating annual sums of PAR (kmol
    m<sup>-2</sup> y<sup>-1</sup>) for each site year, and then taking
    the mean annual rate for each site.
-   **<span style="color:#009688">“Stream\_PAR\_sum”:</span>** Mean
    annual cumulative predicted PAR at the stream surface (kmol
    m<sup>-2</sup> y<sup>-1</sup>). This was calculated by first
    calculating annual sums of predicted PAR at the stream surface (kmol
    m<sup>-2</sup> y<sup>-1</sup>) for each site year, and then taking
    the mean annual rate for each site.
-   **<span style="color:#009688">“gpp\_C\_mean”:</span>** Mean daily
    GPP (g C m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“gpp\_C\_cv”:</span>** CV of daily GPP
-   **<span style="color:#009688">“gpp\_C\_skew”:</span>** Skewness of
    daily GPP
-   **<span style="color:#009688">“gpp\_C\_kurt”:</span>** Kurtosis of
    daily GPP
-   **<span style="color:#009688">“gpp\_C\_amp”:</span>** Amplitude of
    daily GPP
-   **<span style="color:#009688">“gpp\_C\_phase”:</span>** Phase of
    daily GPP (day of year)
-   **<span style="color:#009688">“gpp\_C\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily GPP
-   **<span style="color:#009688">“er\_C\_mean”:</span>** Mean daily ER
    (g C m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“er\_C\_cv”:</span>** CV of daily ER
-   **<span style="color:#009688">“er\_C\_skew”:</span>** Skewness of
    daily ER
-   **<span style="color:#009688">“er\_C\_kurt”:</span>** Kurtosis of
    daily ER
-   **<span style="color:#009688">“er\_C\_amp”:</span>** Amplitude of
    daily ER
-   **<span style="color:#009688">“er\_C\_phase”:</span>** Phase of
    daily ER (day of year)
-   **<span style="color:#009688">“er\_C\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily ER
-   **<span style="color:#009688">"Wtemp\_mean** Mean daily water
    temperature (degreees C)
-   **<span style="color:#009688">“Wtemp\_cv”:</span>** CV of daily GPP
-   **<span style="color:#009688">“Wtemp\_skew”:</span>** Skewness of
    daily water temperature
-   **<span style="color:#009688">“Wtemp\_kurt”:</span>** Kurtosis of
    daily water temperature
-   **<span style="color:#009688">“Wtemp\_amp”:</span>** Amplitude of
    daily water temperature
-   **<span style="color:#009688">“Wtemp\_phase”:</span>** Phase of
    daily water temperature (day of year)
-   **<span style="color:#009688">“Wtemp\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily water temperature
-   **<span style="color:#009688">“Disch\_mean”:</span>** Mean daily
    discharge (m<sup>3</sup> s<sup>-1</sup>)
-   **<span style="color:#009688">“Disch\_cv”:</span>** CV of daily
    discharge
-   **<span style="color:#009688">“Disch\_skew”:</span>** Skewness of
    daily discharge
-   **<span style="color:#009688">“Disch\_kurt”:</span>** Kurtosis of
    daily discharge
-   **<span style="color:#009688">“Disch\_amp”:</span>** Amplitude of
    daily discharge
-   **<span style="color:#009688">“Disch\_phase”:</span>** Phaste of
    daily discharge (day of year)
-   **<span style="color:#009688">“Disch\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily discharge
-   **<span style="color:#009688">“PAR\_mean”:</span>** Mean daily sum
    of incoming above canopy PAR (mol m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“PAR\_cv”:</span>** CV of daily PAR
-   **<span style="color:#009688">“PAR\_skew”:</span>** Skewness of
    daily PAR
-   **<span style="color:#009688">“PAR\_kurt”:</span>** Kurtosis of
    daily PAR
-   **<span style="color:#009688">“PAR\_amp”:</span>** Amplitude of
    daily PAR
-   **<span style="color:#009688">“PAR\_phase”:</span>** Phase of daily
    PAR (day of year)
-   **<span style="color:#009688">“PAR\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily PAR
-   **<span style="color:#009688">“LAI\_mean”:</span>** Mean daily LAI
    (m<sup>2</sup> m<sup>-2</sup>)
-   **<span style="color:#009688">“LAI\_cv”:</span>** CV of LAI
-   **<span style="color:#009688">“LAI\_skew”:</span>** Skewness of LAI
-   **<span style="color:#009688">“LAI\_kurt”:</span>** Kurtosis of LAI
-   **<span style="color:#009688">“LAI\_amp”:</span>** Amplitude of LAI
-   **<span style="color:#009688">“LAI\_phase”:</span>** Phase of LAI
    (day of year)
-   **<span style="color:#009688">“LAI\_ar1”:</span>** Autoregressive
    lag-one correlation coefficient of daily LAI
-   **<span style="color:#009688">“MOD\_ann\_NPP”:</span>** Mean annual
    MODIS NPP (g C m<sup>-2</sup> y<sup>-1</sup>) for the concurrent
    period of record for stream metabolism data at a site. Annual sums
    of NPP (g C m<sup>-2</sup> y<sup>-1</sup>) were available for each
    site year and then the mean was taken to get a mean annual rate for
    each site.
-   **<span style="color:#009688">“ndays”:</span>** Total number of days
    with daily GPP (non gap-filled) for the site in the filtered dataset
-   **<span style="color:#009688">“nyears”:</span>** Total number of
    years for the site in the filtered dataset
-   **<span style="color:#009688">“coverage”:</span>** Total coverage of
    daily GPP (non gap-filled) for the site, calculated as ndays / all
    possible days for all site-years included in the filtered dataset.
    Ranges from 0-1.

<font size="3"><i><span
style="color:#00cc99">lotic\_gap\_filled(.rds/.csv)</span></i></font>

Format: A list of data frames, where each element of the list is a data
frame of timeseries for a single site (note the .csv version is one
single data frame that contains all data for all sites). The names of
each list element correspond to the unique site identifier (Site\_ID)
for a site. Each data frame contains the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Date”:</span>** Date in YYYY-MM-DD
    format
-   **<span style="color:#009688">“U\_ID”:</span>** Unique date
    identifier (format as year + DOY)
-   **<span style="color:#009688">“Year”:</span>** Year
-   **<span style="color:#009688">“DOY”:</span>** Day of year (1 to 365
    or 366)
-   **<span style="color:#009688">“GPP\_raw”:</span>** Stream GPP
    estimates (g O2 m<sup>-2</sup> d<sup>-1</sup>) from Appling et
    al. (2018), raw data
-   **<span style="color:#009688">“ER\_raw”:</span>** Stream ER
    estimates (g O2 m<sup>-2</sup> d<sup>-1</sup>) from Appling et
    al. (2018), raw data
-   **<span style="color:#009688">“GPP”:</span>** Stream GPP estimates
    (g O2 m<sup>-2</sup> d<sup>-1</sup>) from Appling et al. (2018),
    negative values replaced with NA
-   **<span style="color:#009688">“ER”:</span>** Stream ER estimates (g
    O2 m<sup>-2</sup> d<sup>-1</sup>) from Appling et al. (2018),
    positive values replaced with NA
-   **<span style="color:#009688">“K600”:</span>** Model estimate of
    K600, the mean reaeration rate coefficient, scaled to a Schmidt
    number of 600, for this date. Value is the median of the post warmup
    MCMC distribution
-   **<span style="color:#009688">“DO.obs”:</span>** Mean dissolved
    oxygen concentration (mg O2 L<sup>-1</sup>) for the date (4am to
    3:59am)
-   **<span style="color:#009688">“DO.sat”:</span>** Mean theoretical
    saturation concentration (mg O2 L<sup>-1</sup>) for the date (4am to
    3:59am)
-   **<span style="color:#009688">“temp.water”:</span>** Mean water
    temperature (degrees C) for the date (4am to 3:59pm)
-   **<span style="color:#009688">“discharge”:</span>** Mean discharge
    (m<sup>3</sup> s<sup>-1</sup>) for the date (4am to 3:59pm)
-   **<span style="color:#009688">“PAR\_sum”:</span>** Daily sum of
    incoming above canopy PAR (mol m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“Stream\_PAR\_sum”:</span>** Daily sum
    of PAR estimated at the stream surface (mol m<sup>-2</sup>
    d<sup>-1</sup>)
-   **<span style="color:#009688">“LAI\_proc”:</span>** MODIS LAI data
    that has been processed and gap-filled (m<sup>2</sup>
    m<sup>-2</sup>)
-   **<span style="color:#009688">“GPP\_filled”:</span>** Gap-filled
    stream GPP estimates (g O2 m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“ER\_filled”:</span>** Gap-filled
    stream ER estimates (g O2 m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“NEP\_filled”:</span>** Gap-filled
    stream NEP estimates (g O2 m<sup>-2</sup> d<sup>-1</sup>)
-   **<span style="color:#009688">“GPP\_C\_filled”:</span>** Gap-filled
    stream GPP estimates, expressed in carbon (g C m<sup>-2</sup>
    d<sup>-1</sup>)
-   **<span style="color:#009688">“ER\_C\_filled”:</span>** Gap-filled
    stream ER estimates, expressed in carbon (g C m<sup>-2</sup>
    d<sup>-1</sup>)
-   **<span style="color:#009688">“NEP\_c\_filled”:</span>** Gap-filled
    stream NEP estimates, expressed in carbon (g C m<sup>-2</sup>
    d<sup>-1</sup>)
-   **<span style="color:#009688">“Wtemp\_filled”:</span>** Gap-filled
    mean water temperature from the “temp.water” column.
-   **<span style="color:#009688">“Disch\_filled”:</span>** Gap-filled
    mean discharge from the “discharge” column
-   **<span style="color:#009688">“PAR\_filled”:</span>** Gap-filled
    daily sum of incoming above canopy PAR, from the "PAR\_sum column
-   **<span style="color:#009688">“GPP\_norm”:</span>** Z-normalized GPP
-   **<span style="color:#009688">“ER\_norm”:</span>** Z-normalized ER
-   **<span style="color:#009688">“NEP\_norm”:</span>** Z-normalized NEP
-   **<span style="color:#009688">“Wtemp\_norm”:</span>** Z-normalized
    water temperature
-   **<span style="color:#009688">“PAR\_norm”:</span>** Z-normalized
    daily sum of incoming above canopy PAR

#### 2.2 Standardized and filtered FLUXNET

This is the dataset of terrestrial ecosystem fluxes used in Bernhardt et
al. (2022). This data has been filtered based on data availability and
several basic site metrics were calculated for use in analysis.

<font size="3"><i><span
style="color:#00cc99">fluxnet\_site\_info\_filtered(.rds/.csv)</span></i></font>

Format: A single data frame with the following columns:

-   **<span style="color:#009688">“Site\_ID”:</span>** Unique site
    identifier
-   **<span style="color:#009688">“Name”:</span>** Site long name
-   **<span style="color:#009688">“Lat”:</span>** Site Latitude
-   **<span style="color:#009688">“Lon”:</span>** Site Longitude
-   **<span style="color:#009688">“ann\_GPP”:</span>** Mean annual
    cumulative GPP (g C m<sup>-2</sup> y<sup>-1</sup>). This was
    calculated from the annual sums of GPP (g C m<sup>-2</sup>
    y<sup>-1</sup>) for each site year provided by FLUXNET, and then
    taking the mean annual rate for each site.
-   **<span style="color:#009688">“upper\_GPP”:</span>** 95th percentile
    of daily rates of GPP (g C m<sup>-2</sup> d<sup>-1</sup>).
-   **<span style="color:#009688">“ann\_ER”:</span>** Mean annual
    cumulative ER (g C m<sup>-2</sup> y<sup>-1</sup>). This was
    calculated from the annual sums of ER (g C m<sup>-2</sup>
    y<sup>-1</sup>) for each site year provided by FLUXNET, and then
    taking the mean annual rate for each site
-   **<span style="color:#009688">“lower\_ER”:</span>** 5th percentile
    of daily rates of ER (g C m<sup>-2</sup> d<sup>-1</sup>). Since ER
    is negative, you can think of this as equivalent to the 95th
    percentile done for GPP.
-   **<span style="color:#009688">“ndays”:</span>** Total number of days
    with daily GPP (non gap-filled) for the site in the filtered dataset
-   **<span style="color:#009688">“nyears”:</span>** Total number of
    years for the site in the filtered dataset
-   **<span style="color:#009688">“coverage”:</span>** Total coverage of
    daily GPP (non gap-filled) for the site, calculated as ndays / all
    possible days for all site-years included in the filtered dataset.
    Ranges from 0-1.

<font size="3"><i><span
style="color:#00cc99">fluxnet\_filtered\_metabolism(.rds/.csv)</span></i></font>

Format: A list of data frames, where each element of the list is a data
frame of timeseries for a single site (note the .csv version is one
single data frame that contains all data for all sites). The names of
each list element correspond to the unique site identifier (Site\_ID)
for a site. Each data frame contains the following columns:

-   **<span style="color:#009688">“Date”:</span>** Date in YYYY-MM-DD
    format
-   **<span style="color:#009688">“U\_ID”:</span>** Unique date
    identifier (format as year + DOY)
-   **<span style="color:#009688">“Year”:</span>** Year
-   **<span style="color:#009688">“DOY”:</span>** Day of year (1 to 365
    or 366)
-   **<span style="color:#009688">“GPP\_raw”:</span>** FLUXNET annual
    GPP (sum from daily estimates) (g C m<sup>-2</sup> y<sup>-1</sup>)
    “GPP\_NT\_VUT\_REF”, raw data. Gross Primary Production, from
    Nighttime partitioning method, reference version selected from GPP
    versions using a model efficiency approach.
-   **<span style="color:#009688">“ER\_raw”:</span>** FLUXNET annual ER
    (sum from daily estimates) (g C m<sup>-2</sup> y<sup>-1</sup>)
    “RECO\_NT\_VUT\_REF”, raw data. Ecosystem Respiration, from
    Nighttime partitioning method, reference selected from RECO versions
    using a model efficiency approach.
-   **<span style="color:#009688">“GPP”:</span>** FLUXNET annual GPP
    (sum from daily estimates) (g C m<sup>-2</sup> y<sup>-1</sup>)
    “GPP\_NT\_VUT\_REF”, negative values replaced with NA. Gross Primary
    Production, from Nighttime partitioning method, reference version
    selected from GPP versions using a model efficiency approach.
-   **<span style="color:#009688">“ER”:</span>** FLUXNET annual ER (sum
    from daily estimates) (g C m<sup>-2</sup> y<sup>-1</sup>)
    “RECO\_NT\_VUT\_REF”, positive values replaced with NA. Ecosystem
    Respiration, from Nighttime partitioning method, reference selected
    from RECO versions using a model efficiency approach.
-   **<span style="color:#009688">“Net”:</span>** FLUXNET NEE (changed
    to NEP by \* -1) “NEE\_VUT\_REF”. I derived this from data of this
    description:Net Ecosystem Exchange, using Variable Ustar Threshold
    (VUT) for each year, reference selected on the basis of the model
    efficiency.
-   **<span style="color:#009688">“Temp”:</span>** Average air
    temperature from daily data (degrees C)
-   **<span style="color:#009688">“Precip”:</span>** Average precipition
    from daily data (mm)
-   **<span style="color:#009688">“VPD”:</span>** Average vapor pressure
    deficit from daily data (hPa)
-   **<span style="color:#009688">“SW”:</span>** Average incoming
    shortwave radiation from daily data (W m^-2)

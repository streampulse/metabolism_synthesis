#mike vlah (michael.vlah@duke.edu)
#2020-10-31
#datasets and plots for Bernhardt, Savoy et al. metabolism synthesis paper

#NOTE: lips-specific datasets are generated in section 7. be sure to apply
#   any processing changes made in sections 2-3 to those as well.

library(plyr)
library(tidyverse)
library(RColorBrewer)
library(plotrix)
library(nhdR)

setwd('~/git/streampulse/metab_synthesis_2/')

# 0: setup ####

dir.create('figures',
           showWarnings = FALSE)
dir.create('export_datasets',
           showWarnings = FALSE)

#choose colors for plot categories
fnetcolor = 'sienna4'
spcolor = 'cadetblue4'
gppcolor = 'forestgreen'
ercolor = 'sienna'

source('plot_helpers.R')

# 1: load site data that will be used in plots (as opposed to stats) ####
site_data_1 = readRDS('output/synthesis_site_metrics.rds') %>%
    as_tibble() %>%
    select(sitecode = Site_ID, Stream_PAR_sum, Disch_ar1, MOD_ann_NPP)

# OBSOLETE: load and prepare fluxnet data (highlighting discrepancy) ####

fnet_diag = readRDS('output/FLUXNET_yearly_diagnostics.rds') %>%
    as_tibble()

#summarize by site, starting with a dataset that's summarized by site-year
fnet_ann_2 = readRDS('output/FLUXNET_annual_compiled.rds')
fnet_names_2 = names(fnet_ann_2)

for(i in 1:length(fnet_ann_2)){
    fnet_ann_2[[i]]$sitecode = fnet_names_2[i]
}

fnet_full_2 = Reduce(bind_rows, fnet_ann_2) %>%
    as_tibble()

fnet_site_2 = fnet_full_2 %>%
    select(sitecode, GPP, ER, Net, Year) %>%
    left_join(fnet_diag, by = c('sitecode' = 'Site_ID', 'Year')) %>%
    filter(num_days / 365 >= 0.6) %>% #coverage filter
    group_by(sitecode) %>%
    summarize(GPP = mean(GPP, na.rm = TRUE),
              ER = mean(ER, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(sitecode) %>%
    rename(GPP_C_filled = GPP, ER_C_filled = ER) %>%
    mutate(NEP_C_filled = GPP_C_filled + ER_C_filled)

#summarize by site, starting with a dataset that's summarized by site-year-DOY
fnet_list = readRDS('output/FLUXNET_filtered.rds')

fnet_names = names(fnet_list)
for(i in 1:length(fnet_list)){
    fnet_list[[i]]$sitecode = fnet_names[i]
}

fnet_full = Reduce(bind_rows, fnet_list) %>%
    as_tibble()

fnet_site = fnet_full %>%
    select(sitecode, GPP, ER, DOY, Year) %>%
    group_by(sitecode, Year) %>%
    summarize(GPP_ann_sum = sum(GPP, na.rm = TRUE),
              ER_ann_sum = sum(ER, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(sitecode) %>%
    summarize(GPP_site_mean = mean(GPP_ann_sum, na.rm = TRUE),
              ER_site_mean = mean(ER_ann_sum, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(sitecode) %>%
    mutate(NEP_site_mean = GPP_site_mean + ER_site_mean)

#compare the two summaries (CHECK: should these be identical?
#   FLUXNET_annual_compiled.rds and FLUXNET_filtered.rds don't produce the
#   same result when summarized by site.
#   maybe FLUXNET_annual_compiled.rds wasn't gapfilled?)
head(fnet_site, 3)
head(fnet_site_2, 3)
readRDS('output/FLUXNET_site_metrics.rds') %>%
    as_tibble()

# 2: load and prepare fluxnet data ####

fnet_diag = readRDS('output/FLUXNET_yearly_diagnostics.rds') %>%
    as_tibble()

restricted_use_sites = c("RU-Sam", "RU-SkP", "RU-Tks", "RU-Vrk", "SE-St1", "ZA-Kru")

#summarize by site, starting with a dataset that's summarized by site-year
fnet_ann = readRDS('output/FLUXNET_annual_compiled.rds')
fnet_names = names(fnet_ann)

for(i in 1:length(fnet_ann)){
    fnet_ann[[i]]$sitecode = fnet_names[i]
}

fnet_ann = fnet_ann[! fnet_names %in% restricted_use_sites]

fnet_full = Reduce(bind_rows, fnet_ann) %>%
    as_tibble()

# #REPLACED BY NEXT CHUNK. RESULTS VERIFIED IDENTICAL
# fnet_site = fnet_full %>%
#     select(sitecode, GPP, ER, Net, Year) %>%
#     left_join(fnet_diag, by = c('sitecode' = 'Site_ID', 'Year')) %>%
#     filter(num_days / 365 >= 0.6) %>% #coverage filter
#     group_by(sitecode) %>%
#     summarize(GPP = mean(GPP, na.rm = TRUE),
#               ER = mean(ER, na.rm = TRUE)) %>%
#     ungroup() %>%
#     arrange(sitecode) %>%
#     rename(GPP_site_mean = GPP, ER_site_mean = ER) %>%
#     mutate(NEP_site_mean = GPP_site_mean + ER_site_mean)

fnet_site = readRDS('output/FLUXNET_site_metrics.rds') %>%
    as_tibble() %>%
    select(sitecode = Site_ID,
           GPP_site_mean = ann_GPP,
           ER_site_mean = ann_ER) %>%
    filter(! is.na(GPP_site_mean),
           ! is.na(ER_site_mean)) %>%
    mutate(NEP_site_mean = GPP_site_mean + ER_site_mean) %>%
    arrange(sitecode)

#summarize by DOY for lips plots
fnet_byday_list = readRDS('output/FLUXNET_filtered.rds')

fnet_byday_names = names(fnet_byday_list)
for(i in 1:length(fnet_byday_list)){
    fnet_byday_list[[i]]$sitecode = fnet_byday_names[i]
}

fnet_byday_full = Reduce(bind_rows, fnet_byday_list) %>%
    as_tibble()

fnet_lips = fnet_byday_full %>%
    select(sitecode, GPP, ER, DOY, Year) %>%
    group_by(sitecode, DOY) %>%
    summarize(GPP = mean(GPP, na.rm=TRUE), #average metab by day across years
              ER = mean(ER, na.rm=TRUE)) %>%
    # nyear = length(unique(Year))) %>%
    ungroup() %>%
    arrange(sitecode, DOY) %>%
    rename(GPP_C_filled = GPP, ER_C_filled = ER) %>%
    mutate(NEP_C_filled = GPP_C_filled + ER_C_filled)


# 3. load and prepare streampulse data ####

sp_list = readRDS('output/synthesis_gap_filled.rds')

# #REPLACED BY NEXT CHUNK. RESULTS VERIFIED IDENTICAL
#
# #summarize by site
# sp_names = names(sp_list)
# for(i in 1:length(sp_list)){
#
#     #create sitecode column
#     sp_list[[i]]$sitecode = sp_names[i]
#
#     #create columns for ER-K600 correlation significance and max K600
#     sp_list[[i]] = sp_list[[i]] %>%
#         group_by(Year) %>%
#         mutate(
#             ER_K600_R2 = erk_r2(ER, K600),
#             max_K600 = max(K600, na.rm = TRUE)) %>%
#         ungroup()
# }
#
# sp_full = Reduce(bind_rows, sp_list) %>%
#     as_tibble()

#summarize by site
sp_names = names(sp_list)

for(i in 1:length(sp_list)){
    sp_list[[i]]$sitecode = sp_names[i]
}

sp_intermediate = Reduce(bind_rows, sp_list) %>%
    as_tibble()

sp_full = readRDS('output/lotic_yearly_diagnostics.rds') %>%
    as_tibble() %>%
    select(sitecode = Site_ID,
           Year,
           ER_K,
           max_K600 = K600_max) %>%
    right_join(sp_intermediate,
               by = c('sitecode', 'Year')) %>%
    arrange(sitecode, Date) %>%
    select(Date, U_ID, Year, DOY, GPP:discharge, PAR_sum, Stream_PAR_sum,
           LAI_proc, GPP_filled:PAR_norm, sitecode, ER_K600_R2 = ER_K, max_K600)

sp_site = sp_full %>%
    filter( #remove spurious model results
        ER_K600_R2 < 0.6,
        max_K600 < 100) %>%
    select(sitecode, Year, DOY, GPP_C_filled, ER_C_filled) %>%
    group_by(sitecode, Year) %>%
    summarize(GPP_ann_sum = sum(GPP_C_filled, na.rm = TRUE),
              ER_ann_sum = sum(ER_C_filled, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(sitecode) %>%
    summarize(GPP_site_mean = mean(GPP_ann_sum, na.rm = TRUE),
              ER_site_mean = mean(ER_ann_sum, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(sitecode) %>%
    mutate(NEP_site_mean = GPP_site_mean + ER_site_mean) %>%
    left_join(site_data_1, by = 'sitecode')

#summarize by DOY for lips plots
sp_lips = sp_full %>%
    filter( #remove spurious model results
        ER_K600_R2 < 0.6,#) %>%
        max_K600 < 100) %>%
    mutate(GPP_C_filled = case_when(GPP_C_filled < 0 ~ 0,
                                    TRUE ~ GPP_C_filled),
           ER_C_filled = case_when(ER_C_filled > 0 ~ 0,
                                   TRUE ~ ER_C_filled)) %>%
    select(sitecode, Year, DOY, GPP_C_filled, ER_C_filled) %>%
    group_by(sitecode, DOY) %>%
    summarize(GPP_C_filled = mean(GPP_C_filled, na.rm=TRUE), #average metab by day across years
              ER_C_filled = mean(ER_C_filled, na.rm=TRUE)) %>%
    # nyear = length(unique(Year))) %>%
    ungroup() %>%
    arrange(sitecode, DOY) %>%
    mutate(NEP_C_filled = GPP_C_filled + ER_C_filled) %>%
    left_join(site_data_1, by = 'sitecode') #include site data

# TESTING ONLY: Er * K600 R^2 discrepancy ####

source('R/functions/filter_metab.R')
synthesis_filtered <- filter_metab(
    diag = readRDS("output/lotic_yearly_diagnostics.rds"),
    filters = paste0(paste0("num_days >= ", 365 * 0.6), " & ER_K < 0.6", " & K600_max < 100"),
    metab_rds = readRDS("output/lotic_standardized_metabolism.rds"),
    has_NLDAS = TRUE,
    has_MODIS_NPP = TRUE
)
sp_names_2 = names(synthesis_filtered)
for(i in 1:length(synthesis_filtered)){
    synthesis_filtered[[i]]$sitecode = sp_names[i]
    synthesis_filtered[[i]] = synthesis_filtered[[i]] %>%
        group_by(Year) %>%
        mutate(
            ER_K600_R2XX = erk_r2(ER, K600),
            max_K600XX = max(K600, na.rm = TRUE)) %>%
        ungroup()
}
sp_full_2 = Reduce(bind_rows, synthesis_filtered) %>%
    as_tibble()
select(sp_full, sitecode, Year, ER_K600_R2, max_K600) %>% distinct(sitecode, Year, .keep_all = T)
select(sp_full_2, sitecode, Year, ER_K600_R2XX, max_K600XX) %>% distinct(sitecode, Year, .keep_all = T)
distinct(sp_full, sitecode, Year) %>% group_by(sitecode) %>% summarize(n()) %>% as.data.frame()
distinct(sp_full_2, sitecode, Year) %>% group_by(sitecode) %>% summarize(n()) %>% as.data.frame()
group_by(sp_full_2, sitecode) %>% summarize(length(unique(ER_K600_R2XX)))
group_by(sp_full, sitecode) %>% summarize(length(unique(ER_K600_R2)))
lyd = readRDS("output/lotic_yearly_diagnostics.rds") %>% as_tibble()

#here's what we'd lose by filtering wonky model results:
sp_full %>%
    filter(ER_K600_R2 >= 0.6 | max_K600 > 100) %>%
    select(Year, sitecode, ER_K600_R2, max_K600) %>%
    group_by(sitecode, Year) %>%
    summarize(across(everything(), first),
              .groups = 'drop') %>%
    arrange(sitecode, Year) %>%
    mutate(across(all_of(c('ER_K600_R2', 'max_K600')),
                  ~round(., 2))) %>%
    as.data.frame()

# 4: split sites by coverage ####

sp_coverage = tapply(sp_full$GPP_C, sp_full$sitecode, function(x){
    round(sum(! is.na(x)) / length(x), 2)
})

fnet_coverage = tapply(fnet_full$GPP, fnet_full$sitecode, function(x){
    round(sum(! is.na(x)) / length(x), 2)
})

sp_high_cov_sites = names(which(sp_coverage >= 0.8))
sp_high_cov_bool = sp_site$sitecode %in% sp_high_cov_sites

fnet_high_cov_sites = names(which(fnet_coverage >= 0.8))
fnet_high_cov_bool = fnet_site$sitecode %in% fnet_high_cov_sites

# sp_high_cov_sites = sp$sitecode[sp$coverage >= 0.8]
# sp_high_cov_bool = sp$sitecode %in% sp_high_cov_sites
# fnet_high_cov_sites = fnet$sitecode[fnet$coverage >= 0.8]
# fnet_high_cov_bool = fnet$sitecode %in% fnet_high_cov_sites

# 5: (skip if only plotting) generate site data that will be used for stats ####

site_data_2 = readRDS('output/lotic_site_info.rds') %>%
    as_tibble() %>%
    select(sitecode = Site_ID,
           lat = Lat,
           lon = Lon,
           epsg_crs, COMID, VPU)

#get reach proportional distance for as many sites as possible
errflag = FALSE
site_data_2$reach_proportion = NA
ni = nrow(site_data_2)
for(i in 1:ni){

    print(paste(i, ni, sep = '/'))
    tryCatch({
        rp = calc_reach_prop(VPU = site_data_2$VPU[i],
                             COMID = site_data_2$COMID[i],
                             lat = site_data_2$lat[i],
                             long = site_data_2$lon[i],
                             CRS = site_data_2$epsg_crs[i],
                             quiet = TRUE,
                             force_redownload = TRUE)
    }, error = function(e) {print('error'); errflag <<- TRUE} )

    if(errflag) {
        errflag = FALSE
        site_data_2$reach_proportion[i] = NA
        next
    }

    site_data_2$reach_proportion[i] = rp
}

#(save progress)
saveRDS(site_data_2$reach_proportion, 'output/spatial/reach_prop_col.rds')
site_data_2$reach_proportion = readRDS('output/spatial/reach_prop_col.rds')
site_data_2 = filter(site_data_2, ! is.na(COMID))

#construct list of DSN=component pairs to acquire. see NHDPlus docs for more.
setlist = list('NHDPlusAttributes'='PlusFlowlineVAA',
               'NHDPlusAttributes'='ElevSlope')

#retrieve NHDPlusV2 data
nhdplusv2_data = nhdplusv2_bulk(site_data_2, setlist, quiet=TRUE)

#nhd variable names do not have consistent naming conventions. sometimes they're
#all caps; other times camel case. here's a crude way to deal with that.
colnames(nhdplusv2_data) = toupper(colnames(nhdplusv2_data))
nhdplusv2_data = nhdplusv2_data[, ! duplicated(colnames(nhdplusv2_data))]

#choose variables to join; filter dupes
nhdplusv2_data = select(nhdplusv2_data, COMID, STREAMORDE, FROMMEAS, TOMEAS,
                        SLOPE, REACHCODE, AREASQKM, TOTDASQKM, MAXELEVSMO,
                        MINELEVSMO)

#(save progress again)
saveRDS(nhdplusv2_data, 'output/spatial/nhdplusv2_data.rds')
nhdplusv2_data = readRDS('output/spatial/nhdplusv2_data.rds') %>%
    as_tibble() %>%
    group_by(COMID) %>%
    summarize_all(first) %>%
    ungroup()

site_data_2 = left_join(site_data_2, nhdplusv2_data, by='COMID')

#correct catchment area (AREASQKM) based on where each site falls within its reach.
#use this to correct watershed area (TOTDASQKM) and to determine an areal
#correction factor that can be multiplied with any areal summary data.
site_data_2$AREASQKM_corr = round(site_data_2$AREASQKM * site_data_2$reach_proportion, 5)
site_data_2$TOTDASQKM_corr = site_data_2$TOTDASQKM - (site_data_2$AREASQKM - site_data_2$AREASQKM_corr)
site_data_2$areal_corr_factor = site_data_2$TOTDASQKM_corr / site_data_2$TOTDASQKM

#final
# saveRDS(site_data_2, 'output/spatial/site_data2.rds')
site_data_2 = readRDS('output/spatial/site_data2.rds')

# 6: (Figure 1) GPP-ER biplot and dist plots ####

axis_cex = 2.4 #applies to labels and tick values

jpeg(width=10, height=10, units='in', res=300, quality=100, type='cairo',
     filename='figures/gpp_er_biplot_cumulAnnual.jpeg')
# par(mar=c(4.5, 4.5, 2, 2))
par(mar=c(6.5, 6.5, 2, 2))

log_gpp_fnet = log(fnet_site$GPP_site_mean)
log_er_fnet = log(fnet_site$ER_site_mean * -1) * -1
log_gpp_sp = log(sp_site$GPP_site_mean)
log_er_sp = log(sp_site$ER_site_mean * -1) * -1

gpptck = c(1, 10, 100, 1000, 10000)
# gpptck = c(0.1, 1, 10, 100, 1000, 10000)
ertck = rev(c(10000, 1000, 100, 10, 1))

# plot(log_gpp_fnet,
#      log_er_fnet, col=alpha(fnetcolor, alpha=0.5),
plot(log_gpp_fnet[fnet_high_cov_bool],
     log_er_fnet[fnet_high_cov_bool], col=alpha(fnetcolor, alpha=0.5),
     xlab='', ylab='', bg=alpha(fnetcolor, alpha=0.5),
     cex=1.5, cex.lab=axis_cex, cex.axis=axis_cex, ylim=-log(rev(range(ertck))),
     pch=21, yaxt='n', xaxt='n', xlim=log(range(gpptck)), lwd=2)
mtext(expression(paste("Cumulative GPP (gC"~"m"^"-2"*" y"^"-1"*')')),
      1, line=5, cex=axis_cex)
mtext(expression(paste("Cumulative ER (gC"~"m"^"-2"*" y"^"-1"*')')),
      2, line=3.5, cex=axis_cex)
# points(log_gpp_sp, log_er_sp, lwd=2,
#        col=alpha(spcolor, alpha=0.5), cex=1.5, pch=21, bg=alpha(spcolor, alpha=0.5))
points(log_gpp_sp[! sp_high_cov_bool], log_er_sp[! sp_high_cov_bool], lwd=2,
       col=alpha(spcolor, alpha=0.5), cex=1.5, pch=21, bg='transparent')
points(log_gpp_fnet[! fnet_high_cov_bool], log_er_fnet[! fnet_high_cov_bool],
       col=alpha(fnetcolor, alpha=0.5), cex=1.5, pch=21, bg='transparent', lwd=2)
points(log_gpp_sp[sp_high_cov_bool], log_er_sp[sp_high_cov_bool], lwd=2,
       col=alpha(spcolor, alpha=0.5), cex=1.5, pch=21, bg=alpha(spcolor, alpha=0.5))
legend('topright', legend=c('FLUXNET', 'StreamPULSE'), pch=21, bty='n', pt.cex=1.5,
       col=c(alpha(fnetcolor, alpha=0.5), alpha(spcolor, alpha=0.5)), pt.lwd=2,
       pt.bg=c(alpha(fnetcolor, alpha=0.5), alpha(spcolor, alpha=0.5)), x.intersp=2)
legend('topright', legend=c('FLUXNET', 'StreamPULSE'), pch=21, bty='n', pt.cex=1.5,
       col=c(alpha(fnetcolor, alpha=0.5), alpha(spcolor, alpha=0.5)),
       pt.bg='transparent', pt.lwd=2)
all_gpp = c(fnet_site$GPP_site_mean, sp_site$GPP_site_mean)
all_gpp[all_gpp <= 0] = NA
gpprng = range(all_gpp, na.rm=TRUE)
all_er = c(fnet_site$ER_site_mean, sp_site$ER_site_mean)
all_er[all_er >= 0] = NA
errng = range(all_er, na.rm=TRUE)

gpptck_log = log(gpptck)
axis(1, at=gpptck_log, labels=gpptck, cex.axis=axis_cex, padj=0.4)
ertck_log = log(ertck) * -1
axis(2, at=ertck_log, labels=ertck * -1, cex.axis=axis_cex, padj=0.2)

abline(a=0, b=-1, lty=2)

dev.off()

#--- distplots

jpeg(width=8, height=8, units='in', res=300, quality=100, type='cairo',
     filename='figures/gpp_er_distplots.jpeg')

par(mfrow=c(2, 1), mar=c(2,1,1,1), oma=c(0, 0, 0, 0))

#GPP
dens = density(na.omit(log_gpp_fnet))
# dens = density(na.omit(log_gpp_fnet[fnet_high_cov_bool]))
gpp_dens_fnet = tibble(x=dens$x, y=dens$y)
dens = density(na.omit(log_gpp_sp))
# dens = density(na.omit(log_gpp_sp[sp_high_cov_bool]))
gpp_dens_sp = tibble(x=dens$x, y=dens$y)

plot(gpp_dens_fnet$x, gpp_dens_fnet$y, type='n', ann=FALSE, xaxt='n', yaxt='n',
     bty='n', xlim=c(-2,9.5))
axis(1, at=gpptck_log, labels=gpptck, padj=-1.3, tick = TRUE, line=-0.2, tcl=-0.2)
mtext('GPP', 1, line=1)
polygon(x=c(gpp_dens_fnet$x, rev(gpp_dens_fnet$x)),
        y=c(gpp_dens_fnet$y, rep(0, nrow(gpp_dens_fnet))),
        col=alpha(fnetcolor, alpha=0.7),
        border=alpha(fnetcolor, alpha=0.7))
polygon(x=c(gpp_dens_sp$x, rev(gpp_dens_sp$x)),
        y=c(gpp_dens_sp$y, rep(0, nrow(gpp_dens_sp))),
        col=alpha(spcolor, alpha=0.7),
        border=alpha(spcolor, alpha=0.7))

#ER
dens = density(na.omit(log_er_fnet))
# dens = density(na.omit(log_er_fnet[fnet_high_cov_bool]))
er_dens_fnet = tibble(x=dens$x, y=dens$y)
dens = density(na.omit(log_er_sp))
# dens = density(na.omit(log_er_sp[sp_high_cov_bool]))
er_dens_sp = tibble(x=dens$x, y=dens$y)

plot(er_dens_fnet$x, er_dens_fnet$y, type='n', ann=FALSE, xaxt='n', yaxt='n',
     bty='n', xlim=c(2, -9.5))
axis(1, at=ertck_log, labels=ertck, padj=-1.3, tick = TRUE, line=-0.2, tcl=-0.2)
mtext('ER', 1, line=1)
polygon(x=c(rev(er_dens_fnet$x), er_dens_fnet$x),
        y=c(rev(er_dens_fnet$y), rep(0, nrow(er_dens_fnet))),
        col=alpha(fnetcolor, alpha=0.7),
        border=alpha(fnetcolor, alpha=0.7))
polygon(x=c(er_dens_sp$x, rev(er_dens_sp$x)),
        y=c(rev(er_dens_sp$y), rep(0, nrow(er_dens_sp))),
        col=alpha(spcolor, alpha=0.7),
        border=alpha(spcolor, alpha=0.7))

dev.off()

#NEP
nep_fnet = fnet_site$NEP_site_mean
nep_sp = sp_site$NEP_site_mean

jpeg(width=8, height=8, units='in', res=300, quality=100, type='cairo',
     filename='figures/nep_distplot.jpeg')

par(mfrow=c(2, 1), mar=c(2,1,1,1), oma=c(0, 0, 0, 0))

dens = density(na.omit(nep_fnet))
# dens = density(na.omit(nep_fnet[fnet_high_cov_bool]))
nep_dens_fnet = tibble(x=dens$x, y=dens$y)
dens = density(na.omit(nep_sp))
# dens = density(na.omit(nep_sp[sp_high_cov_bool]))
nep_dens_sp = tibble(x=dens$x, y=dens$y)

plot(nep_dens_sp$x, nep_dens_sp$y, type='n', ann=FALSE, xaxt='n', yaxt='n',
     bty='n', xlim=c(-1723, 2116))
axis(1, padj=-1.3, tick = TRUE, line=-0.2, tcl=-0.2)
mtext('NEP', 1, line=1)
polygon(x=c(rev(nep_dens_fnet$x), nep_dens_fnet$x),
        y=c(rev(nep_dens_fnet$y), rep(0, nrow(nep_dens_fnet))),
        col=alpha(fnetcolor, alpha=0.7),
        border=alpha(fnetcolor, alpha=0.7))
polygon(x=c(nep_dens_sp$x, rev(nep_dens_sp$x)),
        y=c(rev(nep_dens_sp$y), rep(0, nrow(nep_dens_sp))),
        col=alpha(spcolor, alpha=0.7),
        border=alpha(spcolor, alpha=0.7))

dev.off()

# 7: (Figures 2 and 3) lips plots and probability densities ####

dir.create('figures/lips', showWarnings = FALSE)
dir.create('figures/probdens', showWarnings = FALSE)

axis_cex = 2.4 #applies to labels and tick values

lips_plot = function(quant_filt=NULL, outfile, use_fnet=FALSE, ylims=NULL,
                     wee=FALSE){

    dset = if(use_fnet) fnet_lips else sp_lips
    dset_site = if(use_fnet) fnet_site else sp_site

    var_quant_filt = NULL
    if(! is.null(quant_filt)){
        quant_comp = strsplit(quant_filt, ' ')[[1]]
        qf = quantile(dset_site[, quant_comp[1], drop=TRUE], na.rm=TRUE,
        probs=as.numeric(quant_comp[3]))
        filt_sites = dset_site %>%
            filter_(paste(quant_comp[1], quant_comp[2], qf)) %>%
            pull(sitecode)
        dset = filter(dset, sitecode %in% filt_sites)

        var_quant_filt = paste0(quant_comp[1], ' ', quant_comp[2], ' ',
                                as.numeric(quant_comp[3]) * 100, '%')
    }

    nsites_included = length(unique(dset$sitecode))

    smry = dset %>%
        select(DOY, GPP_C_filled, ER_C_filled, NEP_C_filled) %>%
        group_by(DOY) %>%
        summarize_all(list(median=~median(., na.rm=TRUE),
                           quant25=~quantile(., na.rm=TRUE)[2],
                           quant75=~quantile(., na.rm=TRUE)[4])) %>%
        ungroup()

    jpeg(width=11, height=11, units='in', filename=outfile, type='cairo',
         res=300, quality=100)
    # pdf(file=outfile, width=10, height=10)

    if(! is.null(ylims)){

        # plot_dims = top_dims = dev.size('in')
        # top_dims[2] = plot_dims[2] * abs(ylims[2]) / sum(abs(ylims))

        gpplim = c(0, ylims[2])
        erlim = c(ylims[1], 0)

    } else {

        # top_dims = dev.size('in')

        gpplim=c(0, 5)
        erlim=c(-5, 0)

        if(use_fnet){
            gpplim[2] = 11
            erlim[1] = -7
        }
    }

    if(wee){
        axis_cex_ = axis_cex * 2
        line_ = -1.5
        par(mfrow=c(2, 1), oma=c(1, 1, 0, 0), mar=c(0, 5, 3, 1), lend=2)
    } else {
        axis_cex_ = axis_cex
        line_ = -2
        par(mfrow=c(2, 1), oma=c(1, 1, 0, 0), mar=c(0, 5, 3, 1), lend=2)
    }
    # fin=top_dims, new=TRUE)

    plot(smry$DOY, smry$GPP_C_filled_median, ylab='', yaxs='i', type='l',
         bty='n', lwd=4, xlab='', ylim=gpplim, xaxs='i', xaxt='n', yaxt='n',
         col='gray30')
    polygon(x=c(smry$DOY, rev(smry$DOY)),
            y=c(smry$GPP_C_filled_quant25, rev(smry$GPP_C_filled_quant75)),
            border=NA, col=alpha('forestgreen', alpha=0.6))
    axislocs = if(max(ylims) >= 10) seq(0, 10, 2) else c(0, 2.5, 5)
    axis(2, las=0, line=0, xpd=NA, tck=-.02, labels=FALSE,
         at=axislocs, cex.axis=axis_cex_, tcl=-0.3)
    # at=round(seq(0, gpplim[2], length.out=5), 1))
    axis(2, las=0, line=-0.5, xpd=NA, tcl=0, col='transparent',
         at=axislocs, cex.axis=axis_cex_)
    # at=round(seq(0, gpplim[2], length.out=5), 1))
    abline(h=0, lty=2, lwd=2, col='gray60')
    medsums = round(colSums(select(smry, contains('median'))), 1)

    if(! wee){
        mtext(expression(paste(bold("gC") ~ bold("m") ^ bold("-2") ~
                                   bold(" d") ^ bold("-1"))), side=2,
              line=line_, outer=TRUE, cex=axis_cex_)
    }

    # if(filter_label){
    #     legend('topright', title='Filters', bty='n', title.col='gray30',
    #            lty=1, seg.len=0.2, lwd=2, legend=c(..., var_quant_filt))
    # }

    legend('right', title='Cumulative\nMedian Sums', bty='n',
           legend=c(paste('GPP:', medsums[1]), paste('ER:', medsums[2]),
                    paste('NEP:', medsums[3])), title.col='gray30')
    legend('left', paste('Sites included:', nsites_included), bty='n')

    # if(! is.null(ylims)){
    #     plot_dims = bottom_dims = dev.size('in')
    #     bottom_dims[2] = plot_dims[2] * abs(ylims[1]) / sum(abs(ylims))
    # } else {
    #     bottom_dims = dev.size('in')
    # }

    par(mar=c(4, 5, 0, 1))#, fin=bottom_dims)
    if(wee){
        padj_ = 0.8
        DOY1 = ''
    } else {
        padj_ = 0.5
        DOY1 = 1
    }

    plot(smry$DOY, smry$ER_C_filled_median, ylab='', yaxs='i', type='l',
         bty='n', lwd=4, xlab='', ylim=erlim, xaxs='i', xaxt='n', yaxt='n')
    polygon(x=c(smry$DOY, rev(smry$DOY)),
            y=c(smry$ER_C_filled_quant25, rev(smry$ER_C_filled_quant75)),
            border=NA, col=alpha('sienna', alpha=0.6))
    if(max(ylims) >= 10){
        axislocs = seq(0, -10, -2)
        axisextra = -11
    } else {
        axislocs = c(0, -2.5, -5)
        axisextra = NULL
    }
    axis(2, las=0, line=0, xpd=NA, tck=-.02, labels=FALSE,
         # at=round(seq(0, erlim[1], length.out=5), 1))
         at=c(axislocs, axisextra), cex.axis=axis_cex_, tcl=-0.3)
    axis(2, las=0, line=-0.5, xpd=NA, tcl=0, col='transparent',
         # at=round(seq(0, erlim[1], length.out=5), 1))
         at=axislocs, cex.axis=axis_cex_)
    axis(1, line=0, tck=-.02, labels=FALSE, at=c(1, seq(60, max(smry$DOY), 60)),
         cex.axis=axis_cex_, tcl=-0.3)
    axis(1, line=-0.5, tcl=0, col='transparent', at=c(1, seq(60, max(smry$DOY), 60)),
         cex.axis=axis_cex_, padj = padj_, tcl=-0.3,
         labels = c(DOY1, seq(60, max(smry$DOY), 60)))
    lines(smry$DOY, smry$NEP_C_filled_median, col='black', lwd=4, xpd=NA, lend=1)

    if(! wee){
        mtext('DOY', side=1, line=3.5, font=2, cex=axis_cex_)
    }
    dev.off()
}

pdf_plot = function(var, outfile){

    jpeg(width=7, height=4, units='in', filename=outfile, type='cairo',
         res=300, quality=100)
    # pdf(outfile, height=4, width=7)

    vv = na.omit(sort(sp_site[[var]]))
    dens = density(vv)
    vq = quantile(vv, probs=c(0.25, 0.75))
    densdf = tibble(x=dens$x, y=dens$y)
    dens25 = dens75 = densdf
    dens25 = dens25[densdf$x <= vq[1], ]
    dens75 = dens75[densdf$x >= vq[2], ]
    plot(densdf$x, densdf$y, type='l', xlab='width', ylab='density', bty='l',
         col='gray50', lwd=2)
    polygon(x=c(dens25$x, rev(dens25$x)),
            y=c(dens25$y, rep(0, nrow(dens25))), col='gray50', border='gray50')
    polygon(x=c(dens75$x, rev(dens75$x)),
            y=c(dens75$y, rep(0, nrow(dens75))), col='gray50', border='gray50')

    dev.off()
}

#overall
lips_ylim = c(-11, 11)
lips_plot(quant_filt=NULL, outfile='figures/lips/lips_overall_sp.jpeg',
          ylims=lips_ylim)
lips_plot(quant_filt=NULL, outfile='figures/lips/lips_overall_fnet.jpeg',
          use_fnet=TRUE, ylims=lips_ylim)

pdf_plot('Disch_ar1', 'figures/probdens/probdens_Qar1.jpeg')
pdf_plot('MOD_ann_NPP', 'figures/probdens/probdens_MODIS.jpeg')
pdf_plot('Stream_PAR_sum', 'figures/probdens/probdens_PAR.jpeg')

#subsets
lips_ylim = c(-5, 5)
lips_plot(quant_filt='Disch_ar1 > 0.75', outfile='figures/lips/lips_Qar1_75.jpeg',
          ylims=lips_ylim, wee=TRUE)
lips_plot(quant_filt='Disch_ar1 < 0.25', outfile='figures/lips/lips_Qar1_25.jpeg',
          ylims=lips_ylim, wee=TRUE)
lips_plot(quant_filt='MOD_ann_NPP > 0.75', outfile='figures/lips/lips_MODIS_75.jpeg',
          ylims=lips_ylim, wee=TRUE)
lips_plot(quant_filt='MOD_ann_NPP < 0.25', outfile='figures/lips/lips_MODIS_25.jpeg',
          ylims=lips_ylim, wee=TRUE)
lips_plot(quant_filt='Stream_PAR_sum > 0.75', outfile='figures/lips/lips_PAR_75.jpeg',
          ylims=lips_ylim, wee=TRUE)
lips_plot(quant_filt='Stream_PAR_sum < 0.25', outfile='figures/lips/lips_PAR_25.jpeg',
          ylims=lips_ylim, wee=TRUE)

# OBSOLETE?: df for emily to build annual rates dist plot ####

terr_aq_cumul_metab = sp_site %>%
    mutate(source = 'streampulse') %>%
    bind_rows(fnet_site) %>%
    select(sitecode, source, GPP_site_mean, ER_site_mean) %>%
    rename(sum_gpp_C=GPP_site_mean, sum_er_C=ER_site_mean) %>%
    mutate(
        source = ifelse(is.na(source), 'fluxnet', source),
        sum_gpp_C = round(sum_gpp_C, 2),
        sum_er_C = round(sum_er_C, 2)) %>%
    arrange(source, sitecode)

write.csv(terr_aq_cumul_metab,
          'export_datasets/fnet_sp_cumul_metab.csv',
          row.names=FALSE)

# 9: (Figure 4) bubble plots ####

dir.create('figures/bubble_plots', showWarnings = FALSE)

bubble_plot = function(xvar, comp, logx=FALSE, outfile){

    jpeg(width=8, height=8, units='in', filename=outfile, type='cairo',
         res=300, quality=100)
    # pdf(outfile, height=8, width=8)
    par(mar=c(5, 5, 4, 6))

    plotcol = case_when(comp == 'GPP_site_mean' ~ gppcolor,
                        comp == 'ER_site_mean' ~ ercolor,
                        comp == 'NEP_site_mean' ~ 'black')

    if(xvar == 'Stream_PAR_sum'){
        xxlab = 'Light Availability (Mean Annual Surface PAR)'
    } else {
        xxlab = 'MODIS NPP'
    }

    gpprng = range(sp_site$GPP_site_mean, na.rm=TRUE)
    if(comp == 'ER_site_mean'){
        xx = abs(sp_site[[comp]])
    } else {
        xx = sp_site[[comp]]
    }
    if(comp != 'NEP_site_mean'){
        rescaled = ((xx - gpprng[1]) /
                        (gpprng[2] - gpprng[1])) * (5 - 1) + 1
    } else {
        xxrng = range(xx, na.rm=TRUE)
        # rescaled = ((xx - xxrng[1]) /
        #                 (xxrng[2] - xxrng[1])) * (3) + 0.2
        rescaled = ((xx - -2000) /
                        (1000 - -2000)) * (4) + 0.5
    }

    xxvar = sp_site[[xvar]]
    if(logx){
        xxvar = log(xxvar)
        if(xvar == 'Stream_PAR_sum'){
            xlm = c(.5, 2.8)
        } else {
            xlm = c(4.3, 7.2)
        }
    } else {
        if(xvar == 'Stream_PAR_sum'){
            xlm = c(2, 15.5)
        } else {
            xlm = range(xxvar, na.rm=TRUE)
        }
    }

    plot(xxvar,
         sp_site$Disch_ar1, pch=21,
    # plot(xxvar[sp_high_cov_bool],
    #      sp_site$Disch_ar1[sp_high_cov_bool], pch=21,
         xlab=xxlab, xaxt='n',
         ylab='Predictability of Flow (Discharge AR-1 Coeff.)',
         col=alpha(plotcol, alpha=0.5), bty='o',
         ylim=c(0.1, 1), xlim=xlm,
         xpd=NA, main='', cex=rescaled, font.lab=2,
         # xpd=NA, main='', cex=rescaled[sp_high_cov_bool], font.lab=2,
         bg=alpha(plotcol, alpha=0.5))
    # points(xxvar[! sp_high_cov_bool],
    #        sp_site$Disch_ar1[! sp_high_cov_bool], pch=21,
    #        col=alpha(plotcol, alpha=0.5), bg='transparent',
    #        cex=rescaled[! sp_high_cov_bool])

    if(comp != 'NEP_site_mean'){
        legend('right', legend=c(expression(paste(0.01)), '', '',
                                 expression(paste(4000))),
               bty='n', pch=21, pt.bg='transparent', x.intersp=1.7,
               pt.cex=c(1, 2, 3, 5), col='gray30', xpd=NA,
               y.intersp=c(1, 2, 1.2, 1.6),
               inset=c(-0.16, 0), title='')
    } else {
        legend('right', legend=c(expression(paste(-2000)), '', '',
                                 expression(paste(1000))),
               bty='n', pch=21, pt.bg='transparent', x.intersp=1.7,
               pt.cex=seq(0.5, 4, length.out=4), col='gray30', xpd=NA,
               y.intersp=c(1, 2, 1.2, 1.6),
               inset=c(-0.18, 0), title='')
    }
    if(comp == 'GPP_site_mean'){
        legend('right', legend=c('','','',''),
               bty='n', pch=21, pt.bg='transparent', x.intersp=1.5,
               pt.cex=c(1, 2, 3, 5), col='transparent', xpd=NA,
               inset=c(-0.15, 0),
               title=expression(paste(bold('Cumul.\nAnnual\nGPP'))))
    } else if(comp == 'ER_site_mean'){
        legend('right', legend=c('','','',''),
               bty='n', pch=21, pt.bg='transparent', x.intersp=1.5,
               pt.cex=c(1, 2, 3, 5), col='transparent', xpd=NA,
               inset=c(-0.15, 0),
               title=expression(paste(bold('Cumul.\nAnnual\nER'))))
    } else {
        legend('right', legend=c('','','',''),
               bty='n', pch=21, pt.bg='transparent', x.intersp=1.5,
               pt.cex=c(1, 2, 3, 5), col='transparent', xpd=NA,
               inset=c(-0.15, 0),
               title=expression(paste(bold('Cumul.\nAnnual\nNEP'))))
    }

    if(logx){
        if(xvar == 'Stream_PAR_sum'){
            tcks = c(1, 2, 4, 8, 16)
        } else {
            tcks = c(200, 400, 800, 1600, 3200)
        }
        tcks_log = log(tcks)
        axis(1, at=tcks_log, labels=tcks)
    } else {
        axis(1)
    }

    dev.off()
}

bubble_plot(xvar='Stream_PAR_sum', comp='GPP_site_mean', logx=FALSE, outfile='figures/bubble_plots/PAR_GPP_linear.jpeg')
bubble_plot(xvar='Stream_PAR_sum', comp='GPP_site_mean', logx=TRUE, outfile='figures/bubble_plots/PAR_GPP_log.jpeg')
bubble_plot(xvar='Stream_PAR_sum', comp='ER_site_mean', logx=FALSE, outfile='figures/bubble_plots/PAR_ER_linear.jpeg')
bubble_plot(xvar='Stream_PAR_sum', comp='ER_site_mean', logx=TRUE, outfile='figures/bubble_plots/PAR_ER_log.jpeg')
bubble_plot(xvar='Stream_PAR_sum', comp='NEP_site_mean', logx=FALSE, outfile='figures/bubble_plots/PAR_NEP_linear.jpeg')
bubble_plot(xvar='Stream_PAR_sum', comp='NEP_site_mean', logx=TRUE, outfile='figures/bubble_plots/PAR_NEP_log.jpeg')
bubble_plot(xvar='MOD_ann_NPP', comp='GPP_site_mean', logx=FALSE, outfile='figures/bubble_plots/MODNPP_GPP_linear.jpeg')
bubble_plot(xvar='MOD_ann_NPP', comp='GPP_site_mean', logx=TRUE, outfile='figures/bubble_plots/MODNPP_GPP_log.jpeg')
bubble_plot(xvar='MOD_ann_NPP', comp='ER_site_mean', logx=FALSE, outfile='figures/bubble_plots/MODNPP_ER_linear.jpeg')
bubble_plot(xvar='MOD_ann_NPP', comp='ER_site_mean', logx=TRUE, outfile='figures/bubble_plots/MODNPP_ER_log.jpeg')
bubble_plot(xvar='MOD_ann_NPP', comp='NEP_site_mean', logx=FALSE, outfile='figures/bubble_plots/MODNPP_NEP_linear.jpeg')
bubble_plot(xvar='MOD_ann_NPP', comp='NEP_site_mean', logx=TRUE, outfile='figures/bubble_plots/MODNPP_NEP_log.jpeg')

# 10: export regression dataset (needs update) ####

coverage_tb = sp %>%
    group_by(sitecode) %>%
    summarize(
        nyear = first(nyear),
        coverage = first(coverage)) %>%
    ungroup()

stats_set = site_data_2 %>%
    full_join(sp_site, by='sitecode') %>%
    left_join(coverage_tb) %>%
    select(sitecode, GPP_site_mean, ER_site_mean,
           nyear, coverage, Disch_ar1, MOD_ann_NPP, Stream_PAR_sum,
           lat, lon, epsg_crs, area_km_corr = TOTDASQKM_corr) %>%
    filter(
        ! is.na(GPP_site_mean),
        ! is.na(ER_site_mean)) %>%
    arrange(sitecode)

write.csv(stats_set, 'export_datasets/streampulse_synthesis_statset.csv',
          row.names=FALSE)

# OBSOLETE: data for emily to explore ####

#accumulate all site data
width = readRDS('~/git/streampulse/metab_synthesis/data/lotic_streamlight_params.rds') %>%
    as_tibble() %>%
    select(sitecode = Site_ID, width = Width)
site_data_A = readRDS('output/synthesis_site_metrics.rds') %>%
    as_tibble() %>%
    rename(sitecode = Site_ID) %>%
    full_join(width, by = 'sitecode')
site_data_B = readRDS('output/spatial/site_data2.rds') %>%
    select(sitecode, lat, lon, stream_order = STREAMORDE, slope = SLOPE,
           ws_area = TOTDASQKM_corr)
site_data = full_join(site_data_A, site_data_B) %>%
    rename_all(function(x) paste0('__', x)) %>%
    rename(sitecode = '__sitecode')

#not really "unmodified". ER_K600_R2 and max_K600 have already been added
unmodified_from_phil = Reduce(bind_rows, sp_list) %>%
    as_tibble()

#missing metab rows removed, ERxK R2 < 0.5, max K < 100, impossible GPP & ER corrected
cleaned_up = unmodified_from_phil %>%
    filter( #remove records with missing metab and/or spurious model results
        ! is.na(GPP),
        ! is.na(ER),
        ER_K600_R2 < 0.5,#) %>%
        max_K600 < 100) %>%
    mutate(GPP = case_when(GPP < 0 ~ 0, #correct impossible metab vals
                           TRUE ~ GPP),
           ER = case_when(ER > 0 ~ 0,
                          TRUE ~ ER)) %>%
    select(-starts_with('U_ID'))

#summary 1 (by site-DOY)
summarized_by_siteDOY = cleaned_up %>%
    group_by(sitecode, DOY) %>%
    summarize(across(everything(), mean, na.rm = TRUE),
              nyear = length(unique(Year))) %>%
    ungroup() %>%
    group_by(sitecode) %>%
    mutate(coverage = n() / 365) %>% #get proportion nonmissing for the aggregate year
    arrange(DOY) %>%
    mutate( #fill gaps, convert from gO2/m^2/d to gC/m^2/d
        GPP_C_filled = O2_to_C(impute_ts(GPP)),
        ER_C_filled = O2_to_C(impute_ts(ER))) %>%
    ungroup() %>%
    select(-GPP, -ER) %>%
    arrange(sitecode, DOY) %>%
    mutate(NEP_C_filled = GPP_C_filled + ER_C_filled) %>%
    left_join(site_data, by = 'sitecode') #include site data

#summary 2 (by site-year)
summarized_by_siteyear = cleaned_up %>%
    group_by(sitecode, Year) %>%
    arrange(DOY) %>%
    mutate(
        GPP = impute_ts(GPP),
        ER = impute_ts(ER)) %>%
    summarize(across(all_of(c('GPP', 'ER', 'K600', 'DO_obs', 'DO_sat',
                              'temp_water', 'LAI_proc', 'discharge',
                              'PAR_sum', 'Stream_PAR_sum',
                              'ER_K600_R2', 'max_K600')),
                     .fn = list(mean = mean),
                     .names = '{col}_{fn}'),
              across(all_of(c('GPP', 'ER', 'LAI_proc', 'discharge', 'PAR_sum',
                              'Stream_PAR_sum')),
                     .fn = list(sum = sum),
                     .names = '{col}_{fn}')) %>%
    mutate(coverage = n() / 365) %>% #get proportion nonmissing for the aggregate year
    ungroup() %>%
    rename(PAR_sum = PAR_sum_sum, Stream_PAR_sum = Stream_PAR_sum_sum,
           PAR_mean = PAR_sum_mean, Stream_PAR_mean = Stream_PAR_sum_mean) %>%
    group_by(sitecode) %>%
    mutate( #fill gaps, convert from gO2/m^2/d to gC/m^2/d
        GPP_C_sum_filled = O2_to_C(GPP_sum),
        ER_C_sum_filled = O2_to_C(ER_sum)) %>%
    ungroup() %>%
    select(-GPP_sum, -ER_sum) %>%
    arrange(sitecode, Year) %>%
    mutate(NEP_C_sum_filled = GPP_C_sum_filled + ER_C_sum_filled) %>%
    left_join(site_data, by = 'sitecode') #include site data

# #summary 3 (by site)
# summarized_by_site = sp_site
summarized_by_site = summarized_by_siteDOY %>%
    group_by(sitecode) %>%
    summarize(across(all_of(c('GPP_C_filled', 'ER_C_filled', 'NEP_C_filled',
                              'K600', 'DO_obs', 'DO_sat',
                              'temp_water', 'LAI_proc', 'discharge',
                              'PAR_sum', 'Stream_PAR_sum')),
                     .fn = list(mean = mean),
                     .names = '{col}_{fn}'),
              across(all_of(c('GPP_C_filled', 'ER_C_filled', 'NEP_C_filled',
                              'LAI_proc', 'discharge', 'PAR_sum',
                              'Stream_PAR_sum')),
                     .fn = list(sum = sum),
                     .names = '{col}_{fn}'),
              across(all_of(starts_with(c('__', 'coverage', 'nyear'))),
                     # across(all_of(c('lat', 'lon', 'stream_order', 'slope', 'ws_area',
                     #                 'coverage', 'nyear')),
                     .fn = first,
                     .names = '{col}')) %>%
    ungroup() %>%
    rename(PAR_sum = PAR_sum_sum, Stream_PAR_sum = Stream_PAR_sum_sum,
           PAR_mean = PAR_sum_mean, Stream_PAR_mean = Stream_PAR_sum_mean) %>%
    arrange(sitecode)

#write outputs
write_csv(summarized_by_siteDOY,
          'export_datasets/sp_summarized_by_site-DOY.csv')
write_csv(summarized_by_siteyear,
          'export_datasets/sp_summarized_by_site-year.csv')
write_csv(summarized_by_site,
          'export_datasets/sp_summarized_by_site.csv')

# junk? ####

# fnet_site = fnet %>%
#     group_by(sitecode) %>%
#     summarize(GPP_ann_mean = mean(GPP_C_filled, na.rm=TRUE),
#               ER_ann_mean = mean(ER_C_filled, na.rm=TRUE),
#               GPP_site_mean = sum(GPP_C_filled, na.rm=TRUE),
#               ER_site_mean = sum(ER_C_filled, na.rm=TRUE),
#               nyear = first(nyear)) %>%
#     ungroup() %>%
#     mutate(
#         NEP_ann_mean = GPP_ann_mean + ER_ann_mean,
#         NEP_site_mean = GPP_site_mean + ER_site_mean) %>%
#     arrange(sitecode)

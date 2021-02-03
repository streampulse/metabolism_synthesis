comid_from_point = function(lat, long, CRS) {
    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    COMID = nhdplusTools::discover_nhdplus_id(ptc)
    if(! length(COMID)) COMID = NA
    return(COMID)
}

vpu_from_point = function(lat, long, CRS) {
    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    VPU = nhdR::find_vpu(ptc)
    return(VPU)
}

#this calculates how far along a reach any given point falls. That way when we pull in
#watershed summary data for a reach, we can adjust it according to how much
#of the total upstream area actually contributes to the point in question.
# A value of 0 means upstream end; 1 means downstream end.
calc_reach_prop = function(VPU, COMID, lat, long, CRS, quiet=FALSE,
                           force_redownload = FALSE){

    if(! quiet){
        message(paste0('The nhdR package downloads NHDPlusV2 components to ',
            nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
            ' Fortunately, each component need only be downloaded once.'))
    }

    fl = nhdR::nhd_plus_load(vpu=VPU, component='NHDSnapshot',
        dsn='NHDFlowline', approve_all_dl=TRUE, force_dl = force_redownload)
    fl_etc = nhdR::nhd_plus_load(vpu=VPU, component='NHDPlusAttributes',
        dsn='PlusFlowlineVAA', approve_all_dl=TRUE, force_dl = force_redownload)

    colnames(fl)[colnames(fl) == 'ComID'] = 'COMID'
    colnames(fl)[colnames(fl) == 'ReachCode'] = 'REACHCODE'
    fl = fl[fl$COMID == COMID,]
    fl = left_join(fl, fl_etc[, c('ComID', 'ToMeas', 'FromMeas')],
        by=c('COMID'='ComID'))

    pt = sf::st_point(c(long, lat))
    ptc = sf::st_sfc(pt, crs=CRS)
    ptct = sf::st_transform(ptc, crs=4269) #CRS for NAD 83
    x = suppressWarnings(nhdplusTools::get_flowline_index(fl, points=ptct))
    out = 1 - x$REACH_meas / 100 #0=upstream end; 1=downstream end

    return(out)
}

#this acquires nhdplusv2 data for a single site by COMID.
#it's just a thin wrapper around nhdR::nhd_plus_load
nhdplusv2_from_comid = function(VPU, COMID, component, DSN, quiet=FALSE) {

    if(! quiet){
        message(paste0('The nhdR package downloads NHDPlusV2 components to ',
            nhdR:::nhd_path(), '. Unfortunately this cannot be changed.',
            ' Fortunately, each component need only be downloaded once.'))
    }

    data = nhdR::nhd_plus_load(vpu=VPU, component=component, dsn=DSN,
        approve_all_dl=TRUE)

    colnames(data)[colnames(data) == 'ComID'] = 'COMID'
    colnames(data)[colnames(data) == 'ReachCode'] = 'REACHCODE'
    data = data[data$COMID == COMID,]

    return(data)
}

#this calls nhdplusv2_from_comid repeatedly to get data for all your sites.
#the dataframe must include COMID and VPU columns
nhdplusv2_bulk = function(site_df, nhdplusv2_sets, quiet=FALSE){

    nhdplus_data = data.frame()
    if(any(is.na(site_df$COMID))) stop('you have missing COMIDs')

    for(j in 1:nrow(site_df)){
        for(i in 1:length(setlist)){
            print(paste(j, nhdplusv2_sets[[i]]))

            if(i == 1 || initerr){
                row_base = try(nhdplusv2_from_comid(
                    VPU = site_df$VPU[j],
                    COMID = site_df$COMID[j],
                    component = names(setlist[i]),
                    DSN = setlist[[i]],
                    quiet = quiet))
                if('try-error' %in% class(row_base) || nrow(row_base) > 1){
                    initerr = TRUE
                    row_base = data.frame(COMID=site_df$COMID[j])
                } else {
                    initerr = FALSE
                }
            } else {
                row_ext = try(nhdplusv2_from_comid(site_df$VPU[j],
                    site_df$COMID[j], names(setlist[i]), setlist[[i]],
                    quiet=quiet))
                if(! 'try-error' %in% class(row_ext) && nrow(row_ext) == 1){
                    row_base = left_join(row_base, row_ext)
                }
            }

        }

        if(nrow(row_base) > 1){
            row_base = data.frame(COMID=site_df$COMID[j])
        }

        nhdplus_data = rbind.fill(nhdplus_data, row_base)
    }

    return(nhdplus_data)
}

query_streamcat_datasets = function(keyword=NULL){

    ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
        'NHDPlusLandscapeAttributes/StreamCat/States/')

    url_list = getURL(ftpdir, dirlistonly=TRUE)
    url_list = strsplit(url_list, split='\n')[[1]]

    if(! is.null(keyword)){
        url_list = url_list[grep(keyword, url_list, ignore.case=TRUE)]
    }

    return(url_list)
}

#this function acquires streamcat data for a single site by NHDPlusV2 COMID.
streamcat_from_comid = function(USstate, COMID, dataset){

    ftpdir = paste0('ftp://newftp.epa.gov/EPADataCommons/ORD/',
        'NHDPlusLandscapeAttributes/StreamCat/States/')
    zip_name = paste0(dataset, '_', USstate, '.zip')

    csv_name = gsub('.zip', '.csv', zip_name)
    temp = tempfile()
    download.file(paste0(ftpdir, zip_name), temp)
    data = read.csv(unz(temp, csv_name), stringsAsFactors=FALSE)
    data = data[data$COMID == COMID,]

    return(data)
}

#this calls streamcat_from_comid repeatedly to get data for all your sites
#the dataframe must include COMID and region columns, where "region" refers to
#each state's 2-letter abbreviation.
streamcat_bulk = function(site_df, streamcat_sets){

    streamcat_data = data.frame()
    if(any(is.na(site_df$COMID))) stop('you have missing COMIDs')

    for(j in 1:nrow(site_df)){
        for(i in 1:length(streamcat_sets)){
            print(paste(j, streamcat_sets[i]))

            if(i == 1 || initerr){
                row_base = try(streamcat_from_comid(site_df$region[j],
                    site_df$COMID[j], streamcat_sets[i]))
                if('try-error' %in% class(row_base) || nrow(row_base) > 1){
                    initerr = TRUE
                    row_base = data.frame(COMID=site_df$COMID[j])
                } else {
                    initerr = FALSE
                }
            } else {
                row_ext = try(streamcat_from_comid(site_df$region[j],
                    site_df$COMID[j], streamcat_sets[i]))
                if(! 'try-error' %in% class(row_ext) && nrow(row_ext) == 1){
                    row_base = left_join(row_base, row_ext)
                }
            }

        }

        if(nrow(row_base) > 1){
            row_base = data.frame(COMID=site_df$COMID[j])
        }

        streamcat_data = rbind.fill(streamcat_data, row_base)
    }

    return(streamcat_data)
}

filter_and_impute = function(diagnostics, models, terr=FALSE, ...){
    # diagnostics=fake_diag; models=fnet_list

    filt = filter_metab(diag=diagnostics, 'ER_K <= 1', metab_rds=models)
    imp = synthesis_gapfill(filt, PQ=1.25, block=Inf, pmiss=99, terr)

    return(imp)
}

consolidate_list = function(daily_summaries){

    daily_summaries = Map(function(x){
        if(is.null(x)) return()
        x = select(x, DOY, GPP_C_filled, ER_C_filled, NEP_C_filled)
        return(x)
    }, daily_summaries)

    smry = Reduce(function(x, y){
        out = bind_rows(x, y)
        return(out)
    }, daily_summaries)

    return(smry)
}

consolidate_list_withnames = function(daily_summaries){

    daily_summaries = Map(function(x){
        if(is.null(x)) return()
        x = select(x, sitename, DOY, GPP_C_filled, ER_C_filled, NEP_C_filled)
        return(x)
    }, daily_summaries)

    smry = Reduce(function(x, y){
        out = bind_rows(x, y)
        return(out)
    }, daily_summaries)

    return(smry)
}

phil_to_mike_format = function(dset, mod_dset, arrange=TRUE){

    nwis_ind = substr(dset$Site_ID, 1, 4) == 'nwis'
    dset$Site_ID[nwis_ind] = gsub('_', '-', dset$Site_ID[nwis_ind])
    dset$site = unname(sapply(dset$Site_ID, function(x) strsplit(x, '_')[[1]][2]))
    dset$site[is.na(dset$site)] = dset$Site_ID[is.na(dset$site)]
    dset$Site_ID = NULL

    #join region column from model outputs (which for streampulse == US state)
    dset = dset %>%
        left_join(mod_dset, 'site') %>%
        mutate(sitecode=paste(region, site, sep='_')) %>%
        filter(! is.na(sitecode)) %>%
        select(-year) %>%
        distinct()

    if(arrange){
        dset = arrange(dset, sitecode)
    }

    return(dset)
}

gO2_to_gC = function(x, PQ){
    x = x * (1 / (2 * 15.9994)) * (1 / PQ) * 12.0107
    return(x)
}

temp20_std = function(resp, temp){
    #resp = R20 * 1.047 ^ (temp - 20) #solving for R20 here
    logterm = log(1.047)
    R20 = exp(log(resp) / (logterm * temp - logterm * 20))
    return(R20)
}

is_erk_correlated = function(er, k){

    erk_p = try(summary(lm(er ~ k))$coefficients[2,4])
    if('try-error' %in% class(erk_p)) erk_p = NA

    is_correlated = erk_p <= 0.5

    return(is_correlated)
}

erk_r2 = function(er, k){

    r2 = try(summary(lm(er ~ k))$adj.r.squared)
    if('try-error' %in% class(r2)) r2 = NA

    return(r2)
}

impute_ts = function(x){

    na_loc = is.na(x)
    if(all(! na_loc)) return(x)

    x = window(x)

    #first value can't be NA, so fill it
    if(is.na(x)[1]){
        x[1] = zoo::na.locf(x,
                            option = "nocb",
                            na.remaining = "rev")[1]
    }

    #fit structural time series if possible, linear interp if not
    # tryCatch({
    tryCatch({
        x_interp = StructTS(x, type = 'trend')
    }, error=function(e){
        x_interp = StructTS(x, type = "level")
    }, warning=function(w){
        x_interp = StructTS(x, type = "level")
    })

    # }, error = function(e){
    #     NULL
    # })

    #smooth the imputed series; replace missing values in the original series
    x_smooth = tsSmooth(x_interp)
    x[na_loc] = x_smooth[na_loc, 1]

    return(x)
}

O2_to_C = function(x_O2){

    x_C = x_O2 * (1 / (2 * 15.9994)) * (1 / 1.25) * 12.0107

    return(x_C)
}

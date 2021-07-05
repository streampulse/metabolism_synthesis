
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

    try_result = try(nhdR::nhd_plus_get(vpu = VPU,
                                        component = "NHDSnapshot",
                                        force_unzip = TRUE))

    if(inherits(try_result, 'try-error') || try_result != 0){
        nhdR::nhd_plus_get(vpu = VPU, component = "NHDSnapshot", force_dl = TRUE)
        nhdR::nhd_plus_get(vpu = VPU, component = "NHDSnapshot", force_unzip = TRUE)
    }

    try_result = try(nhdR::nhd_plus_get(vpu = VPU,
                                        component = "NHDPlusAttributes",
                                        force_unzip = TRUE))

    if(inherits(try_result, 'try-error') || try_result != 0){
        nhdR::nhd_plus_get(vpu = VPU, component = "NHDPlusAttributes", force_dl = TRUE)
        nhdR::nhd_plus_get(vpu = VPU, component = "NHDPlusAttributes", force_unzip = TRUE)
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


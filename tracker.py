#!/usr/bin/python
## Tsung-Lin Hsieh (hsiehtl@princeton.edu)

import numpy as np
import xarray as xr
import scipy.ndimage.measurements as measure
from scipy.interpolate import interp2d
import itertools
import sys

def tracker(rain, vort, vref, rain_percentile, dist_threshold, size_threshold, latlim, interp50=True, return_snapshot=False):
    if interp50:
        ds50 = xr.open_dataset('HiRAM_land_static.nc') # 50km model grid
        ds50 = ds50.rename({'grid_yt':'lat', 'grid_xt':'lon'})
        
    lat_dim = rain.dims[1]
    lon_dim = rain.dims[2]
    rain = rain.rename({lat_dim:'lat', lon_dim:'lon'})
    vort = vort.rename({lat_dim:'lat', lon_dim:'lon'})
    vref = vref.rename({lat_dim:'lat', lon_dim:'lon'})

    latbound = (-latlim-5, latlim+5) # values of lat, not the indices
    ilatbound = (np.searchsorted(vort.lat, latbound[0]), np.searchsorted(vort.lat, latbound[1]))
    if interp50:
        ilatbound_lores = (np.searchsorted(ds50.lat, latbound[0]), np.searchsorted(ds50.lat, latbound[1]))
        lon, lat = np.meshgrid(ds50.lon, ds50.lat[ilatbound_lores[0]:ilatbound_lores[1]])
    else: # 50km models
        lon, lat = np.meshgrid(vort.lon, vort.lat[ilatbound[0]:ilatbound[1]])

    tracking = []
    tmp = []
    completed = []

    itrange = range(len(vort.time))
    # itrange = range(100*4, 300*4) # fast

    for it in itrange:
        currtime = vort.indexes['time'][it]

        vort_trop = vort[it,ilatbound[0]:ilatbound[1],:]*np.sign(vort.lat[ilatbound[0]:ilatbound[1]]) # multiplied by sign(lat) so that maximum() selects cyclonic vorticity
        vref_trop = vref[it,ilatbound[0]:ilatbound[1],:]
        rain_trop = rain[it,ilatbound[0]:ilatbound[1],:]
        # rain_trop = rain[it*2+1,ilatbound[0]:ilatbound[1],:] # for 3hrly data

        if interp50:
            vort_trop = interp2d(vort_trop.lon, vort_trop.lat, vort_trop)(ds50.lon, ds50.lat[ilatbound_lores[0]:ilatbound_lores[1]])
            vref_trop = interp2d(vref_trop.lon, vref_trop.lat, vref_trop)(ds50.lon, ds50.lat[ilatbound_lores[0]:ilatbound_lores[1]])
            rain_trop = interp2d(rain_trop.lon, rain_trop.lat, rain_trop)(ds50.lon, ds50.lat[ilatbound_lores[0]:ilatbound_lores[1]])

        rain_threshold = np.percentile(rain_trop, rain_percentile) # recalculate for every snapshot, to account for seasonal cycle
        rain_large = rain_trop > rain_threshold # array of 1's and 0's
        labels, num_features = measure.label(rain_large)
        idx_features = range(1, num_features+1)

        meanlon = measure.mean(input=lon, labels=labels, index=idx_features) # mean longitude of an identified TC
        meanlat = measure.mean(input=lat, labels=labels, index=idx_features)
        maxvort = measure.maximum(input=vort_trop, labels=labels, index=idx_features)
        maxvref = measure.maximum(input=abs(vref_trop), labels=labels, index=idx_features)
        maxrain = measure.maximum(input=rain_trop, labels=labels, index=idx_features)

        sumlon = measure.sum(input=lon, labels=labels, index=idx_features) # used to compute size
        sizes = sumlon/meanlon

        for i in range(num_features): # i: index of identified cluster in the current time step
            if sizes[i] < size_threshold: # remove small clusters for denoising
                continue

            hasancestor = False
            for j in range(len(tracking)): # j: index of sublists in "tracking"
                ## tracking[j] = [(cftime,lat,lon,vort,size,vref,rain,year,month,day,hour), ..., (cftime,lat,lon,vort,size,vref,rain,year,month,day,hour)]
                lonprev = tracking[j][-1][2] # meanlon of a tracked TC at previous time step
                latprev = tracking[j][-1][1]
                if (meanlon[i] - lonprev)**2 + (meanlat[i] - latprev)**2 < dist_threshold**2 or \
                   (360 - abs(meanlon[i] - lonprev))**2 + (meanlat[i] - latprev)**2 < dist_threshold**2: # if within dist_threshold
                    tracking[j].append( (currtime, 
                                         meanlat[i], meanlon[i], maxvort[i], sizes[i], maxvref[i], maxrain[i],
                                         currtime.year, currtime.month, currtime.day, currtime.hour) )
                    tmp.append( tracking.pop(j) ) # move the sublist to another list, so that a TC will not have two descendants in one snapshot (note: if a TC breaks into two, the second child will be considered a brand new TC)
                    hasancestor = True
                    break # because a current TC cannot have two ancestors (note: if two TCs merge, the second one will be considered dead)
            if hasancestor == False: # if the identified TC does not have an ancestor
                tmp.append( [(currtime, 
                              meanlat[i], meanlon[i], maxvort[i], sizes[i], maxvref[i], maxrain[i],
                              currtime.year, currtime.month, currtime.day, currtime.hour)] )

        while tracking: # while tracking not empty, i.e. those TCs in "tracking" have no descendants
            completed.append( tracking.pop() ) # move sublist to completed
            if len(completed[-1]) < 4: #i remove short tracks
                completed.pop()

        tracking = tmp
        tmp = []

    while tracking: # while tracking not empty
        completed.append( tracking.pop() ) # move sublist to completed
        if len(completed[-1]) < 4: #i remove short tracks
            completed.pop()

    if return_snapshot:
        return completed, vort, rain, rain_threshold, itrange
    else:
        return completed
    
def track2netcdf(raw_track_list, yrranges):
    ## concat lists for one run
    tracks_all_yrranges = itertools.chain(*raw_track_list)
    all_yrrange = range(yrranges[0][0], yrranges[-1][-1]+1)
    ## fill the list "completed" with nan's
    completed_filled = list(itertools.zip_longest(*tracks_all_yrranges, fillvalue=[np.nan]*len(raw_track_list[0][0][0])))
    completed_array = np.array(completed_filled)
    ## sort by time
    itrack = np.argsort(completed_array[0, :, 0])
    completed_sorted = completed_array[:, itrack, :]

    yrlist = np.array([cftime.year for cftime in completed_sorted[0, :, 0]])
    iyrbeg = [np.sum(yrlist < yr) for yr in all_yrrange] # indices from which a new year starts
    ntrack_per_year = np.diff(iyrbeg+[len(yrlist)]) # number of tracks in a year
    # print(ntrack_per_year)

    n_year = len(all_yrrange)
    n_track = np.max(ntrack_per_year)
    n_lifetime = completed_sorted.shape[0]

    ds = xr.Dataset()
    variables = ['cftime', 'lat', 'lon', 'vort', 'size', 'vref', 'rain', 'year', 'month', 'day', 'hour']
    for ivar, variable in zip(range(len(variables)), variables):
        if variable == 'cftime' or variable == 'year': # don't save these variables to NetCDF
            continue 

        tmp = np.full((n_year, n_track, n_lifetime), np.nan)
        for iyr in range(len(all_yrrange)):
            data_curr = completed_sorted[:, iyrbeg[iyr]:iyrbeg[iyr]+ntrack_per_year[iyr], ivar] # track data in current year
            tmp[iyr, :ntrack_per_year[iyr], :] = np.transpose(data_curr)

        ds.__setitem__(variable, xr.DataArray(tmp, coords=[all_yrrange, range(n_track), range(n_lifetime)], dims=['year', 'track', 'lifetime']))

    return ds, all_yrrange

def run_50km(model_output_6hrly):
    ## tracker parameters
    rain_percentile = 99.5
    dist_threshold = 2.9 # degrees
    size_threshold = 3.1 # grid cells
    latlim = 30
    interp50 = False
    
    yrranges = [range(yr, yr+5) for yr in range(111, 130, 5)] # process each 5-year chuck at a time
    print(yrranges)
    
    raw_tracks = {}

    for yrrange in yrranges:
        print(yrrange)
        
        rain_in = xr.concat([xr.open_dataset(model_output_6hrly+'%04d0101.atmos_4xdaily.nc'%yr).precip for yr in yrrange], dim='time')
        vort_in = xr.concat([xr.open_dataset(model_output_6hrly+'%04d0101.atmos_4xdaily.nc'%yr).vort850 for yr in yrrange], dim='time')
        vref_in = xr.concat([xr.open_dataset(model_output_6hrly+'%04d0101.atmos_4xdaily.nc'%yr).v_ref for yr in yrrange], dim='time')

        raw_tracks[str(yrrange)] = tracker(rain_in, vort_in, vref_in, rain_percentile, dist_threshold, size_threshold, latlim, interp50=interp50)
        
    raw_track_list = [raw_tracks[str(yrrange)] for yrrange in yrranges]
    ds, all_yrrange = track2netcdf(raw_track_list, yrranges)
    encoding = {k: {'dtype': 'float32', 'zlib': True, 'complevel': 1} for k in ds.variables}
    ds.to_netcdf('seed_tracks.nc')
    
if __name__ == "__main__":
    model_output_6hrly = sys.argv[1]
    run_50km(model_output_6hrly)
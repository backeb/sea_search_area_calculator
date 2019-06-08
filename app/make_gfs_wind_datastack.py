# coding: utf-8
"""Purpose: make datastack of GFS gridded 48 hour wind forecasts and 48 hour hindcast and write to netcdf file

Dependencies: python 3.6, datetime, timedelta, pydap.client, numpy, xarray, os

Creation Date: 27 Jul 2018

Author: bjornb <bbackeberg@csir.co.za> Bjorn Backeberg

Log: <Date: Author - Comment>
    27 Jul 2018: bjornb - added code to get 48 hour forecast and 48 hour hindcast from gfs using xarray and pynio
                        - ISSUE: xarray / pynio bug: https://github.com/NCAR/pynio/issues/22
    02 Aug 2018: bjornb - rewrote code to use pydap instead of pynio (less neat but more robust)
    ___TODO____:        - implement roll over meridian, e.g. lon=ds.lon.values-180; nlon=len(lon); data=np.roll(data, nlon/2, axis=1)
"""


from datetime import datetime, timedelta
from pydap.client import open_url
import numpy as np
import xarray as xr
import pandas as pd
import os

startTime = datetime.now()

# here we set the domain we want 
# [lon-lower-left-corner, lon-upper-right-corner, lat-lower-left-corner, lat-upper-right-corner]
domain = [0, 60, -50, -10]


# here we fetch the latest GFS forcest and write ugrd10m and vgrd10m to u and v
ref_date = datetime.now()
url = "http://nomads.ncep.noaa.gov:9090/dods/gfs_0p25_1hr/gfs{date}/gfs_0p25_1hr_{hour}"
#url = "http://nomads.ncep.noaa.gov/dods/gfs_0p25_1hr/gfs{date}/gfs_0p25_1hr_{hour}"

for z in ["18z", "12z", "06z", "00z"]:
    url_z = url.format(date=ref_date.strftime('%Y%m%d'), hour=z)
    try:
        dataset = open_url(url_z)
        print("adding latest wind forcast " + url_z + " to datastack")
        break
    except:
        pass

# here we subset and load region of interest
lon = dataset['lon'];
lat = dataset['lat']; 
lonIndx = np.squeeze(np.array(np.where((lon >= domain[0]) & (lon <= domain[1]))))
latIndx = np.squeeze(np.array(np.where((lat >= domain[2]) & (lat <= domain[3]))))
ugrd10m = dataset['ugrd10m'][:,slice(latIndx[0],latIndx[-1]),slice(lonIndx[0],lonIndx[-1])]
vgrd10m = dataset['vgrd10m'][:,slice(latIndx[0],latIndx[-1]),slice(lonIndx[0],lonIndx[-1])]
time = dataset['time']
timeUnits = time.attributes['units']
del lon, lat, dataset, url_z

# here we open a xarray.Dataset to which we write ugrd10m and vgrd10m
wds = xr.Dataset()
decode_time = xr.coding.times.decode_cf_datetime # need this guy when setting the time coord in xarray.Dataset
# here we write ugrd10m and vgrd10m to the xarray.Dataset
wds['ugrd10m'] = xr.DataArray(data=ugrd10m.data[0],
        name=ugrd10m.name, dims=ugrd10m.dimensions,
        coords={'time': decode_time(time[:], timeUnits), 'lat': ugrd10m.data[2], 'lon': ugrd10m.data[3]},
        attrs=ugrd10m.attributes)
wds['vgrd10m'] = xr.DataArray(data=vgrd10m.data[0],
        name=vgrd10m.name, dims=vgrd10m.dimensions,
        coords={'time': decode_time(time[:], timeUnits), 'lat': vgrd10m.data[2], 'lon': vgrd10m.data[3]},
        attrs=vgrd10m.attributes)
del ugrd10m, vgrd10m, time


# here we fetch the last 48 hours in 6-hour timesteps
# update ref_date
ref_date = pd.to_datetime(wds.time[0].values) - timedelta(hours=6)

# now we go back and fetch the last 48 hours of GFS forecasts and add to u and v
print("Fetching and adding previous 48 hours of forecasts")
for k in range(0,8):
    url_z = url.format(date=ref_date.strftime('%Y%m%d'), hour=ref_date.strftime('%Hz'))
    dataset = open_url(url_z)
    print ("adding "+url_z+" to datastack")
    ugrd10m = dataset['ugrd10m'][slice(0,6),slice(latIndx[0],latIndx[-1]),slice(lonIndx[0],lonIndx[-1])]
    vgrd10m = dataset['vgrd10m'][slice(0,6),slice(latIndx[0],latIndx[-1]),slice(lonIndx[0],lonIndx[-1])]
    time = dataset['time'][slice(0,6)]
    del dataset, url_z
    # here we write ugrd10m and vgrd10m to the tds xarray.Dataset
    # here we open another xarray.Dataset we use to add to wds
    tds = xr.Dataset()
    tds['ugrd10m'] = xr.DataArray(data=ugrd10m.data[0],
        name=ugrd10m.name, dims=ugrd10m.dimensions,
        coords={'time': decode_time(time[:], timeUnits), 'lat': ugrd10m.data[2], 'lon': ugrd10m.data[3]},
        attrs=ugrd10m.attributes)
    tds['vgrd10m'] = xr.DataArray(data=vgrd10m.data[0],
        name=vgrd10m.name, dims=vgrd10m.dimensions,
        coords={'time': decode_time(time[:], timeUnits), 'lat': vgrd10m.data[2], 'lon': vgrd10m.data[3]},
        attrs=vgrd10m.attributes)
    del ugrd10m, vgrd10m, time
    # here we add tds to the front of wds 
    wds = xr.concat([tds, wds], dim = 'time') 
    del tds
    # update ref_date
    ref_date = pd.to_datetime(wds.time[0].values) - timedelta(hours=6)

del k

# here we write the xarray.Dataset to a netcdf file
fname = "gfs_winds.nc"
try:
    os.remove(fname)
except OSError:
    pass
print("writing datastack to netcdf4 :: "+fname)
wds.to_netcdf(fname, 'w', 'NETCDF4')

print(datetime.now() - startTime)

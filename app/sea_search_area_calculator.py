# coding: utf-8
"""Purpose: calculate search area based on input winds

Dependencies: python 2.7, gpxpy.geo, math, numpy, datetime, pynio

Creation Date: 15 Sep 2017

Author: bjornb <bbackeberg@csir.co.za> Bjorn Backeberg


Log: <Date: Author - Comment>
    28 Sep 2017: bjornb - added code to plot drift vector on basemap
    03 Oct 2017: bjornb - added code get search area radius, note dependency on gpxpy library 
                 to install do "conda install -c conda-forge gpxpy" or "pip install gpxpy --user"
    09 Oct 2017: bjornb - added code to get coordinates of the square around the search area
    10 Oct 2017: bjornb - rewrote into functions
    01 Nov 2017: bjornb - added deltaT, i.e. multiply drift vector by deltaT in hours and recalculate search area
    20 Mar 2018: bjornb - rewrote code to do "dynamic" drift vector calculation and include divergence when tdelta > 1
    11 Jul 2018: bjornb - added opendap call to get real-time winds from GFS, uses xarray and pynio
                          Note: xarray/pynio engine is buggy when working on large datasets but works fine when calling point locations
                                see issue: https://github.com/NCAR/pynio/issues/22
                          to install pynio do the following:
                            `conda create -n pynioEnv python=2.7 anaconda`
                            `source activate pynioEnv`
                            `conda install numpy pandas jupyter basemap xarray netcdf4 dask`
                            `conda install -c ncar pynio`
    02 Aug 2018: bjornb - added boolean to get_winds() to use local gfs_winds.nc datastack if available, else call from GFS opendap
"""
from __future__ import print_function

def get_resultant_coords(lon_start, lat_start, dist, brng):
    """function to calculate new lon/lat coordinate given a starting point, distance and bearing
     
    USAGE
        lon_end, lat_end = get_resultant_coords(lon_start, lat_start, drift, brng)

    INPUT
        lon_start   =	starting longitude in degrees
        lat_start   =	starting latitude in degrees
        dist        =	dist to new position in km
        brng        =	bearing, direction going towards (in 0-360 degrees)
    
    OUTPUT
        lon_end	    =	end longitude in degrees
        lat_end	    =	end latitude in degrees
    """
    import gpxpy.geo
    import math
    R = gpxpy.geo.EARTH_RADIUS/1000 # Radius of the Earth in km
    brng_rad = math.radians(brng) # convert bearing in radians
    
    # convert starting lat/lon to radians
    lat1 = math.radians(lat_start)
    lon1 = math.radians(lon_start)
    
    # get new positions
    lat2 = math.asin( math.sin(lat1)*math.cos(dist/R) + math.cos(lat1)*math.sin(dist/R)*math.cos(brng_rad))
    lon2 = lon1 + math.atan2(math.sin(brng_rad)*math.sin(dist/R)*math.cos(lat1),math.cos(dist/R)-math.sin(lat1)*math.sin(lat2))
    
    # convert back to degrees
    lat_end = math.degrees(lat2)
    lon_end = math.degrees(lon2)
    return lon_end, lat_end

def get_dist(lon1, lat1, lon2, lat2):
    """function that calculates distance between two points
    
    USAGE
        distance = get_dist(lon1, lat1, lon2, lat2)
    
    INPUT
        lon1  =   start longitude
        lat2  =   start latitude
        lon2  =   end longitude
        lat2  =   end latitude
    """
    from math import sin, cos, sqrt, atan2, radians
    import gpxpy.geo

    R = gpxpy.geo.EARTH_RADIUS/1000 # Radius of the Earth in km
    
    # convert to input lon/lat to radians 
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    distance = R * c

    return distance

def get_compass_bearing(lon1, lat1, lon2, lat2):
    """function that calculates compass bearing from point 1 to point 2
    
    USAGE
        compass_bearing = get_compass_bearing(lon1, lat1, lon2, lat2)
    
    INPUT
        lon1  =   start longitude
        lat2  =   start latitude
        lon2  =   end longitude
        lat2  =   end latitude
    """
    import math
    import numpy as np

    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)
    diffLong = math.radians(lon2 - lon1)
    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1) * math.cos(lat2) * math.cos(diffLong))
    initial_bearing = math.atan2(x, y)
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

def get_search_area_coords(lon_pts, lat_pts, sar):
    """function to calculate search area coordinates given a lon/lat points of the drift vectors
            the 4 corners of the search area in lon[0, 1, 2, 3] and lat[0, 1, 2, 3] are always given as:
            lon[0],lat[0] - 'A' the northern most location
            lon[1],lat[1] - 'B' the eastern most location
            lon[2],lat[2] - 'C' the southern most location
            lon[3],lat[3] - 'D' the western most location
    
    USAGE
        lon_sa, lat_sa = get_search_area_coords(lon_pts,lat_pts)
    
    INPUT
        lon_pts	=   longitude points of the drift vector
        lat_pts	=   latitude points of the drift vector
        sar     =   search area radius
    
    OUTPUT
        lon_sa	=   longitudes of search area
        lat_sa	=   latitudes of search area
    """
    from math import sin, cos, radians
    import numpy as np

    # get compass bearing from DATUM to Last Known Position
    if len(lon_pts.shape) > 1:
        compass_bearing = get_compass_bearing(lon_pts[0,-1], lat_pts[0,-1], lon_pts[0,0], lat_pts[0,0])
        # get centre point of furthest long axis
        lonA, latA = get_resultant_coords(lon_pts[0,-1], lat_pts[0,-1], SAR, compass_bearing-180)
        # get centre point of right hand side axis
        lonB, latB = get_resultant_coords(lon_pts[1,-1], lat_pts[1,-1], SAR, compass_bearing-90)
    else:
        compass_bearing = get_compass_bearing(lon_pts[-1], lat_pts[-1], lon_pts[0], lat_pts[0])
        # get centre point of furthest long axis
        lonA, latA = get_resultant_coords(lon_pts[-1], lat_pts[-1], SAR, compass_bearing-180)
        # get centre point of right hand side axis
        lonB, latB = get_resultant_coords(lon_pts[-1], lat_pts[-1], SAR, compass_bearing-90)

    # get distance AB
    distAB = get_dist(lonA, latA, lonB, latB)

    # get angle alpha to calculate length of adjacent and opposite
    compassAB = get_compass_bearing(lonA, latA, lonB, latB)
    alpha = 90-(compass_bearing-compassAB)

    # get 1/2 length of long axis
    long_axis = cos(radians(alpha))*distAB
    # get 1/2 length of short axis
    short_axis = 0.5*(sin(radians(alpha))*distAB+sar)

    # get coordinates of box around search area
    lon = []
    lat = []
    compass_bearing = compass_bearing - 180
    lon1, lat1 = lonA, latA
    for x in range(0, 4):
        compass_bearing = compass_bearing + 90
        if compass_bearing >= 360:
            compass_bearing = compass_bearing - 360
        if x == 0:
            lon2, lat2 = get_resultant_coords(lon1, lat1, long_axis, compass_bearing)
            lon.append(lon2)
            lat.append(lat2)
            lat1 = lat2
            lon1 = lon2
        elif x == 1 or x==3:
            lon2, lat2 = get_resultant_coords(lon1, lat1, 2*short_axis, compass_bearing)
            lon.append(lon2)
            lat.append(lat2)
            lat1 = lat2
            lon1 = lon2
        else:
            lon2, lat2 = get_resultant_coords(lon1, lat1, 2*long_axis, compass_bearing)
            lon.append(lon2)
            lat.append(lat2)
            lat1 = lat2
            lon1 = lon2
     
    # sort search area coordinates so that A is northern most location, and B is western most location
    lat = np.array(lat)
    lon = np.array(lon)
    north = np.where(lat == lat.max())
    east = np.where(lon == lon.max())
    south = np.where(lat == lat.min())
    west = np.where(lon == lon.min())
    lon_sa = [lon[north], lon[east], lon[south], lon[west]]
    lat_sa = [lat[north], lat[east], lat[south], lat[west]]

    lon_sa = np.array(lon_sa)
    lat_sa = np.array(lat_sa)
    
    return lon_sa, lat_sa

def get_time_delta(tLKP, tETA):
    """function to calculate time in hours between input start time and input end time
    
    USAGE
        tdelta = get_time_delta(tLKP, tETA)
    
    INPUT
        tLKP    =   time of last known position, in 'HH:MM' (as string)
        tETA    =   estimated time of arrival of rescue unit, in 'HH:MM' (as string)
    
    OUTPUT
        tdelta  =   time delta in hours
    """
    from datetime import datetime
    
    # define format of input
    FMT = '%H:%M'

    # calcualte time in hours
    tdelta = datetime.strptime(tETA, FMT) - datetime.strptime(tLKP, FMT)
    tdelta = round(float(tdelta.seconds)/60/60) # convert to hours
    
    return tdelta

def get_drift_and_divergence(casualty_type):
    if casualty_type == 'Person in water, state unknown':
        M, N, D = 0.01, 0.08, 40
    elif casualty_type == 'Person in water, with lifejacket':
        M, N, D = 0.02, 0.0, 45
    elif casualty_type == 'Person in water, verticle':
        M, N, D = 0.01, 0.08, 25
    elif casualty_type == 'Person in water, sitting / huddled':
        M, N, D = 0.02, 0.01, 25
    elif casualty_type == 'Person in water, floating on back':
        M, N, D = 0.02, 0.08, 40
    elif casualty_type == 'Liferafts, no ballast pockets, general type':
        M, N, D = 0.05, 0.03, 38
    elif casualty_type == 'Liferafts, no ballast pockets, no canopy, no drogue':
        M, N, D = 0.06, 0.20, 32
    elif casualty_type == 'Liferafts, no ballast pockets, with canopy, with drogue':
        M, N, D = 0.03, 0.0, 38
    elif casualty_type == 'Liferafts, shallow ballast pockets, with canopy, capsized':
        M, N, D = 0.02, -0.10, 12
    elif casualty_type == 'Liferafts, 4 to 6 man, with canopy, with drogue':
        M, N, D = 0.03, 0.04, 20
    elif casualty_type == 'Liferafts, 15 to 25 man, with canopy, with drogue':
        M, N, D = 0.04, 0.08, 15
    elif casualty_type == 'Aviation raft, 4 to 6 man, with canopy, no drogue':
        M, N, D = 0.04, 0.08, 15
    elif casualty_type == 'Sea kayak, with person':
        M, N, D = 0.01, 0.26, 20
    elif casualty_type == 'Homemake wood raft':
        M, N, D = 0.02, 0.18, 25
    elif casualty_type == 'Homemake wood raft, with sail':
        M, N, D = 0.08, 0.18, 45
    elif casualty_type == 'Surfboard with person':
        M, N, D = 0.02, 0.0, 20
    elif casualty_type == 'Windsurfer with person, sail and mast in the water':
        M, N, D = 0.03, 0.1, 16
    elif casualty_type == 'Mono hull, keel, medium displacement':
        M, N, D = 0.04, 0.0, 65
    elif casualty_type == 'Enclosed lifeboat':
        M, N, D = 0.04, -0.08, 30
    elif casualty_type == 'Vessel with outboard motors no drogue':
        M, N, D = 0.07, 0.04, 35
    elif casualty_type == 'Flat bottomed boat, Boston whaler':
        M, N, D = 0.04, 0.04, 30
    elif casualty_type == 'V hull boat':
        M, N, D = 0.03, 0.08, 25
    elif casualty_type == 'Sport fisher, centre open console':
        M, N, D = 0.06, 0.09, 30
    elif casualty_type == 'Commercial fishing vessel type unknown':
        M, N, D = 0.04, 0.02, 65
    elif casualty_type == 'Commercial fishing vessel longline, stern or net':
        M, N, D = 0.04, 0.0, 65
    elif casualty_type == 'Coastal freighter':
        M, N, D = 0.03, 0.0, 65
    elif casualty_type == 'Fishing vessel general debris':
        M, N, D = 0.02, 0.0, 15
    elif casualty_type == 'Cubic meter bait box, loading unknown':
        M, N, D = 0.02, 0.02, 30
    else:
        print('no such casualty type exists')
    return M, N, D


def calc_spd_drctn(u, v):
    """function to calculate speed and direction from u and v components. Return speed in knots
    
    USAGE
        spd, drctn = calc_spd_drctn(u, v)
    
    INPUT
        u = u-component velocity
        v = v-component velocity
        
    OUTPUT
        spd   = speed
        drctn = direction in degrees from north
    """
    import numpy as np
    
    spd = np.sqrt(u**2 + v**2) * 1.94384 # convert to knots
    drctn_trig_to = np.arctan2(u/spd, v/spd)
    drctn_trig_to_degrees = drctn_trig_to * 180/np.pi
    drctn_trig_from_degrees = drctn_trig_to_degrees + 180
    drctn = drctn_trig_from_degrees
    del drctn_trig_to, drctn_trig_to_degrees, drctn_trig_from_degrees
    
    return spd, drctn

def get_winds(lonIn, latIn, timeIn):
    """"function to get real time winds 10m above surface from NCEP for point location
            List of available real time data: http://nomads.ncep.noaa.gov/

        USAGE
            u, v = get_winds(lonIn, latIn, timeIn)

        INPUT
            lonIn  = longitude point location for which to extract the winds
            latIn  = latitude point location for which to extract the winds
            timeIn = time for which to extract the winds, will pull nearest forecast time
                     format = HH:MM

        OUTPUT
            u = u-component winds at 10m above surface for specified lonIn, latIn, timeIn
            v = v-component winds at 10m above surface for specified lonIn, latIn, timeIn
    """
    import numpy as np
    import xarray as xr
    from datetime import datetime
    import os.path

    dateIn = datetime.now()
    dateIn = dateIn.replace(hour = int(timeIn[0:2]), minute = int(timeIn[3:5]), second = 0, microsecond  = 0)

    # try get winds from local OCIMS datastack else pull from NCEP 
    ds = None
    if os.path.isfile("gfs_winds.nc"):
        ds = xr.open_dataset("gfs_winds.nc")\
                .sel(lon = lonIn, lat = latIn, method = "nearest")\
                .sel(time = dateIn.strftime('%Y-%m-%dT%H:%M:%S'), method = "nearest")
        print("Fetched winds for "+str(ds.time.values)+" at "+str(ds.lon.values)+", "+str(ds.lat.values))  
    else:
        print("Locally downloaded datastack not available, pulling data from NCEP")
        for z in ["18z", "12z", "06z", "00z"]:
            url = "http://nomads.ncep.noaa.gov:9090/dods/gfs_0p25_1hr/gfs"+dateIn.strftime('%Y%m%d')+"/gfs_0p25_1hr_"+z
            try:
                ds = xr.open_dataset(url,  engine = 'pynio')\
                                    .sel(lon = lonIn, lat = latIn, method = "nearest")\
                                    .sel(time = dateIn.strftime('%Y-%m-%dT%H:%M:%S'), method = "nearest")
                print(url)
                print("Fetched winds for "+str(ds.time.values)+" at "+str(ds.lon.values)+", "+str(ds.lat.values))  
                break
            except:
                pass

    u = ds.ugrd10m.values
    v = ds.vgrd10m.values
    #time = ds.time

    return u, v


if __name__ == '__main__':
    from numpy import array, zeros
    from datetime import datetime, timedelta
    
    # Last Known Position and Time
    lat_LKP = array(-34.358225)
    lon_LKP = array(18.500566)
    time_LKP = '16:30'
    
    # Rescue Unit ETA
    time_ETA = '19:30'

    CasualtyType = 'Sea kayak, with person'

    # get drift and divergence for casualty type
    M, N, D = get_drift_and_divergence(CasualtyType)
    
    # get time delta
    tdelta = get_time_delta(time_LKP, time_ETA)
    
    # initiatlise some variables
    tdrft = 0.0 # record total drift distance in nautical miles
    if tdelta == 1:
        lon_PTS = zeros([int(tdelta)+1]); lon_PTS[0] = lon_LKP
        lat_PTS = zeros([int(tdelta)+1]); lat_PTS[0] = lat_LKP
    else:
        lon_PTS = zeros([3, int(tdelta)+1]); lon_PTS[:,0] = lon_LKP
        lat_PTS = zeros([3, int(tdelta)+1]); lat_PTS[:,0] = lat_LKP
    
    # CALCULATE DRIFT DISTANCE DYNAMICALLY IN LOOP FOR LENGTH OF TDELTA
    for i in range(0,int(tdelta)):
    
        # sort out time for get_winds()
        timeIn = (datetime.strptime(time_LKP, '%H:%M')+timedelta(hours=i)).strftime('%H:%M')
        # get wind speed and direction - currently pulling from GFS
        # u- and v-component winds in m/s
        if tdelta == 1:
            u10, v10 = get_winds(lon_PTS[i], lat_PTS[i], timeIn)
        else:
            u10, v10 = get_winds(lon_PTS[0, i], lat_PTS[0, i], timeIn)
        # wind speed in knots and direction in degrees coming from
        wnd_spd, wnd_drctn = calc_spd_drctn(u10, v10)
        # wind direction is always coming from, i.e. in this case southerly (coming from the south)
        drctn = wnd_drctn - 180 # convert to going towards

        # drift per hour in nautical miles
        drft = wnd_spd * M + N
        # convert distance from nautical miles to km
        drft = drft * 1.852 # distance in km
    
        # record total drift distance
        tdrft = tdrft + drft
        
        """TODO get ocean currents"""
        """TODO add ocean current as 1:1 drift component to wind drift"""

        if tdelta == 1:
            lon_PTS[i+1], lat_PTS[i+1] = get_resultant_coords(lon_PTS[i], lat_PTS[i], drft, drctn)
        else:
            lon_PTS[0, i+1], lat_PTS[0, i+1] = get_resultant_coords(lon_PTS[0, i], lat_PTS[0, i], drft, drctn)
            lon_PTS[1, i+1], lat_PTS[1, i+1] = get_resultant_coords(lon_PTS[1, i], lat_PTS[1, i], drft, drctn + D)
            lon_PTS[2, i+1], lat_PTS[2, i+1] = get_resultant_coords(lon_PTS[2, i], lat_PTS[2, i], drft, drctn - D)
        
    # define the DATUM
    if tdelta > 1:
        lon_DATUM, lat_DATUM = lon_PTS[0,-1], lat_PTS[0,-1]
    else:
        lon_DATUM, lat_DATUM = lon_PTS[-1], lat_PTS[-1]
        
    """TODO add get_wave call here and include in plot_SAR_map"""

    # define the search area radius as 1/3 of the distance of the total drift
    SAR = tdrft/3
        
    # get search area coordinates
    lon_SA, lat_SA = get_search_area_coords(lon_PTS, lat_PTS, SAR)
    
    print("Suggested search area")
    print("DATUM:   Longitude = "+str(lon_DATUM)+", Latitude = "+str(lat_DATUM))
    print("A:       Longitude = "+str(lon_SA[0])+", Latitude = "+str(lat_SA[0]))
    print("B:       Longitude = "+str(lon_SA[1])+", Latitude = "+str(lat_SA[1]))
    print("C:       Longitude = "+str(lon_SA[2])+", Latitude = "+str(lat_SA[2]))
    print("D:       Longitude = "+str(lon_SA[3])+", Latitude = "+str(lat_SA[3]))



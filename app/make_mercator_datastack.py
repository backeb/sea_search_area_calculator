# coding: utf-8
"""Purpose: download latest hourly ocean forecasts from CMEMS GLOBAL-ANALYSIS-FORECAST-PHY-001-024-HOURLY-T-U-V-SSH

Dependencies: python 2.7, motu-client (do 'pip install motu-client'), datetime, timedelta

Creation Date: 21 Aug 2018

Author: bjornb <backeb@gmail.com> Bjorn Backeberg

Log: <Date: Author - Comment>
    21 Aug 2018: bjornb - adapted script from Mostafa Bakhoday-Paskyabi <Mostafa.Bakhoday@nersc.no>
"""
from __future__ import print_function
import os
from datetime import datetime, timedelta

startTime = datetime.now()

# here we define all the bits and pieces that get put in the 'runcommand' variable used to call the motu-client
path2motuClient = '/Users/bjornb/miniconda3/envs/motuEnv/lib/python2.7/site-packages'

usrname = 'bbackeberg'
passwd = 'iaTmwJ7D'

# [lon-lower-left-corner, lon-upper-right-corner, lat-lower-left-corner, lat-upper-right-corner]
domain = [0, 60, -50, -10]

startDate = datetime.now() - timedelta(hours = 48)
startDate = startDate.replace(minute = 30, second = 0, microsecond  = 0)
endDate = datetime.now() + timedelta(hours = 120)
endDate = endDate.replace(minute = 30, second = 0, microsecond  = 0)

# thetao = Temperature in degrees C, zos = SSH in m, uo = Eastward velocity in m/s, vo = Northward velocity in m/s
varList = ['thetao', 'zos', 'uo', 'vo']

# NOTE only surface fields available hourly
depths = [0.493, 0.4942]

path2saveData = os.getcwd()+'/'
fname = 'mercator_ocean.nc'

# create the runcommand string
runcommand = 'python '+path2motuClient+'/motu-client.py --quiet'+ \
        ' --user '+usrname+' --pwd '+passwd+ \
        ' --motu http://nrt.cmems-du.eu/motu-web/Motu'+ \
        ' --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS'+ \
        ' --product-id global-analysis-forecast-phy-001-024-hourly-t-u-v-ssh'+ \
        ' --longitude-min '+str(domain[0])+' --longitude-max '+str(domain[1])+ \
        ' --latitude-min '+str(domain[2])+' --latitude-max '+str(domain[3])+ \
        ' --date-min "'+str(startDate.strftime('%Y-%m-%d %H:%M:%S'))+'" --date-max "'+str(endDate.strftime('%Y-%m-%d %H:%M:%S'))+'"'+ \
        ' --depth-min '+str(depths[0])+' --depth-max '+str(depths[1])+ \
        ' --variable '+varList[0]+' --variable '+varList[1]+' --variable '+varList[2]+' --variable '+varList[3]+ \
        ' --out-dir '+path2saveData+' --out-name '+fname

# run the runcommand, i.e. download the data specified above
print('fetching latest mercator ocean forecast from CMEMS and making datastack')
os.system(runcommand)

print(datetime.now() - startTime)


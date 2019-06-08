```
Date: Author - Comment  
28 Sep 2017: backeb - added code to plot drift vector on basemap  
03 Oct 2017: backeb - added code get search area radius, note dependency on gpxpy library, to install do   
                      `conda install -c conda-forge gpxpy`  
                      or   
                      `pip install gpxpy --user`  
09 Oct 2017: backeb - added code to get coordinates of the square around the search area  
10 Oct 2017: backeb - rewrote into functions  
01 Nov 2017: backeb - added deltaT, i.e. multiply drift vector by deltaT in hours and recalculate search area  
20 Mar 2018: backeb - rewrote code to do "dynamic" drift vector calculation and include divergence when tdelta >1  
11 Jul 2018: backeb - added opendap call to get real-time winds from GFS, uses xarray and pynio  
                      Note: xarray/pynio engine is buggy when working on large datasets but works fine when  
                      calling point locations. see issue: https://github.com/NCAR/pynio/issues/22  
                      to install pynio do the following:  
                      ```conda create -n pynioEnv python=2.7 anaconda  
                      source activate pynioEnv  
                      conda install numpy pandas jupyter basemap xarray netcdf4 dask  
                      conda install -c ncar pynio```   
02 Aug 2018: backeb - added boolean to get_winds() to use local gfs_winds.nc datastack if available, else call  
                      from GFS opendap  
08 Jun 2019: backeb - UPATE: pynio now works in python 3 - dependencies updated  
                    - Fixed bug in get_winds function and made notebook more interactive  
                    - added widget for input parameters  
                    - writing positions out to pandas dataframe  
```

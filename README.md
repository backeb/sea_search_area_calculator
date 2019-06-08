# Sea Search Area Calculator
Author(s): [Bjorn Backeberg](mailto:backeb@gmail.com) (backeb) <br>
Creation date: 15-Sep-2017 <br>
Last updated:  08-Jun-2019 <br>

## Purpose
Use real-time wind and ocean current (coming soon) forecasts to calculate search areas coordinates for persons/objects lost at sea.

This notebook is based on the search area calculator I originally developed for the [OCIMS](https://www.ocims.gov.za/) "Operations at Sea" decision support tool in collaboration with the [NSRI](https://www.nsri.org.za/).

## Dependencies
Manage python installation through [conda](https://docs.conda.io/en/latest/miniconda.html)

Currently using Python 3.7.1

`conda install numpy pandas jupyter basemap xarray netcdf4 dask`  
`conda install -c conda-forge gpxpy`  
`conda install -c conda-forge pynio`  
`conda install -c conda-forge pygrib`  
`conda install -c conda-forge ipyleaflet`  

## Directory tree
```
|- CHANGELOG.md
|- README.md
|- Docs:                                    Relevant documentation
|- app:                                     Main code base
    |- make_gfs_wind_datastack.py:          makes wind forecast datastack for sea_search_area_calculator.py
    |- make_mercator_datastack.py:          makes ocean current forecast datastack for sea_search_area_calculator.py
    |- plot_search_area_map.py:             plot output from sea_search_area_calculator.py
    |- sea_search_area_calculator.py:       calculate search area based on wind and ocean current input data
|- notebooks:                               Development notebooks
    |- sea_search_area_calculator.ipynb:    interactive notebook for sea_search_area_calculator.py
```

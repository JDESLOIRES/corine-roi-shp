# Instructions

The purpose of this small module is to get a shapefile from Corine Land Cover (2018) given a bounding box distance and a GPS coordinate.

Corine Land Cover https://land.copernicus.eu/pan-european/corine-land-cover/clc2018?tab=download

# Description

Given the distance (e.g. 500 meters) around a GPS coordinate point in WGS CRS, it will return a geopandas.GeoDataFrame with the subset of Corine.
The script example.py shows how to call the module.

# Environment

To use the libraries required to run the script, you can directly create the environment from the .yml file :
```python
cd corine-roi-shp
conda env create -f environment.yml
conda activate corine-env
```
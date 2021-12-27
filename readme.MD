# Instructions

The purpose of this small module is to get a shapefile from Corine Land Cover (2018) given a bounding box distance and a GPS coordinate.

Corine Land Cover https://land.copernicus.eu/pan-european/corine-land-cover/clc2018?tab=download

# Description

Given the radius (e.g. 500 meters)around a GPS coordinate point in WGS CRS :

- _**preprocess_corine()**_ should be used the first time as it. It converts raw data of Corine into the spatial reference EPSG:4326 and delienates your area of interest. Corine Land Cover is pretty wide; it inscludes even european islands off Africa (below the picture).
- _**add_labels()**_ just add to your shapefile additional columns with respect to an input csv file that includes the definition of the Land Cover categories (from 1 to 44). 
  In this file, we have the column _Code_ (Corine Land Cover) and _libelle_fr_ (definition in french) and groups (custom second level of grouping). 
This function also computes the relative surface of each Code and groups.
  
# Environment

To use the libraries required to run the script, you can directly create the environment from the .yml file :
```python
cd corine-roi-shp
conda env create -f environment.yml
conda activate corine-env
'''
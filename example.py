import matplotlib.pyplot as plt
import geopandas as gpd
from gdal_process import CLCRadius

###########################################################################################
#Corine Land Cover
#Get contours of France

france = gpd.read_file('./France_shapefile/fr_10km.shp')
france['ID'] = 1
extent = france.dissolve('ID').to_crs({'init': 'epsg:4326'})

#Get the distance for the bounding box
radius = 1000

#Take the coordinate that will be use to extract CLC
lon = -3.81
lat = 48.68

#Define the pipeline
pipeline = CLCRadius(
    radius = radius)

#Do it only if it is the first time :  project into degrees + subset the .tif file given the extent of France
#Corine Land Cover https://land.copernicus.eu/pan-european/corine-land-cover/clc2018?tab=download
#here, the folder with .tif corine file was located in the current git repository
pipeline._change_crs(
    extent = extent,
    path_raster=r'./CLC_2018/u2018_raster100m/U2018_CLC2018_V2020_20u1.tif')

#Save a shapefile with raw data from the CLC labels as vector file
gdf_corine = pipeline.vectorize_corine_around_point(lon, lat)
gdf_corine = gdf_corine[gdf_corine.Code>0]
gdf_corine.plot('Code', legend = True)
plt.show()

output = gpd.read_file( './CLC_2018/area/area.shp')

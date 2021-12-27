import rasterio
import geopandas as gpd

from osgeo import gdal, ogr
import numpy as np

import os
import pandas as pd
from shapely.geometry import Point

import utils
import shutil
import earthpy.spatial as es
import pathlib
############################################################################################

class CLCRadius:

    def __init__(self, radius):
        '''
        Parameter to initialize the Corine land cover polygonized
        Args:
            extent (shp file) : extent that crop the Corine land cover (here, France)
            radius: distance (int) for the bounding box output (in meters)
            origin_folder_path (str) : path where the python module is saved
        '''
        self.radius = radius

    @staticmethod
    def get_subdirectory_file(path):
        return '/'.join(path.split('/')[:-1])

    def _change_crs(self,
                    extent,
                    path_raster =r'./CLC_2018/u2018_raster100m/U2018_CLC2018_V2020_20u1.tif'
                    ):
        '''
        From the original file download on Corine Land Cover, perform the following preprocessing:
            - Change the projection to match with degrees coordinates (epsg:4326)
            - Crop the image given an extent (ex. Contours of France)
        Args:
            path_raster (str) : path where the CLC is saved
            extent (geopandas.GeoDataFrame) : dataframe with a single row (polygon of extent)
            output_raster (str) : path where the output file will be saved
        '''
        print('Change crs from original Corine .tif file from 3035 to 4326')
        os.chdir(str(pathlib.Path().resolve()))
        input = gdal.Open(path_raster)
        path_in = self.get_subdirectory_file(path_raster)

        if not os.path.exists(path_raster):
            # change CRS into degrees
            gdal.Warp(os.path.join(path_in, 'tempo_4326.tif'),
                      input,
                      srcSRS="+init=epsg:3035",
                      dstSRS="+init=epsg:4326")

            # Compress the results
            topts = gdal.TranslateOptions(format='GTiff',
                                          outputType=gdal.GDT_UInt16,
                                          creationOptions=['COMPRESS=LZW', 'GDAL_PAM_ENABLED=NO'],
                                          bandList=list(range(1, 2)))  # gdal.GDT_Byte

            gdal.Translate(os.path.join(path_in,'U2018_CLC2018_V2020_20u1_4326.tif'),
                           os.path.join(path_in,'tempo_4326.tif'), options=topts)
            os.remove(os.path.join(path_in,'tempo_4326.tif'))
        if not extent is None:
            # crop the image given the extent of a country (here, France)
            es.crop_all([os.path.join(path_in,'U2018_CLC2018_V2020_20u1_4326.tif')],
                        path_in, extent,
                        overwrite=True, all_touched=True, verbose=True)

    @staticmethod
    def _radius_around_points(lon, lat, radius=500, save=True):
        '''
        Create a bounding box of n meters given a longitude and latitude
        Args:
            lon (float): longitude
            lat (float) : latitude
            radius (int) : distance
            save (boolean) : save the output file

        Returns:
            gpd.GeoDataFrame with bounding box corresponding to the radius around the coordinate given

        '''
        print('Draw radius')
        df = pd.DataFrame(dict(X_COORD=[lon], Y_COORD=[lat]))
        df['coordinates'] = [Point(x, y) for x, y in zip(df.X_COORD, df.Y_COORD)]

        df = gpd.GeoDataFrame(df, geometry=df.coordinates)

        df['circus_' + str(radius)] = df.geometry.apply(lambda x: utils.create_polygon_radius(x, radius))

        df = gpd.GeoDataFrame(df[['X_COORD', 'Y_COORD', 'circus_' + str(radius)]],
                              geometry=df['circus_' + str(radius)], crs="EPSG:4326")

        df['Index'] = range(1, df.shape[0] + 1)
        df = df.drop(['circus_' + str(radius)], axis=1)

        if save:
            output = ('./radius_' + str(radius))
            if not os.path.exists(output):
                os.makedirs(output)
            df.to_file(os.path.join(output,'radius_' + str(radius) + '.shp'))

        return df

    def _mask_corine_around_point(self, lon, lat,
                                  path_corine_wgs84 =r'./CLC_2018/u2018_raster100m/U2018_CLC2018_V2020_20u1_4326_crop.tif'):
        '''
        Get raw data labels from Corine Land Cover (from 1 to 44 classes)
        Args:
            lon (float): longitude
            lat (float): latitude
            output_folder: folder path where the shapefile will be saved

        Returns:
            gpd.GeoDataFrame with polygons corresponding for each Code of Corine Land Cover from the ROI

        '''
        print('Mask Corine around radius')

        path_in = self.get_subdirectory_file(path_corine_wgs84)
        with rasterio.open(path_corine_wgs84) as src0:
            meta = src0.meta
            meta['nodata'] = 0
            meta['dtype'] = 'uint16'
            meta['count'] = 1

        if str(meta['crs']) == 'EPSG:3035':
            raise ValueError('You must convert CLC into EPSG:4326 using preprocess_corine()')

        # get a shapefile with each traps given a zone around
        self._radius_around_points(lon, lat, self.radius, save=True)

        # convert this shapefile into raster
        out = utils.Vector_to_Raster(Dir='./',
                                     vector_path='./radius_' + str(self.radius) + '/radius_' + str(self.radius) + '.shp',
                                     reference_file=path_corine_wgs84,
                                     attribute='Index',
                                     file_name='radius_' + str(self.radius))

        # remove intermediary files : only out interests us
        os.remove('./radius_' + str(self.radius) + '.tif')

        # Apply a mask if a pixel is outself.radius
        clc = gdal.Open(path_corine_wgs84).ReadAsArray()
        mask = (out == 0)

        x = np.ma.array(clc,
                        dtype=np.int16,
                        mask=(mask).astype(bool),
                        fill_value=0)
        results = x.filled()
        del clc


        # Write results
        with rasterio.open(os.path.join(path_in, 'tempo.tif'), 'w', **meta) as dst:
            dst.write_band(1, results.astype(np.uint16))

        # Compress the results
        topts = gdal.TranslateOptions(format='GTiff',
                                      outputType=gdal.GDT_UInt16,
                                      creationOptions=['COMPRESS=LZW', 'GDAL_PAM_ENABLED=NO'],
                                      bandList=list(range(1, 2)))  # gdal.GDT_Byte

        gdal.Translate(os.path.join(path_in, 'tempo' + '_' + str(self.radius) + '.tif'),
                       os.path.join(path_in, 'tempo.tif'), options=topts)
        os.remove(os.path.join(path_in, 'tempo.tif'))
        ############################################################

    def _get_intersection_area(self, output):
        print('Get intersection polygon with area')
        ###########################################################################################################
        # Do the intersection between Land Cover and GPS coordinates
        df = gpd.read_file(os.path.join('./radius_' + str(self.radius),
                                        './radius_' + str(self.radius) + '.shp'))
        shutil.rmtree('./radius_' + str(self.radius))

        output.crs = 'epsg:4326'
        join_left_df = gpd.overlay(df, output, how='intersection')

        return join_left_df

    def vectorize_corine_around_point(self, lon, lat,
                                      path_corine_wgs84=r'./CLC_2018/u2018_raster100m/U2018_CLC2018_V2020_20u1_4326_crop.tif'):
        print('Vectorize output file')
        path_in = self.get_subdirectory_file(path_corine_wgs84)
        self._mask_corine_around_point(lon, lat,
                                       path_corine_wgs84= path_corine_wgs84)

        if not os.path.exists('./tempo_shp'):
            os.makedirs('./tempo_shp')

        #Return a shapefile corresponding to the subset of Corine Land Cover
        sourceRaster = gdal.Open(os.path.join(os.getcwd(), os.path.join(path_in, 'tempo' + '_' + str(self.radius) + '.tif')))
        band = sourceRaster.GetRasterBand(1)

        driver = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(os.path.join('./tempo_shp', 'tempo_shp.shp')):
            driver.DeleteDataSource(os.path.join(os.path.join('./tempo_shp', 'tempo_shp.shp')))
        
        outDatasource = driver.CreateDataSource(os.path.join('./tempo_shp', 'tempo_shp.shp'))
        outLayer = outDatasource.CreateLayer('tempo_shp', srs=None)
        newField = ogr.FieldDefn('DN', ogr.OFTInteger)
        outLayer.CreateField(newField)
        gdal.Polygonize(band, None, outLayer, 0, [], callback=None)
        outDatasource.Destroy()
        
        os.remove(os.path.join(path_in, 'tempo' + '_' + str(self.radius) + '.tif'))

        #Do an intersection with the input radius to get only its intersected polygons
        output = gpd.read_file(os.path.join('./tempo_shp', 'tempo_shp.shp'))
        shutil.rmtree('./tempo_shp')
        gpd_output = self._get_intersection_area(output)
        gpd_output = gpd.GeoDataFrame(gpd_output, geometry=output.geometry)
        gpd_output = gpd_output.sort_values('Index')
        gpd_output = gpd_output.rename({'DN': 'Code'}, axis=1)
  
        return gpd_output

    def add_labels(self,
                   gpd_output,
                   regrouping_file = './CLC_2018/regroupement_classe_corinne_01_06_21.csv'):

        '''
        Add labels from the created .csv file with columns Code (CLC2018) and groups (name of the groups)
        Args:
            output_folder (str) : Folder where the shapefile is created
            regrouping_file (str) : csv file where we have labels

        Returns:
            gpd.GeoDataFrame with diversity indicators (relative surfaces) by code and groups
        '''

        # Aggregate the labels (44) from Corine into specific given classes define in the files 'regroupement_classe_corinne_01_06_21'
        regrouping = pd.read_csv(regrouping_file, encoding="ISO-8859-1", sep=',')
        regrouping.Code = regrouping.Code.astype(str)
        gpd_output.Code = gpd_output.Code.astype(str)
        gpd_output = pd.merge(gpd_output, regrouping[['libelle_fr', 'Code', 'groups']],
                              on='Code', how='left')

        # Compute relative surf from each class
        surf_agg = gpd_output[['Index', 'groups', 'surf_obj', 'surf_rel']].groupby(
            ['Index', 'groups']).agg('sum')

        surf_agg.reset_index(inplace=True)
        surf_agg.columns = ['Index', 'groups', 'surf_obj_g', 'surf_rel_g']

        gpd_output = pd.merge(gpd_output, surf_agg, on=['Index', 'groups'], how='left')
        gpd_output.crs = 'epsg:4326'

        return gpd_output


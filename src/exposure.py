import os
import pandas as pd
import numpy as np
import geopandas as gpd
import pygeos
from osgeo import ogr,gdal
from tqdm import tqdm
from pygeos import from_wkb
from rasterstats import zonal_stats
import rasterio
import rioxarray


def overlay_hazard_data(province):
    building_path = os.path.join("C:\\Data",'OSM','feather')
    forecast_path = os.path.join("C:\\Data",'forecast')
    
    building_province_path = os.path.join(building_path,'{}.ft'.format(province))
    
    # read buildings
    df_buildings = pd.read_feather(building_province_path)
    df_buildings.geometry = pygeos.from_wkb(df_buildings.geometry.values)
    
    # read hazard data in polygon format
    rds = rioxarray.open_rasterio(hazard_file)
    rds.name = "data"
    df_ds = rds.squeeze().to_dataframe().reset_index()
    df_ds['geometry'] = [pygeos.points(x) for x in list(zip(df_ds['x'],df_ds['y']))]
    df_ds['geometry'] = pygeos.get_parts(pygeos.voronoi_polygons(pygeos.multipoints(df_ds.geometry)))
    df_ds = df_ds.drop(['y','x','band','spatial_ref'],axis=1)
   
    #overlay data
    tqdm.pandas()
    tree = pygeos.STRtree(df_ds.geometry)
    df_buildings['windspeed'] = df_buildings.geometry.progress_apply(lambda x: df_ds.iloc[tree.query(x,predicate='intersects')[0]]['data'])
    
    df_buildings['province'] = province
    #return file
    return df_buildings


def overlay_hazard_data_rasterstats(province):
    building_path = os.path.join("C:\\Data",'OSM','feather')
    forecast_path = os.path.join("C:\\Data",'forecast')
    
    building_province_path = os.path.join(building_path,'{}.ft'.format(province))
    
    df_buildings = pd.read_feather(building_province_path)
    df_buildings.geometry = pygeos.from_wkb(df_buildings.geometry.values)
    
    gdf = gpd.GeoDataFrame(df_buildings.copy())
    
    with rasterio.open(os.path.join(forecast_path,'GeoTIFF',"f033_NL.tif")) as src:
        affine = src.transform
        array = src.read(1)
    
    tqdm.pandas()
    gdf['max_speed'] = gdf.geometry.progress_apply(lambda x : zonal_stats(x, array, affine=affine,all_touched=True,stats="max",nodata=9999)[0]['max'])
    
    return gdf
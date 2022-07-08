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
import pyproj

from multiprocessing import Pool,cpu_count

def reproject_assets(df_ds,current_crs="epsg:4326",approximate_crs = "epsg:28992"):
    """[summary]

    Args:
        df_ds ([type]): [description]
        current_crs (str, optional): [description]. Defaults to "epsg:3857".
        approximate_crs (str, optional): [description]. Defaults to "epsg:4326".

    Returns:
        [type]: [description]
    """    

    geometries = df_ds['geometry']
    coords = pygeos.get_coordinates(geometries)
    transformer=pyproj.Transformer.from_crs(current_crs, approximate_crs,always_xy=True)
    new_coords = transformer.transform(coords[:, 0], coords[:, 1])
    
    return pygeos.set_coordinates(geometries.copy(), np.array(new_coords).T) 

def get_use(x,omschrijving,tree):
    try:
        return omschrijving[tree.query(x,predicate='intersects')][0]
    except:
        return 'overig landgebruik'
    
def update_building_information(province):
    building_path = os.path.join("C:\\Data",'OSM','feather')
    cbs_path = os.path.join("C://","Data","Bodemgebruik")

    building_province_path = os.path.join(building_path,'{}.ft'.format(province))

    # read df_buildings
    df_buildings = pd.read_feather(building_province_path)
    df_buildings.geometry = pygeos.from_wkb(df_buildings.geometry.values)
    df_buildings.geometry = reproject_assets(df_buildings)

    BBG2015 = gpd.read_file(os.path.join(cbs_path,'BBG2015.shp'))

    print('BBG loaded for {}'.format(province))

    landuse_to_keep = ['Spoorweg', 'Hoofdweg', 'Vliegveld', 'Woongebied',
           'Detailhandel en horeca', 'Openbare voorziening',
           'Sociaal-culturele voorziening', 'Bedrijfsterrein', 'Stortplaats',
           'Wrakkenopslagplaats', 'Begraafplaats', 'Delfstofwinplaats',
           'Bouwterrein', 'Semi verhard overig terrein', 'Park en plantsoen',
           'Sportterrein', 'Volkstuin', 'Dagrecreatief terrein',
           'Verblijfsrecreatie', 'Glastuinbouw', 'Overig agrarisch gebruik']

    BBG2015 = BBG2015.loc[BBG2015.Omschrijvi.isin(landuse_to_keep)].reset_index()    

    df_BBG = pd.DataFrame(BBG2015[['Omschrijvi','geometry']].copy())
    df_BBG.geometry = pygeos.from_shapely(df_BBG.geometry)

    province_bbox = pygeos.envelope(pygeos.geometrycollections(df_buildings.geometry.values))
    province_BBG = df_BBG.loc[pygeos.intersects(df_BBG.geometry,province_bbox)].reset_index()
    
    del BBG2015,df_BBG

    omschrijving = province_BBG.Omschrijvi.values

    tree = pygeos.STRtree(province_BBG.geometry)

    df_buildings['centroid'] = pygeos.centroid(df_buildings.geometry)
        
    tqdm.pandas(desc=province)
    df_buildings['BBG'] = df_buildings['centroid'].progress_apply(lambda x: get_use(x,omschrijving,tree))

    df_buildings.geometry = reproject_assets(df_buildings,current_crs="epsg:28992",approximate_crs = "epsg:4326")
    df_buildings = df_buildings.drop(['centroid'],axis=1)
    df_buildings.geometry = pygeos.to_wkb(df_buildings.geometry)

    df_buildings.to_feather(building_province_path)

    return df_buildings

def overlay_hazard_data(province):
    building_path = os.path.join("C:\\Data",'OSM','feather')
    forecast_path = os.path.join("C:\\Data",'forecast')
    
    building_province_path = os.path.join(building_path,'{}.ft'.format(province))
    
    hazard_file = os.path.join(forecast_path,'GeoTIFF',"f033_NL.tif")

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
    tqdm.pandas(desc=province)
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
    
    tqdm.pandas(desc=province)
    gdf['max_speed'] = gdf.geometry.progress_apply(lambda x : zonal_stats(x, array, affine=affine,all_touched=True,stats="max",nodata=9999)[0]['max'])
    
    return gdf


if __name__ == '__main__':       

    osm_files = os.path.join("C:\\Data",'OSM','feather')
    provinces = [x.split('.')[0] for x in os.listdir(osm_files)]

    # # in parallel get exposure information
    with Pool(processes=cpu_count()) as pool: 
         pool.map(update_building_information,provinces,chunksize=1)           
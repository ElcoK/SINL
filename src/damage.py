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

def overlay_hazard_data(province,model='ECMWF',fc_time='041'):
    building_path = os.path.join("C:\\Data",'OSM','feather')
    forecast_path = os.path.join("C:\\Data",'forecast')
    
    building_province_path = os.path.join(building_path,'{}.ft'.format(province))
    
    # read buildings
    df_buildings = pd.read_feather(building_province_path)
    df_buildings.geometry = pygeos.from_wkb(df_buildings.geometry.values)
    
    hazard_file = os.path.join(forecast_path,'GeoTIFF',"{}_f{}_NL.tif".format(model,fc_time))
        
    # read hazard data in polygon format
    rds = rioxarray.open_rasterio(hazard_file)
    rds.name = "data"
    df_ds = rds.squeeze().to_dataframe().reset_index()
    df_ds['geometry'] = [pygeos.points(x) for x in list(zip(df_ds['x'],df_ds['y']))]
    df_ds['geometry'] = pygeos.get_parts(pygeos.voronoi_polygons(pygeos.multipoints(df_ds.geometry)))
    df_ds = df_ds.drop(['y','x','band','spatial_ref'],axis=1)
   
    #overlay data
    tqdm.pandas(desc='hazard overlay {}'.format(province))
    tree = pygeos.STRtree(df_ds.geometry)
    df_buildings['fg10_{}_{}'.format(model,fc_time)] = df_buildings.geometry.progress_apply(lambda x: df_ds.iloc[tree.query(x,predicate='intersects')[0]]['data'])
    
    df_buildings['province'] = province
    #return file
    return df_buildings

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

def get_m2(df):
    
    df.geometry = reproject_assets(df)
    
    return pygeos.area(df.geometry)

def damage_estimate(province,model='GFS',fc_time='018'):
    
    damage_path = os.path.join("C:\\Data",'forecast','damage')

    damage_province_path = os.path.join(damage_path,'{}_{}_{}.ft'.format(province,model,fc_time))
    
    overlaid = overlay_hazard_data(province,model=model,fc_time=fc_time)
    overlaid['m2'] = get_m2(overlaid)
    
    damage_functions = pd.read_excel(os.path.join('..','data','wind_damage_functions.xlsx'),sheet_name='building_class_global',usecols="A:E")
    max_damage = pd.read_excel(os.path.join('..','data','wind_damage_functions.xlsx'),sheet_name='building_class_global',usecols="A:E")
    
    curve = 'RM1L'
    
    max_damage = pd.read_excel(os.path.join('..','data','MaxDAM_CBS.xlsx'),sheet_name='Sheet1',usecols="A:C")
    max_damage.columns = ['BBG','Rotterdam','Top of curve (60%)']
    max_dam_dict = dict(zip(max_damage['BBG'],max_damage['Top of curve (60%)']))
    
    max_values = overlaid.BBG.apply(lambda x : max_dam_dict[x]).values
    
    overlaid['damage_{}'.format(curve)] = (np.interp(overlaid['fg10_{}_{}'.format(model,fc_time)].values,damage_functions['Wind velocity (kph)'].values,damage_functions[curve].values)*overlaid['m2'].values*max_values)
    
    print('Total damage for {} for {} at forecast {}: {}'.format(province,model,fc_time,overlaid['damage_{}'.format(curve)].sum()))

    overlaid['province'] = province

    overlaid.geometry = pygeos.to_wkb(overlaid.geometry)

    overlaid.to_feather(damage_province_path)

    #return file
    return overlaid


if __name__ == '__main__':       

    osm_files = os.path.join("C:\\Data",'OSM','feather')
    provinces = [x.split('.')[0] for x in os.listdir(osm_files)]

    #damage_estimate('friesland')

    # # in parallel get exposure information
    with Pool(processes=cpu_count()) as pool: 
         pool.map(damage_estimate,provinces,chunksize=1)      
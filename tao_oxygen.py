import os
import pyproj
import rasterio
import numpy as np
import pandas as pd
import rasterio.mask
from tqdm import tqdm
import geopandas as gpd
from shapely.ops import transform


def main(image_path, in_fp_grid, out_dir_shp):
    df = gpd.read_file(in_fp_grid)
    df = df.to_crs(32617)
    data_geo = []
    name_classes = {0:'Backgroud',
                    1:'Mangrove/Wetlands',
                    2:'Forest',
                    3:'Build up/Human settlement',
                    4:'Water',
                    5:'Grass/Scrobland'}

    n = 1
    for i in tqdm(df.iterrows()):
        shape_leng = i[1][0]
        shape_area = i[1][1]
        num_celda = i[1][2]
        geometry = i[1][-1] 
        with rasterio.open(image_path,BIGTIFF='YES') as src:
            out_image, _ = rasterio.mask.mask(src, [geometry], crop=True)
        value,count = np.unique(out_image, return_counts=True)
        if 0 in value:
            if len(value)>1:
                index = 0
                value = np.delete(value, index)
                count = np.delete(count, index)
        max_count_id = np.argmax(count)
        max_value_id = value[max_count_id]
        classes = name_classes[max_value_id]
        total_pixel = np.sum(count)
        dict_obj = dict(zip(value, count))
        
        statics = {}
        for j in dict_obj:
            percent = dict_obj[j]/total_pixel
            statics.update({name_classes[j]:str(percent*100)+'%'})
        country = 'Panama'
        location = 'Punta Patino,Darien'
        protected_area_name = 'Reserva Ancon'
        # ecosystem = "(Forest, Mangrove/wetlands, Grass/stubble, Water, Build up/Human settlement)"
        year = os.path.basename(image_path).split('_')[0]
        quarter = os.path.basename(image_path).split('_')[1]
        coverage_date = year+'_'+quarter
        ids = 'PAN-'+ ('000000000'+ str(n))
        n+=1
        
        wgs84 = pyproj.CRS('EPSG:32617')
        utm = pyproj.CRS('EPSG:4326')
        project = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform
        centroid = transform(project, geometry.centroid)
        
        latitud = centroid.y
        longitud = centroid.x
        data_geo.append([ids, shape_leng, shape_area, num_celda, country,
                        location, protected_area_name, latitud, longitud,
                        classes, geometry, statics, coverage_date])
    new_df = pd.DataFrame(data_geo, columns=['ID', 'Shape_Leng', 'Shape_Area', 'Num_Celda', 'Country', 
                                            'Location', 'Protected_area_name', 'Latitud', 'Longitud', 
                                            'Ecosystem', 'geometry', 'properties', 'Date'])

    gdf = gpd.GeoDataFrame(new_df, crs=df.crs)
    os.makedirs(out_dir_shp, exist_ok=True)

    gdf = gdf.to_crs(32617)
    data_geo = []
    dict_class = {'Backgroud':0,
                'Mangrove/Wetlands':27.28,
                'Forest':27.28,
                'Build up/Human settlement':0,
                'Water':1.51,
                'Grass/Scrobland':11.84}

    for i in gdf.iterrows():
        total_oxy = 0 
        for j in i[1]['properties']:
            total_oxy+= (dict_class[j] * float(i[1]['properties'][j].split('%')[0])) 
        data_geo.append([i[1]['ID'], i[1]['Shape_Leng'], i[1]['Shape_Area'], i[1]['Num_Celda'],
                        i[1]['Country'], i[1]['Location'], i[1]['Protected_area_name'], i[1]['Latitud'], 
                        i[1]['Longitud'], i[1]['Ecosystem'], i[1]['properties'], i[1]['Date'], i[1]['geometry'], total_oxy])

    new_df = pd.DataFrame(data_geo, columns=['ID', 'Shape_Leng', 'Shape_Area', 'Num_Celda', 'Country', 
                                            'Location', 'Protected_area_name', 'Latitud', 'Longitud', 
                                            'Ecosystem', 'properties', 'date', 'geometry','Oxygen'])
    gdf = gpd.GeoDataFrame(new_df, crs=gdf.crs)
    out_fp_shp = os.path.join(out_dir_shp, os.path.basename(image_path).replace('.tif','_final.geojson'))
    gdf.to_file(out_fp_shp)

if __name__=='__main__':
    image_path = r"D:\HIEUUUUUU\TMP\run_panama\origin\RS_1.tif"
    in_fp_grid = r"D:\HIEUUUUUU\TMP\run_panama\origin\AOI_30k\AOI_30k.shp"
    out_dir_shp = r'D:\HIEUUUUUU\TMP\run_panama\origin\Phase2_tmp'
    main(image_path, in_fp_grid, out_dir_shp)
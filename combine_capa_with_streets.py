import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask
import json
from shapely.geometry import shape
import numpy as np
import os

def get_centroid(geom_obj):
    """
    Return a Shapely centroid for any GeoJSON geometry.
    """
    g = shape(geom_obj)          # convert to Shapely
    return g.centroid            # works for all geometry types

def find_nearest_temp_to_geometry(geom_obj, points):
    """
    Given one street geometry (any type) and a list of
    (x, y, temp) tuples, return the temperature of
    the closest traverse point.
    """
    centroid = get_centroid(geom_obj)   # always safe now

    # squared-distance search
    min_dist = float("inf")
    nearest_temp = None
    for x, y, temp in points:
        dist = (centroid.x - x)**2 + (centroid.y - y)**2
        if dist < min_dist:
            min_dist = dist
            nearest_temp = temp
    return nearest_temp

def extract_temp_from_raster_buffered(geometry, raster_path, buffer_distance=15):
    """
    Extract temperature from raster using a buffered geometry for small features
    buffer_distance: buffer in meters (default 15m = 1.5 pixels at 10m resolution)
    """
    if not os.path.exists(raster_path):
        print(f"Raster file not found: {raster_path}")
        return None
        
    try:
        with rasterio.open(raster_path) as src:
            # Convert geometry to shapely
            geom_shape = shape(geometry)
            
            # Transform geometry to raster CRS if needed
            if str(src.crs) != 'EPSG:4326':
                try:
                    transformed_geom = transform_geom('EPSG:4326', src.crs, geometry)
                    geom_shape = shape(transformed_geom)
                except Exception as transform_error:
                    print(f"CRS transformation failed: {transform_error}")
                    return None
            
            # Buffer the geometry to ensure it covers at least a few pixels
            buffered_geom = geom_shape.buffer(buffer_distance)
            
            # Check if buffered geometry intersects with raster bounds
            raster_bounds = src.bounds
            raster_polygon = shape({
                'type': 'Polygon',
                'coordinates': [[
                    [raster_bounds.left, raster_bounds.bottom],
                    [raster_bounds.right, raster_bounds.bottom], 
                    [raster_bounds.right, raster_bounds.top],
                    [raster_bounds.left, raster_bounds.top],
                    [raster_bounds.left, raster_bounds.bottom]
                ]]
            })
            
            if not buffered_geom.intersects(raster_polygon):
                print(f"Buffered geometry outside raster bounds")
                return None
            
            # Create mask for the buffered geometry
            mask = geometry_mask([buffered_geom], 
                               transform=src.transform, 
                               invert=True, 
                               out_shape=src.shape)
            
            # Check if mask has any True values
            if not np.any(mask):
                print(f"Buffered geometry mask is still empty")
                return None
            
            # Read raster data
            raster_data = src.read(1)
            
            # Extract values within buffered geometry
            values = raster_data[mask]
            
            # Filter out invalid values
            valid_values = values[~np.isnan(values)]
            if src.nodata is not None:
                valid_values = valid_values[valid_values != src.nodata]
            
            if len(valid_values) > 0:
                return np.mean(valid_values)
            else:
                print(f"No valid temperature values found within buffered geometry")
                return None
                
    except Exception as e:
        print(f"Error processing raster {raster_path}: {e}")
        return None

def tag_street_temperatures(streets_geojson, raster_paths, traverse_files):
    """
    Tag streets with temperature data for morning, afternoon, and evening
    """
    
    # Load street data
    if streets_geojson['type'] == 'Feature':
        features = [streets_geojson]
    else:
        features = streets_geojson['features']
    
    time_periods = ['morning', 'afternoon', 'evening']
    
    processed_count = 0
    error_count = 0
    
    for feature in features:
        processed_count += 1
        if processed_count % 100 == 0:
            print(f"Processed {processed_count} features...")
            
        # Ensure properties dictionary exists
        if 'properties' not in feature:
            feature['properties'] = {}
            
        for period in time_periods:
            temp_property = f'temp_{period}'
            temp = None
            
            if period in raster_paths and raster_paths[period]:
                temp = extract_temp_from_raster_buffered(feature['geometry'], raster_paths[period])
                if temp is not None:
                    feature['properties'][temp_property] = round(float(temp), 1)
                    print(f"Successfully extracted {period} temp from raster: {temp}")
                    continue
                else:
                    print(f"Raster extraction failed for {period}, trying traverse...")

            # Fallback to traverse points if raster extraction failed OR wasn't available
            if temp is None and period in traverse_files and traverse_files[period]:
                try:
                    with open(traverse_files[period], 'r') as f:
                        traverse_data = json.load(f)
                    
                    # Extract traverse points
                    traverse_points = []
                    for trav_feature in traverse_data['features']:
                        coords = trav_feature['geometry']['coordinates']
                        temp_val = trav_feature['properties'].get('t_f', None)
                        if temp_val is not None:
                            traverse_points.append((coords[0], coords[1], temp_val))
                    
                    # Find nearest temperature
                    if traverse_points:
                        temp = find_nearest_temp_to_geometry(feature['geometry'], traverse_points)
                        if temp is not None:
                            feature['properties'][temp_property] = round(float(temp), 1)
                            print(f"Successfully extracted {period} temp from traverse: {temp}")
                        else:
                            print(f"No valid temperature found for {period}")
                    else:
                        print(f"No valid traverse points for {period}")
                            
                except Exception as e:
                    print(f"Error processing traverse file for {period}: {e}")
                    error_count += 1
    
    print(f"Processing complete. {processed_count} features processed, {error_count} errors encountered.")
    return streets_geojson

def load_geojson_data():
    """Load the campus streets and city streets geojson data"""
    
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Load whole city streets
    city_streets_path = os.path.join(script_dir, 'data/Bloomington_Road_Centerlines.geojson')
    try:
        with open(city_streets_path, 'r') as f:
            whole_city = json.load(f)
        print(f"Loaded {len(whole_city['features'])} city street features")
    except FileNotFoundError:
        print(f"City streets file not found: {city_streets_path}")
        whole_city = {"type": "FeatureCollection", "features": []}
    
    # Load campus/walking paths
    campus_streets_path = os.path.join(script_dir, 'data/Walks.geojson')
    try:
        with open(campus_streets_path, 'r') as f:
            campus_streets = json.load(f)
        print(f"Loaded {len(campus_streets['features'])} campus/walk features")
    except FileNotFoundError:
        print(f"Campus streets file not found: {campus_streets_path}")
        campus_streets = {"type": "FeatureCollection", "features": []}
    
    return campus_streets, whole_city

def save_tagged_data(tagged_campus, tagged_city, script_dir):
    """Save the tagged street data to new geojson files"""
    
    campus_output_path = os.path.join(script_dir, 'tagged_campus_streets.geojson')
    city_output_path = os.path.join(script_dir, 'tagged_city_streets.geojson')
    
    with open(campus_output_path, 'w') as f:
        json.dump(tagged_campus, f, indent=2)
    print(f"Saved tagged campus streets to: {campus_output_path}")
    
    with open(city_output_path, 'w') as f:
        json.dump(tagged_city, f, indent=2)
    print(f"Saved tagged city streets to: {city_output_path}")

campus_streets, whole_city = load_geojson_data()
script_dir = os.path.dirname(os.path.abspath(__file__))

raster_paths = {
    'morning': os.path.join(script_dir, 'data/Heat_Watch_Campaign/rasters_chw_bloomington_indiana_081424/bloomington-indiana_am_temp_f.tif'),
    'afternoon': os.path.join(script_dir, 'data/Heat_Watch_Campaign/rasters_chw_bloomington_indiana_081424/bloomington-indiana_af_temp_f.tif'),
    'evening': os.path.join(script_dir, 'data/Heat_Watch_Campaign/rasters_chw_bloomington_indiana_081424/bloomington-indiana_pm_temp_f.tif'),
}

traverse_files = {
    'morning': os.path.join(script_dir, 'data/Heat_Watch_Campaign/traverses_chw_bloomington_indiana_081424/bloomington-indiana_am_trav.geojson'),
    'afternoon': os.path.join(script_dir, 'data/Heat_Watch_Campaign/traverses_chw_bloomington_indiana_081424/bloomington-indiana_af_trav.geojson'),
    'evening': os.path.join(script_dir, 'data/Heat_Watch_Campaign/traverses_chw_bloomington_indiana_081424/bloomington-indiana_pm_trav.geojson'),
}

tagged_campus_streets = tag_street_temperatures(campus_streets, raster_paths, traverse_files)
tagged_city_streets = tag_street_temperatures(whole_city, raster_paths, traverse_files)

# Save the results
save_tagged_data(tagged_campus_streets, tagged_city_streets, script_dir)

print("Temperature tagging complete!")
print(f"Campus streets: {len(tagged_campus_streets['features'])} features")
print(f"City streets: {len(tagged_city_streets['features'])} features")


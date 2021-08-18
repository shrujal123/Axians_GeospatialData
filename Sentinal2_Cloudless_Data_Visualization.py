!pip install -q geopandas geemap

import ee
import geemap
import geopandas as gpd
import pandas as pd
import datetime

ee.Authenticate()
ee.Initialize()

saoMiguel = ee.Geometry.Point([-25.3425, 37.7532])
boundary = geemap.kml_to_ee('./study_area.kml').geometry()


def add_cloud_bands(img):
    # Get s2cloudless image, subset the probability band.
    cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

    # Condition s2cloudless by the probability threshold value.
    is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')

    # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))

def add_shadow_bands(img):
    # Identify water pixels from the SCL band.
    not_water = img.select('SCL').neq(6)

    # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    SR_BAND_SCALE = 1e4
    dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

    # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject(**{'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    # Identify the intersection of dark pixels with cloud shadow projection.
    shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    # Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

def add_cld_shdw_mask(img):
    # Add cloud component bands.
    img_cloud = add_cloud_bands(img)

    # Add cloud shadow component bands.
    img_cloud_shadow = add_shadow_bands(img_cloud)

    # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

    # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    is_cld_shdw = (is_cld_shdw.focal_min(2).focal_max(BUFFER*2/20)
        .reproject(**{'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'))

    # Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)

def apply_cld_shdw_mask(img):
    # Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
    not_cld_shdw = img.select('cloudmask').Not()

    # Subset reflectance bands and update their masks, return the result.
    return img.select('B.*').updateMask(not_cld_shdw)

  
AOI = saoMiguel
START_DATE = '2018-01-01'
END_DATE = '2021-07-30'
CLOUD_FILTER = 70
CLD_PRB_THRESH = 60
NIR_DRK_THRESH = 0.15
CLD_PRJ_DIST = 2
BUFFER = 10


def get_s2_sr_cld_col(aoi, start_date, end_date):
    # Import and filter S2 SR.
    s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_SR')
        .filterBounds(aoi)
        .filterDate(start_date, end_date))

    # Import and filter s2cloudless.
    s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        .filterBounds(aoi)
        .filterDate(start_date, end_date))

    # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': s2_sr_col,
        'secondary': s2_cloudless_col,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))

sen2 = get_s2_sr_cld_col(AOI, START_DATE, END_DATE)

sen2 = sen2.map(add_cld_shdw_mask).map(apply_cld_shdw_mask)


def visualize(collection, dates, time_delta, ndvi=False, ndwi=False):
  # Create the map
  Map = geemap.Map(center=[37.7529, -25.3453], zoom=13)
  Map.add_kml('study_area.kml',layer_name='Study Area')
  # Add the images as layers
  delta = pd.Timedelta(days=time_delta)
  for date in dates:
    print(f"\rProcessing {date.date()}", end='')
    image = ee.Image(collection
                    .filterDate(date, date + delta)
                    .median()
                    .clip(boundary)
                    )
    layer_name = str(date).split(' ')[0]
    Map.addLayer(image, trueColour, layer_name, shown=False)
    if ndvi:
      Map.addLayer(image, ndvi_params, layer_name + ' NDVI', shown=False)
    if ndwi:
      Map.addLayer(image, ndwi_params, layer_name + ' NDWI', shown=False)

  Map.addLayerControl()
  return Map


dates = pd.read_excel('./proposed_dataset_list.xlsx', usecols=['DATE'], parse_dates=True)['DATE'].sort_values()
print(f"Start: {dates.iloc[0].date()}, End: {dates.iloc[-1].date()}")
dates.iloc[-1] - dates.iloc[0]


Map = geemap.Map(center=[37.7529, -25.3453], zoom=13)
Map.add_kml('study_area.kml',layer_name='Study Area')

delta = pd.Timedelta(days=1)
for date in dates:
  print(f"\rProcessing {date.date()}", end='')
  image = ee.Image(sen2
                   .filterDate(date - delta, date + delta)
                   .sort("CLOUD_COVERAGE_ASSESSMENT")
                   .first()
                   .clip(boundary)
                   )
  layer_name = str(date).split(' ')[0]
  Map.addLayer(image, trueColour, layer_name, shown=False)

Map.addLayerControl()
Map



dates_BW = pd.date_range("2018-01-01", "2021-08-01", freq="SM").to_list()
len(dates_BW)
visualize(sen2, dates_BW, time_delta=15)


dates_M = pd.date_range("2018-01-01", "2021-08-01", freq="MS").to_list()
len(dates_M)
visualize(sen2, dates_M, time_delta=30)


dates_Q = pd.date_range("2018-01-01", "2021-08-01", freq="QS").to_list()
len(dates_Q)
visualize(sen2, dates_Q, time_delta=90)

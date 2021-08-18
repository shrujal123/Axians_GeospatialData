[a link](https://github.com/shrujal123/geoinfra/blob/main/Sentinal2_Cloudless_Data_Visualization.py)

collection = ee.ImageCollection(sen2).filterDate("2021-07-15", "2021-08-01")
point = {'type':'Point', 'coordinates':[-25.3511, 37.7623]};
#poly_geometry = geemap.geojson_to_ee('/content/drive/MyDrive/RoW_powerline_corridors_webm.geojson')
poly_geometry = geemap.shp_to_ee('./RoW_powerline_corridors.shp')

info = collection.getRegion(poly_geometry, 50).getInfo()

header = info[0]
data = np.array(info[1:])

iTime = header.index('time')
time = [datetime.datetime.fromtimestamp(i/1000) for i in (data[0:,iTime].astype(int))]

band_list = ['B1',u'B2']

iBands = [header.index(b) for b in band_list]
yData = data[0:,iBands].astype(np.float)

#Calculate NDVI
red = yData[:,0]
nir = yData[:,1]
ndvi = (nir - red) / (nir + red)

df = pd.DataFrame(data = ndvi, index =list(range(len(ndvi))), columns = ['NDVI'])
df = df.interpolate()
df['Date'] = pd.Series(time, index = df.index)
df = df.set_index(df.Date)
df.index = pd.to_datetime(df.index)
df['NDVI'] = df['NDVI'].fillna(0)

df.info()

sns.set(rc={'figure.figsize':(15, 6)})
df['NDVI'].plot(linewidth = 0.5)


sd = seasonal_decompose(df['NDVI'], model='additive', freq=352)
sd.seasonal.plot()
sd.trend.plot()
sd.resid.plot()
plt.legend(['Seasonality', 'Trend', 'Residuals'])

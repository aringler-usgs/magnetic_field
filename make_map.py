#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from obspy.core import UTCDateTime
import utils
import pickle
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.font_manager
from obspy.clients.fdsn import Client

path = '/home/aringler/inst_test/magseis/noise_data/'
net ='*'
chan ='LHZ'
client = Client('IRIS')
stime = UTCDateTime('2019-241T00:00:00')
etime = UTCDateTime('2019-242T00:00:00')
permin = 40.
permax = 400.



def setupmap(central_lon, central_lat,handle):
    #handle = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    handle.set_extent(extent)

    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN)
    handle.add_feature(cfeature.COASTLINE)
    handle.add_feature(cfeature.BORDERS, linestyle=':')
    handle.add_feature(cfeature.LAKES, alpha=0.5)
    handle.add_feature(cfeature.RIVERS)
    handle.add_feature(cfeature.STATES, edgecolor='gray')
    return handle

inv = utils.get_dataless(net, chan, stime, etime, client)

# lats, lons,noise = [], [], []
# for net in inv:
    # for sta in net:
        # for chans in sta:
            # print(path)
            # print(sta.code)
            # try:
            # #if True:
                # per, psd = utils.get_psd_stats(path, sta.code, stime, etime, 10.)
                # print(per)
                # print(psd)
                
                # psdm = psd[(per >= permin) & (per <= permax)]
                # noise.append(np.mean(psdm))
            
                # lats.append(chans.latitude)
                # lons.append(chans.longitude)
            # except:
                # print('Problem with: ' + sta.code)
                

# f = open('Test_data_241','wb')
# pickle.dump([lats, lons, noise], f)
# f.close()

f = open('Test_data_241','rb')
lats, lons, noise = pickle.load(f)

f = open('Test_data_244','rb')
lats2, lons2, noise2 = pickle.load(f)


minval = np.mean(noise) -2*np.std(noise)
maxval = np.mean(noise) + 2*np.std(noise)

lats = np.array(lats)
lons = np.array(lons)
noise = np.array(noise)
noise2 = np.array(noise2)

diff = noise2 - noise


print(len(lats))
print(len(lons))
print(len(noise))

fig= plt.figure(figsize=(12,12))


boxcoords=[min(lats) -1., min(lons)-1., max(lats) +1. , max(lons) + 1.]
extent=[boxcoords[1], boxcoords[3], boxcoords[0], boxcoords[2]]
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])            
ax = plt.subplot(1,1,1, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
sc = ax.scatter(lons, lats, c = noise2, transform=ccrs.PlateCarree(), vmin=minval, vmax=maxval)
cbar = fig.colorbar(sc, orientation='vertical', shrink=0.7)
plt.savefig('Alaskamap_244.png', format='PNG')
plt.show()

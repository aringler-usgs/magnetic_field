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

fig= plt.figure(figsize=(12,16))


boxcoords=[min(lats) -1., min(lons)-1., max(lats) +1. , max(lons) + 1.]
extent=[boxcoords[1], boxcoords[3], boxcoords[0], boxcoords[2]]
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])            
ax = plt.subplot(3,2,1, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
sc = ax.scatter(lons, lats, c = noise, transform=ccrs.PlateCarree(), vmin=-185, vmax=-175)
cbar = fig.colorbar(sc, orientation='horizontal', shrink=0.7)
cbar.ax.set_xlabel('Power (dB)', fontsize=18)
plt.text(-0.1, 1., '(a)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
ax = plt.subplot(3,2,3, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
sc = ax.scatter(lons, lats, c = noise2, transform=ccrs.PlateCarree(), vmin=-185, vmax=-175)
cbar = fig.colorbar(sc, orientation='horizontal', shrink=0.7)
cbar.ax.set_xlabel('Power (dB)', fontsize=18)
plt.text(-0.1, 1., '(c)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
ax = plt.subplot(3,2,5, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
sc = ax.scatter(lons, lats, c = diff, transform=ccrs.PlateCarree(), vmin=0, vmax=10)
cbar = fig.colorbar(sc, orientation='horizontal', shrink=0.7)
cbar.ax.set_xlabel('Difference (dB)', fontsize=18)
plt.text(-0.1, 1., '(e)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
stime = UTCDateTime('2019-241T00:00:00')
st = client.get_waveforms("IU", "COLA", "*",
                           "LF*", stime, stime + 24*60*60, attach_response = True)

print(st)
st.detrend('constant')
st.merge(fill_value=0)
ax = plt.subplot(3,2,2)
for tr in st:
	t=tr.times()/(24*60*60)
	ax.plot(tr.times()/(24*60*60), np.sqrt((tr.data/(4.1943*10))**2), label=tr.id, alpha=0.5)
plt.ylim((0., 800.))
plt.ylabel('RMS (nT)')
plt.xlim((min(t), max(t)))
plt.xlabel('Time (hr)')
plt.text(-0.15, 1., '(b)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
stime = UTCDateTime('2019-244T00:00:00')
st = client.get_waveforms("IU", "COLA", "*",
                           "LF*", stime, stime + 24*60*60, attach_response = True)

print(st)
st.detrend('constant')
st.merge(fill_value=0)
ax = plt.subplot(3,2,4)
for tr in st:
	t=tr.times()/(24*60*60)
	ax.plot(tr.times()/(24*60*60), np.sqrt((tr.data/(4.1943*10))**2), label=tr.id.replace('.', ' '), alpha=0.5)
plt.ylim((0., 800.))
plt.ylabel('RMS (nT)')
plt.xlim((min(t), max(t)))
plt.xlabel('Time (hr)')
plt.text(-0.15, 1., '(d)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
          fancybox=True, shadow=True, ncol=1, fontsize=18)









plt.savefig('Alaskamap.png', format='PNG', dpi=200)
plt.savefig('Alaskamap.pdf', format='PDF', dpi=400)
plt.show()

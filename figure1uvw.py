#!/usr/bin/env python
from obspy.core import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from scipy.optimize import fmin
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

# This code generates figure 1 of the Ringler et al. magnetic field paper but with the UVW transfer

stime = UTCDateTime('2019-243T00:00:00')
etime = stime + 1*24*60*60

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

# X is E
def rotate_xyz_uvw(data1, data2, data3):
    # create new trace objects with same info as previous
    rotatedU = data1.copy()
    rotatedV = data2.copy()
    rotatedW = data3.copy()
    rotatedU.data = -np.sqrt(2./3.)*data3.data + np.sqrt(1./3.)*data1.data
    rotatedV.data = np.sqrt(1./6.)*data3.data + np.sqrt(1./2.)*data2.data + np.sqrt(1./3.)*data1.data
    rotatedW.data = np.sqrt(1./6.)*data3.data - np.sqrt(1./2.)*data2.data + np.sqrt(1./3.)*data1.data
    st2 = Stream()
    st2=Stream(traces=[rotatedU, rotatedV, rotatedW])
    return st2

def rotate_uvw_xyz(data1, data2, data3):
    # create new trace objects with same info as previous
    # from U, V, W
    rotatedZ = data1.copy()
    rotatedN = data2.copy()
    rotatedE = data3.copy()
    rotatedE.data = - np.sqrt(2./3.)*data1.data + np.sqrt(1./6.)*data2.data + np.sqrt(1./6.)*data3.data
    rotatedN.data = np.sqrt(1./2.)*data2.data - np.sqrt(1./2.)*data3.data
    rotatedZ.data = np.sqrt(1./3.)*data1.data + np.sqrt(1./3.)*data2.data + np.sqrt(1./3.)*data3.data
    st2 = Stream()
    st2=Stream(traces=[rotatedZ, rotatedN, rotatedE])
    return st2

net, sta, loc = 'IU', 'COLA', '00'

pmax = 1./800.
pmin = 1./40.

client = Client()
st = client.get_waveforms(net, sta, loc,
                          "LH*", stime, etime, attach_response = True)


# Magnetic field data is in T
inv = client.get_stations(network=net, station=sta, location=loc,
						 starttime=stime,endtime=etime, level="response")

st.remove_response(inventory=inv, output='ACC')
# Rotate to uvw
st2 = rotate_xyz_uvw(st[1], st[2], st[0])
for idx, ch in enumerate(['LHU', 'LHV', 'LHW']):
    st2[idx].stats.channel = ch

st += st2
# When using some TA sites you need to chop off a sample
#for tr in st:
#    tr.data = tr.data[1:]

st += client.get_waveforms("IU", "COLA", "*",
                           "LF*", stime, etime, attach_response = True)

st.detrend('constant')
st.merge(fill_value=0)
for tr in st.select(location='40'):
    tr.data /= 41.943
    tr.data /= 10**9


# Filter and then estimate coefficents
st.filter('bandpass', freqmin=pmax, freqmax=pmin)
st.taper(0.05)

for chan in ['LHU', 'LHV','LHW', 'LHZ', 'LH1', 'LH2']:
    def resi_fun(x):
        data = st.select(location=loc, channel=chan)[0].data.copy()
        for idx, tr in enumerate(st.select(location='40')):
            data -= x[idx]*tr.data.copy()
        return np.sum(data**2)/len(data)

    sol = fmin(resi_fun, [1., 1., 1.])

    def correct(x):
        data = st.select(location=loc, channel=chan)[0].data.copy()
        for idx, tr in enumerate(st.select(location='40')):
            data -= x[idx]*tr.data.copy()
        return data
    print(chan)
    print('Here are the coefficients: ' + str(sol))
    print('Here is the initial residual: ' + str(resi_fun([0., 0., 0.])))
    print('Here is the new residual: ' + str(resi_fun(sol)))

    reduction = 100.*(resi_fun(sol) - resi_fun([0., 0., 0.]))/resi_fun([0., 0., 0.])
    print(chan + ' variance reduction: ' + str(reduction))
    print('Here is the magnitude: ' + str(np.sqrt(np.sum((sol)**2))))

#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from scipy.optimize import fmin
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
stime = UTCDateTime('2019-243T00:00:00')
etime = stime + 1*24*60*60

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
#mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

chans = ['LHZ']
locs = ['10', '00']
fig = plt.figure(2, figsize=(14,10))
for idx, loc in enumerate(locs):
    plt.subplot(3,1,idx+1)
    for chan in chans:

        #net, sta, loc, chan = 'TA', 'H22K', '*', 'LHZ'
        net, sta = 'IU', 'COLA'

        pmax = 1./800.
        pmin = 1./40.

        client = Client()
        st = client.get_waveforms(net, sta, loc,
                                   chan, stime, etime, attach_response = True)


        stR = st.copy()
        stR.detrend('constant')
        st += client.get_waveforms("IU", "COLA", "*",
                                   "LF*", stime, etime, attach_response = True)

        st.detrend('constant')
        st.merge(fill_value=0)
        st2 = st.select(location='40').copy()
        for tr in st2:
            tr.data /= 41.943
        paz = {'poles': [-0.037+0.037j, -0.037-0.037j], 'zeros': [0j], 'gain':1., 'sensitivity': 1.}

        # The magnetic field is now interms of m/s/s/nT normalized
        for tr in st.select(location='40'):
            tr.simulate(paz_simulate=paz)
        st2 = st.select(location='40').copy()

        # Filter and then estimate coefficents
        st.filter('bandpass', freqmin=pmax, freqmax=pmin)
        #st.normalize()
        st2.filter('bandpass', freqmin=pmax, freqmax=pmin)

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

        print(sol)
        print(resi_fun([0., 0., 0.]))
        print(resi_fun(sol))

        reduction = 100.*(resi_fun(sol) - resi_fun([0., 0., 0.]))/resi_fun([0., 0., 0.])
        print(chan + ' variance reduction')
        print(reduction)



        tr = st.select(location=loc, channel=chan)[0]
        times = tr.times()/(60*60)
        plt.subplot(3,1,idx + 1)
        if idx == 0:
            plt.text(-3, 400, '(a)')
        else:
            plt.text(-3, 400, '(b)')
        plt.plot(times, tr.data,label = (tr.id.replace('IU.COLA','')).replace('.',' '), alpha=0.5)
        plt.plot(times, correct(sol),label = 'Magnetic Field Corrected', alpha=0.5)
        plt.xlim((min(times), max(times)))
        plt.ylim((-400., 400.))
        plt.legend(loc=8, ncol=2)
        plt.ylabel('Counts')
plt.subplot(3,1,3)
st2.normalize()
for tr in st2.select(location='40'):
    plt.plot(times, tr.data,label = tr.id.replace('IU.COLA.40.',''), alpha=0.5)
    plt.xlim((min(times), max(times)))
    plt.legend(loc=8, ncol=3)
    plt.text(-3, 1, '(c)')
    plt.ylabel('Scaled (T)')
plt.xlabel('Time (Hr)')
plt.savefig('Timeseries' + chan + '.png', format='PNG')
plt.show()







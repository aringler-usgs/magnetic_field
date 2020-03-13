#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from scipy.optimize import fmin
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

# This code generates figure 1 of the Ringler et al. magnetic field paper

stime = UTCDateTime('2019-243T00:00:00')
etime = stime + 1*24*60*60

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

chans = ['LHZ']
locs = ['10', '00']
fig = plt.figure(2, figsize=(14,10))
for idx, loc in enumerate(locs):
    plt.subplot(3,1,idx+1)
    for chan in chans:

        #net, sta, loc, chan = 'TA', 'B20K', '*', 'LHZ'
        net, sta = 'IU', 'COLA'

        pmax = 1./800.
        pmin = 1./40.

        client = Client()
        st = client.get_waveforms(net, sta, loc,
                                   chan, stime, etime, attach_response = True)

        # When using some TA sites you need to chop off a sample
        #for tr in st:
        #    tr.data = tr.data[1:]
        stR = st.copy()
        stR.detrend('constant')
        st += client.get_waveforms("IU", "COLA", "*",
                                   "LF*", stime, etime, attach_response = True)

        
        st.detrend('constant')
        st.merge(fill_value=0)
        #st2 = st.select(location='40').copy()
        for tr in st.select(location='40'):
            tr.data /= 41.943
            tr.data /= 10**9

        # Magnetic field data is in T
        inv = client.get_stations(network=net, station=sta, location=loc,
                            channel=chan, starttime=stime,endtime=etime, level="response")

        #resp = inv.get_response(net +"." + sta + "." + loc + "." + chan, stime )
        #resp_paz = resp.get_paz()
        #paz = {'poles':resp_paz.poles, 'zeros':resp_paz.zeros[1:], 'gain': resp_paz.normalization_factor, 'sensitivity': resp_paz*(2**26/40.)}

        # The magnetic field is now interms of m/s/s/nT normalized
        #for tr in st.select(location='40'):
        #    tr.simulate(paz_simulate=paz)
        #st2 = st.select(location='40').copy()
        st.select(channel=chan).remove_response(inventory=inv, output='ACC')
        # Filter and then estimate coefficents
        st.filter('bandpass', freqmin=pmax, freqmax=pmin)
        st.taper(0.05)
        #st2.filter('bandpass', freqmin=pmax, freqmax=pmin)
        print(st.select(location='40'))
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

        print('Here are the coefficients: ' + str(sol))
        print('Here is the initial residual: ' + str(resi_fun([0., 0., 0.])))
        print('Here is the new residual: ' + str(resi_fun(sol)))

        reduction = 100.*(resi_fun(sol) - resi_fun([0., 0., 0.]))/resi_fun([0., 0., 0.])
        print(chan + ' variance reduction: ' + str(reduction))
        print('Here is the magnitude: ' + str(np.sqrt(np.sum((sol)**2))))


        tr = st.select(location=loc, channel=chan)[0]
        times = tr.times()/(60*60)
        plt.subplot(3,1,idx + 1)
        if idx == 0:
            plt.text(-3.5, 15, '(a)')
        else:
            plt.text(-3.5, 15, '(b)')
        plt.plot(times, tr.data*10**9,label = (tr.id.replace('IU.COLA','')).replace('.',' '), alpha=0.5)
        plt.plot(times, correct(sol)*10**9,label = 'Magnetic Field Corrected', alpha=0.5)
        plt.xlim((min(times), max(times)))
        plt.ylim((-15., 15.))
        plt.legend(loc=8, ncol=2)
        plt.ylabel('Acc. $(nm/s^2)$')
plt.subplot(3,1,3)
#st2.normalize()
# Here we plot the magnetic field data
for tr in st.select(location='40'):
    plt.plot(times, tr.data*10**9,label = tr.id.replace('IU.COLA.40.',''), alpha=0.5)
    plt.xlim((min(times), max(times)))
    plt.legend(loc=8, ncol=3)
    plt.ylim((-300., 300.))
    plt.text(-3.5, 300, '(c)')
    plt.ylabel('Magnetic Field ($nT$)')
plt.xlabel('Time August 31, 2019 (UTC Hr)')
plt.savefig('figure1' +  '.png', format='PNG')
#plt.show()

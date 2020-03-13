#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from scipy.optimize import fmin
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
stime = UTCDateTime('2019-271T00:00:00')
etime = stime + 1*24*60*60

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

chans = ['LH1', 'LH2', 'LHZ']
locs = ['10', '00']
fig = plt.figure(2, figsize=(14,10))
for idx, loc in enumerate(locs):
    plt.subplot(2,1,idx+1)
    for chan in chans:

        #net, sta, loc, chan = 'TA', 'H22K', '*', 'LHZ'
        net, sta = 'IU', 'COLA'

        pmax = 1./2000.
        pmin = 1./10.

        client = Client()
        st = client.get_waveforms(net, sta, loc,
                                   chan, stime, etime, attach_response = True)
        #st = read("/msd/" + net + '_' + sta + '/' + str(ctime.year) + "/" +
        #          str(ctime.julday).zfill(3) + "/"+ loc + "_" + chan + "*")

        stR = st.copy()
        stR.detrend('constant')
        st += client.get_waveforms("IU", "COLA", "*",
                                   "LF*", stime, etime, attach_response = True)
        #st += read("/msd/" + net + '_' + sta + '/' + str(ctime.year) + "/" +
        #          str(ctime.julday).zfill(3) + "/40_LF*")
        st.detrend('constant')
        st.merge(fill_value=0)
        if loc == '10':
            paz = {'poles': [-0.037+0.037j, -0.037-0.037j], 'zeros': [0j], 'gain':1., 'sensitivity': 1.}
        else:
            paz = {'poles': [-0.01234+0.01234j, -0.01234-0.01234j], 'zeros': [0j], 'gain':1., 'sensitivity': 1.}
        #inv = read_inventory("/APPS/metadata/RESPS/RESP." + net + "." + sta + "." +
        #                    "00." + chan)
        # The magnetic field is now interms of m/s/s/nT normalized
        for tr in st.select(location='40'):
            tr.simulate(paz_simulate=paz)

        # Filter and then estimate coefficents
        st.filter('bandpass', freqmin=pmax, freqmax=pmin)
        #st.normalize()


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
        times = tr.times()/(24*60*60)
        # plt.subplot(2,1,1)
        # plt.plot(times, tr.data,label = tr.id, alpha=0.5)
        # plt.plot(times, correct(sol),label = 'Magnetic Field Corrected', alpha=0.5)
        # plt.xlim((min(times), max(times)))
        # plt.legend()
        # plt.subplot(2,1,2)
        # for tr in st.select(location='40'):
            # plt.plot(times, tr.data,label = tr.id, alpha=0.5)
            # plt.xlim((min(times), max(times)))
        # plt.legend()

        # plt.savefig('Timeseries' + chan + '.png', format='PNG')
        # plt.show()
        # plt.clf()
        # plt.close()

        NFFT=2*2*4096


        def correct2(x):
            data = stR.select(location=loc, channel=chan)[0].data.copy()
            print(data)
            for idx, tr in enumerate(st.select(location='40')):
                data -= x[idx]*tr.data.copy()
            return data



        data2 = correct2(sol)

        print(data2)


        tr = stR.select(location=loc, channel=chan)[0]
        f, p = signal.welch(tr.data, fs = tr.stats.sampling_rate,
                             nperseg=NFFT, noverlap=256)
        print(tr)
        amp, f= tr.stats.response.get_evalresp_response(tr.stats.delta, NFFT,
                                                        output='ACC')

        f, pcor =  signal.welch(data2, fs = tr.stats.sampling_rate,
                             nperseg=NFFT, noverlap=256)

        # (m/s^2)^2/Hz
        p /= np.abs(amp)**2
        p=10.*np.log10(p)
        pcor /= np.abs(amp)**2
        pcor = 10.*np.log10(pcor)




        plt.semilogx(1./f, pcor, label='Corr ' + tr.stats.station + ' ' + tr.stats.location + ' ' + chan)
        plt.semilogx(1./f,p, label='Raw ' + tr.stats.station + ' ' + tr.stats.location + ' ' + chan, alpha=0.7, linestyle='dotted', linewidth=3.5)
    per, nlnm = get_nlnm()
    per, nhnm = get_nhnm()
    plt.semilogx(per, nlnm, color='k', linewidth=2)
    plt.semilogx(per, nhnm, color='k', linewidth=2, label='NLNM/NHNM')
    #plt.axvspan(1./pmin, 1./pmax, facecolor='g', alpha=0.1)
    plt.fill_between(1./f, p, pcor,color='.5')
    plt.xlabel('Period (s)')
    plt.ylabel('PSD (dB rel. 1 $(m/s^2)^2/Hz$)', fontsize=16)
    plt.xlim((2., 1000.))
    plt.ylim((-190, -90))
    plt.legend(loc=9, ncol=2, fontsize=16)
    if idx == 0:
        plt.text(1., -90, '(a)')
    else:
        plt.text(1., -90, '(b)')
plt.savefig('figure2.png', format='PNG', dpi=200)
plt.show()

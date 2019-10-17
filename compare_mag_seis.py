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

#net, sta, loc, chan = 'TA', 'H22K', '*', 'LHZ'
net, sta, loc, chan = 'IU', 'COLA', '10', 'LHZ'

pmax = 1./2000.
pmin = 1./10.

client = Client()
st = client.get_waveforms(net, sta, loc,
                           chan, stime, etime, attach_response = True)

#st[0].data = st[0].data[:-1]
stR = st.copy()
stR.detrend('constant')
st += client.get_waveforms("IU", "COLA", "*",
                           "LF*", stime, etime, attach_response = True)




print(st)
st.detrend('constant')
st.merge(fill_value=0)

paz = {'poles': [-0.037+0.037j, -0.037-0.037j], 'zeros': [0j], 'gain':1., 'sensitivity': 1.}

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

print(reduction)


fig = plt.figure(1, figsize=(12,12))
tr = st.select(location=loc, channel=chan)[0]
times = tr.times()/(24*60*60)
plt.subplot(2,1,1)
plt.plot(times, tr.data,label = tr.id, alpha=0.5)
plt.plot(times, correct(sol),label = 'Magnetic Field Corrected', alpha=0.5)
plt.xlim((min(times), max(times)))
plt.legend()
plt.subplot(2,1,2)
for tr in st.select(location='40'):
	plt.plot(times, tr.data,label = tr.id, alpha=0.5)
	plt.xlim((min(times), max(times)))
plt.legend()

plt.savefig('Timeseries.png', format='PNG')
plt.show()
plt.clf()
plt.close()

NFFT=2*2*4096
per, nlnm = get_nlnm()
per, nhnm = get_nhnm()

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

amp, f= tr.stats.response.get_evalresp_response(tr.stats.delta, NFFT,
                                                output='ACC')
                                                
f, pcor =  signal.welch(data2, fs = tr.stats.sampling_rate,
                     nperseg=NFFT, noverlap=256)                                               
                               
# (m/s^2)^2/Hz
p /= np.abs(amp)**2
p=10.*np.log10(p)
pcor /= np.abs(amp)**2
pcor = 10.*np.log10(pcor)


fig = plt.figure(2, figsize=(12,12))
plt.semilogx(f,p, label='Raw')
plt.semilogx(f, pcor, label='Corr')
plt.semilogx(1./per, nlnm, color='k', linewidth=2)
plt.semilogx(1./per, nhnm, color='k', linewidth=2, label='NLNM/NHNM')
plt.axvspan(pmin, pmax, facecolor='g', alpha=0.5)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB rel. 1 $(m/s^2)^2/Hz$)')
plt.xlim((1./(1000.), 1./2.))
plt.legend()
plt.show()




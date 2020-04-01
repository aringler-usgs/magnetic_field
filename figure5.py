#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
from obspy.signal.spectral_estimation import PPSD, get_nlnm, get_nhnm
from obspy.signal.invsim import paz_to_freq_resp
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font',size=20)
npzs = glob.glob('*.npz')

# Code for figure 5 in the manuscript.  This does the plotting you need to grab the data and do the PDF estimates.  See the PDF code.

# sensitivity is nT to V/m/s/s
#[-2.26845802e-12 -3.99273739e-13  4.83281356e-15]


def con_pow(pow, a):
    newpow = pow
    newpow = 10**(newpow/20.)
    newpow *= a
    newpow /= 10**9
    return 20.*np.log10(newpow)











per, NLNM = get_nlnm()
per, NHNM = get_nhnm()


fig = plt.figure(1,figsize=(10,14))

for idx, npz in enumerate(npzs):
    if 'QSPA' in npz:
        continue

    plt.subplot(2,1,1)
    ppsd = PPSD.load_npz(npz)

    per10, pow10 = ppsd.get_percentile(10.)
    per50, pow50 = ppsd.get_percentile(90.)
    if 'mvec' not in vars():
        mvec = pow10
        Mvec = pow50
    else:
        mvec = np.minimum(pow10,mvec)
        Mvec = np.maximum(pow50, Mvec)
    if idx == 0:
        plt.semilogx(per10, pow10, color='C0', alpha=0.5, label='10th Percentile')
        plt.semilogx(per50, pow50, color='C1', alpha=0.5, label='90th Percentile')
    else:
        plt.semilogx(per10, pow10, color='C0', alpha=0.5)
        plt.semilogx(per50, pow50, color='C1', alpha=0.5)
    plt.xlabel('Period (s)')
    plt.xlim((2.,1000.))


plt.subplot(2,1,1)
plt.text(0.8, 65, '(a)')
plt.semilogx(per10, mvec, linewidth=2, color='k', label='MIN/MAX')
plt.semilogx(per50, Mvec, linewidth=2, color='k')
plt.ylim((-45, 60))
plt.ylabel('Power (dB rel. 1 $(nT/s)^2/Hz$)', fontsize=16)
plt.legend()
plt.subplot(2,1,2)


vals = [0.01, 0.1, 1., 2.]


norm = mpl.colors.LogNorm(vmin=min(vals), vmax=max(vals))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.viridis_r)
cmap.set_array([])
plt.subplot(2,1,2)
for i, ele in enumerate(vals):
    npm = con_pow(mvec, ele)
    npM = con_pow(Mvec, ele)
    plt.semilogx(per10, npm, color=cmap.to_rgba(ele))
    plt.semilogx(per50, npM, color= cmap.to_rgba(ele))



plt.semilogx(per, NLNM, linewidth=2, color='k')
plt.semilogx(per,NHNM, linewidth=2, color='k', label='NHNM/NLNM')
plt.xlabel('Period (s)')
plt.xlim((2.,1000.))
plt.ylabel('Power (dB rel. 1 $(m/s^2)/Hz$)', fontsize=16)
plt.colorbar(cmap, ticks = vals, orientation='horizontal', pad=0.2, label='Sensitivity ($m/T{s^2}$)')
plt.text(0.8, -42, '(b)')
#plt.ylim((-50, 75))
plt.savefig('TEST.png', format='PNG')
plt.show()

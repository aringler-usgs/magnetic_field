#!/usr/bin/env python
from obspy.core import read, UTCDateTime, Stream
from obspy.signal.spectral_estimation import PPSD
import matplotlib.pyplot as plt


stas = ['ANMO', 'CASY', 'COLA', 'QSPA', 'SBA', 'SFJD']
chans = ['LF1', 'LF2', 'LFZ']

def comp_ppsd(sta, debug=True):
    net = "IU"
    for chan in chans:
        print('On:' + sta + ' ' + chan)
        stime = UTCDateTime('2016-280T00:00:00')
        etime = UTCDateTime('2019-280T00:00:00')
        ctime = stime
        if sta == 'QSPA':
            sen = 3.43*10**9
        else:
            sen = 41.943
        paz = {'gain': 1., 'poles':[], 'zeros':[], 'sensitivity': sen}
        while ctime < etime:
            if debug:
                print(ctime)
            try:
            #if True:
                st = Stream()
                nctime = ctime
                for addday in range(5):
                    nctime = ctime + addday*24*60*60
                    st += read('/msd/' + net + '_' + sta +'/' + str(nctime.year) +
                 '/' + str(nctime.julday).zfill(3) + '/*' + chan + '*' )
                if debug:
                    print(st)
            except:
                ctime +=5*24*60*60
                continue
            st.merge(fill_value=0)
            if 'ppsd' not in vars():
                ppsd = PPSD(st[0].stats, paz, db_bins=(-100,100,1.), ppsd_length=2**14,
                        period_smoothing_width_octaves = 0.50, special_handling="ringlaser")
            ppsd.add(st)
            ctime +=5*24*60*60
        ppsd.save_npz('PDF_DATA_' + net +'_' + sta + '_' + chan + ".npz")
        ppsd.plot(filename='PPSD_' + net +'_' + sta + '_' + chan + '.PNG',
                  show_noise_models=False, period_lim=(2., 1000.))
        del ppsd
    return


#comp_ppsd("SFJD")

from multiprocessing import Pool
pool = Pool(6)
pool.map(comp_ppsd, stas)

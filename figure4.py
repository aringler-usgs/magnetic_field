#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
#fig= plt.figure(1,figsize=(5,10))

# code for figure 4.  This uses the computed PSDs from 3

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font',size=12)
fig, ax = plt.subplots(figsize=(8,10))
stas, vals1, vals2 = [], [], []

f = open('day_ALL_data','r')
	

for line in f:
	line = line.replace(' ','')
	stas.append(line.split(',')[0])
	val1 = float(line.split(',')[1])
	if val1 > 3*10**-8:
		val1 = 3*10**-8
	val2 = float(line.split(',')[2])
	if val2 > 3*10**-8:
		val2 = 3*10**-8
	vals1.append(val1)
	vals2.append(val2)
f.close()
ax.semilogx(vals1,range(len(stas)),'.')
ax.semilogx(vals2,range(len(stas)),'.')

plt.xlabel('RMS Acceleration ($m/s^2$)')
plt.ylabel('Station code')
ax.set_yticks(range(len(stas)))
ax.set_yticklabels(stas, fontsize=6)
for tick in ax.yaxis.get_major_ticks()[1::2]:
    tick.set_pad(21)
#plt.yticks([])
idx = 0.
mv = min(vals1)
Mv = max(vals2)
step = (Mv - mv)/10. 
print(step)
loc = mv
for trip in zip(stas,vals1,vals2):
    loc += step
    loc = loc % Mv
    plt.semilogx([trip[1],trip[2]],[idx, idx], alpha=.3, color='k')
    #plt.text(np.log10(loc), idx, trip[0], fontsize=8)
    idx += 1
plt.xlim((mv,Mv))
plt.ylim((0,len(stas)))
plt.savefig('New_plt_rms.png', format='PNG', dpi=400)
plt.show()

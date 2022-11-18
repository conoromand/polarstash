import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.cm as cm
import h5py

plt.rcdefaults()
plt.rc('lines', linewidth=1.5)
plt.rc('axes', labelsize=25)
plt.rc('xtick', labelsize=23)
plt.rc('ytick', labelsize=23)
plt.rc('text', usetex=True)
plt.rc("legend", fontsize=18)

filename = "outputs/nph1.0e+07_L46large_t5e6"
		
DAY = 86400.

ls = ['-','-','-','-','-']
lw = [2.5,2,1.5,1,0.5]
label = ['1 core','1 core restarted','4 cores','4 cores restarted']

# Number of observers

fig2=plt.figure(figsize=(12,7))
ax2 = fig2.add_subplot(1,1,1)
ax2.set_xlabel('Time since explosion (days)')
ax2.set_ylabel('Lbol (erg / s)')
ax2.set_yscale('log')

#### COLORS

a = h5py.File(filename+'.hdf5', 'r')

lbol = a['observables']['lbol']
stokes = a['observables']['stokes']
wave = np.array(a['observables']['wave'])
time = np.array(a['observables']['time'])/DAY

Nobs = stokes.shape[0]

cmap = plt.get_cmap('jet_r',Nobs)
norm = mpl.colors.Normalize(vmin=0,vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

for j in range(0,Nobs):
	ax2.plot(time,lbol[j],color=cmap(j),lw=2,zorder=1,label='$\\cos\\theta_{\\rm obs}=%.1f$'%(1/(Nobs-1)*j))

# SUZUKI & MAEDA

if 'L46large' in filename:
	a = genfromtxt('lc_L46large_S21.txt')
	t = a[:,0] / DAY
	lbol_th = [[0 for i in range(len(t))] for j in range(11)]
	lbol_nonth = [[0 for i in range(len(t))] for j in range(11)]

	for i in range(11):
		lbol_th[i] = a[:,1+i]
		lbol_nonth[i] = a[:,12+i]
		
		ax2.plot(t,lbol_th[i],color='grey',lw=0.5,zorder=0)

	# DEPOSITION
	tdep = np.arange(0,1200,1)
	ldep = 1e46 / (1+tdep/(1e6/DAY))**2

	ax2.set_xlim(20,350)
	ax2.set_ylim(1e42,1e46)

else:
	tdep = np.arange(0,1200,1)
	ldep = 1e48 / (1+tdep/(1e4/DAY))**2

	ax2.set_xlim(0,60)
	ax2.set_ylim(2e41,2e45)

ax2.plot(tdep,ldep,ls='--')

subplots_adjust(hspace=0.05)
subplots_adjust(wspace=0.05)

ax2.legend(loc=1,frameon=False)

plt.savefig('plots/' + filename[11:] + '_lcBol.pdf')

fig2.show()


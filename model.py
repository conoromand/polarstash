import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pylab import *
import h5py

CLIGHT=2.99792458e10
DAY=86400

plt.rcdefaults()
plt.rc('lines', linewidth=1.5)
plt.rc('axes', labelsize=28)
plt.rc('xtick', labelsize=23)
plt.rc('ytick', labelsize=23)
plt.rc('text', usetex=True)
plt.rc("legend", fontsize=12)

DAY = 86400.

t_ref = 1e5
t = 5 * DAY

filename = "outputs/nph1.0e+07_L48large"
			
a = h5py.File(filename+'.hdf5', 'r')

times = np.array(a['observables']['time']) / DAY
vel = np.array(a['model']['velocity']) / CLIGHT
rho = np.array(a['model']['density'])
T = np.array(a['model']['temperature'])

# Scale rho and select T at requested time
rho = rho * (t_ref / t)**3  # reference time for rho is 1 DAY
T = T * (t_ref / t)
kabs = 6.64e22 * 8**3 / 16**2 * rho / pow(T,3.5)
ksc = 0.2

# Plot

fig = plt.figure(figsize=(16,8))
ax1 = fig.add_subplot(121,aspect=1)
ax1.set_ylabel('$v_z/c$')
ax3 = fig.add_subplot(122,aspect=1)
ax3.set_yticklabels([])
[i.set_xlabel('$v_y/c$') for i in [ax1,ax3]]

no_bkg = np.where(rho>1e-30)

s1=ax1.pcolormesh(vel,vel,log10(rho),cmap='jet',vmin=amin(log10(rho[no_bkg])),vmax=amax(log10(rho)))
s3=ax3.pcolormesh(vel,vel,log10(T),cmap='jet',vmin=amin(log10(T[no_bkg])),vmax=5)#amax(log10(T)))
#s1=ax1.pcolormesh(vel,vel,log10(kabs),cmap='jet',vmin=-10,vmax=0)#,vmax=amax(log10(rho)))
#s3=ax3.pcolormesh(vel,vel,log10(kabs/ksc),cmap='RdBu_r',vmin=-3,vmax=3)#,vmax=amax(log10(T)))

cbar1=plt.colorbar(s1,cax=fig.add_axes([0.125, 0.9, 0.35, 0.025]),orientation='horizontal')
cbar1.set_label('Log$\\rho$ (g cm$^{-3})$',size=25,labelpad=-70)
cbar3=plt.colorbar(s3,cax=fig.add_axes([0.55, 0.9, 0.35, 0.025]),orientation='horizontal')
cbar3.set_label('Log$T$ (K)',size=25,labelpad=-70)

fig.show()

#fig.savefig("/Users/mbull/Desktop/lc.pdf",bbox_inches='tight')

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from pylab import *
import h5py
from scipy import signal
from rebin import rebin

plt.rcdefaults()
plt.rc('lines', linewidth=1.5)
plt.rc('axes', labelsize=25)
plt.rc('xtick', labelsize=23)
plt.rc('ytick', labelsize=23)
plt.rc('text', usetex=True)
plt.rc("legend", fontsize=10)

DAY = 86400.
CLIGHT = 2.99792458e10

filename = "outputs/nph1.0e+07_L48large"

t1 = 30
t2 = t1*1.2

w_map = 10000
iobs_map = 2
window_width=10

####### SPECTRA

a = h5py.File(filename+'.hdf5', 'r')

stokes = a['observables']['stokes']

Nobs = stokes.shape[0]
Ntime = stokes.shape[1]
Nwave = stokes.shape[2]

cmap = plt.get_cmap('jet_r',Nobs)
norm = mpl.colors.Normalize(vmin=0,vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

nwave = np.array(a['observables']['wave'])
time = np.array(a['observables']['time'])/DAY

#cumsum_wave = np.cumsum(np.insert(nwave, 0, 0)) 
#wave = (cumsum_wave[window_width:] - cumsum_wave[:-window_width]) / window_width
wave = rebin(nwave, factor=(window_width), func=np.mean)

Iall = [0 for i in range(Nobs)] 
qperc = [0 for i in range(Nobs)] 
uperc = [0 for i in range(Nobs)] 
qave = [0 for i in range(Nobs)] 
uave = [0 for i in range(Nobs)] 

ind_t1 = np.argmin(fabs(t1-time))
ind_t2 = np.argmin(fabs(t2-time))

for j in range(0,Nobs):

	I = np.nan_to_num(np.mean(stokes[j,ind_t1:ind_t2,:,0],axis=0))+1E-30
	Q = np.nan_to_num(np.mean(stokes[j,ind_t1:ind_t2,:,1],axis=0))
	U = np.nan_to_num(np.mean(stokes[j,ind_t1:ind_t2,:,2],axis=0))

	#cumsumI=np.cumsum(np.insert(I, 0, 0))
	#cumsumq=np.cumsum(np.insert(Q / I * 100., 0, 0))
	#cumsumu=np.cumsum(np.insert(U / I * 100., 0, 0))
	
	#Iall[j] = (cumsumI[window_width:] - cumsumI[:-window_width]) / window_width
	#qperc[j] = (cumsumq[window_width:] - cumsumq[:-window_width]) / window_width
	#uperc[j] = (cumsumu[window_width:] - cumsumu[:-window_width]) / window_width
	
	Iall[j]=rebin(np.insert(I, 0, 0), factor=(window_width), func=np.mean)
	qperc[j]=rebin(np.insert(Q / I * 100., 0, 0), factor=(window_width), func=np.mean)
	uperc[j]=rebin(np.insert(U / I * 100., 0, 0), factor=(window_width), func=np.mean)

	qave[j] = sum(Q)/sum(I)*100.
	uave[j] = sum(U)/sum(I)*100.

	print(qave[j],uave[j])


####### MAPS

st = a['observables']['stokes_map']

Nobs_map = st.shape[0]
Ntime_map = st.shape[1]
Nwave_map = st.shape[2]
nx_map = st.shape[3]

v_map = a['observables']['vel_map']
time_map = a['observables']['time_map']
wave_map = a['observables']['wave_map']

vmax = max(fabs(v_map[:]))

ind_t = int( (0.5*(t1+t2)-time_map[0]) / (time_map[1]-time_map[0]) )

ind_wave_map = int( (w_map-wave_map[0]) / (wave_map[1]-wave_map[0]) )

Im = st[iobs_map,ind_t,ind_wave_map,:,:,0]
Qm = st[iobs_map,ind_t,ind_wave_map,:,:,1]
Um = st[iobs_map,ind_t,ind_wave_map,:,:,2]

Imap = np.reshape(Im, (-1, int(nx_map)))
Qmap = np.reshape(Qm, (-1, int(nx_map)))
Umap = np.reshape(Um, (-1, int(nx_map)))

qmap = Qmap#/Imap*100.
umap = Umap#/Imap*100.

Imap[np.isinf(Imap)] = 0
qmap[np.isnan(qmap)] = 0
umap[np.isnan(umap)] = 0


#### PLOT

pol_max = 10
lmin = 500
lmax = 10000

fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(321)
ax1.set_xscale('log')
ax1.set_xlim(lmin,lmax)
ax1.set_ylim(0,1.1)
ax1.set_ylabel('Flux (normalized)')

ax2 = fig.add_subplot(323)
ax2.set_xscale('log')
ax2.set_xlim(lmin,lmax)
ax2.set_ylim(-pol_max,pol_max)
ax2.set_ylabel('Q ($\%$)')

ax3 = fig.add_subplot(325)
ax3.set_xscale('log')
ax3.set_xlim(lmin,lmax)
ax3.set_ylim(-pol_max,pol_max)
ax3.set_xlabel('Wavelength (\AA)')
ax3.set_ylabel('U ($\%$)')

ax2.axhline(y=0,ls='--',color='black')
ax3.axhline(y=0,ls='--',color='black')

color = ['black',"purple","green",'blue','red','cyan']

lab = ['Equator','Pole']

for j in range(0,Nobs):
	
	if j==iobs_map:
		lw = 3
		zord = 1#1
	else:
		lw = 1
		zord = 0

	ax1.plot(wave,Iall[j]/max(np.ravel(Iall)),color=color[j],ds="steps",lw=lw
		,zorder=zord,label='$\\cos\\theta_{\\rm obs}=%.1f$'%(1/(Nobs-1)*j))
	
	ax2.plot(wave,qperc[j],color=color[j],ds="steps",lw=lw,zorder=zord)
	
	ax3.plot(wave,uperc[j],color=color[j],ds="steps",lw=lw,zorder=zord)
	#ax3.axhline(y=uave[j],ls='-',lw=2,color=color[j])y = signal.savgol_filter(Iall, 53, 3)



[ax.axvspan(xmin=wave_map[ind_wave_map],xmax=wave_map[ind_wave_map+1],color='grey',alpha=0.3) for ax in [ax1,ax2,ax3]]


ax1.legend(ncol=1,loc=1,frameon=False)

#### MAP


ax1_map = fig.add_subplot(322,aspect=1)
ax2_map = fig.add_subplot(324,aspect=1)
ax3_map = fig.add_subplot(326,aspect=1)
[ax.set_xlim(v_map[0],v_map[-1]) for ax in [ax1_map,ax2_map,ax3_map]]
[ax.set_ylim(v_map[0],v_map[-1]) for ax in [ax1_map,ax2_map,ax3_map]]

maxi = max(np.ravel(Imap))
s1=ax1_map.pcolormesh(v_map,v_map,Imap,cmap='hot',rasterized=True)#,vmin=42,vmax=44) #min(logImap[ind_no_bkg])
cbar1=plt.colorbar(s1,cax=fig.add_axes([0.9,0.66,0.025,0.21]))

maxq = max(max(np.ravel(qmap)),fabs(min(np.ravel(qmap))))
s2=ax2_map.pcolormesh(v_map,v_map,qmap,vmin=-maxq,vmax=maxq,cmap='RdBu_r',rasterized=True)
cbar2=plt.colorbar(s2,cax=fig.add_axes([0.9,0.395,0.025,0.21]))

maxu = max(max(np.ravel(umap)),fabs(min(np.ravel(umap))))
s3=ax3_map.pcolormesh(v_map,v_map,umap,vmin=-maxu,vmax=maxu,cmap='RdBu_r',rasterized=True)
cbar3=plt.colorbar(s3,cax=fig.add_axes([0.9,0.13,0.025,0.21]))



subplots_adjust(hspace=0.25)
subplots_adjust(wspace=0.05)

#ax.legend(loc=4,frameon=False)
plt.savefig('plots/' + filename[11:] + '_polar_D' + str(t1) + '.pdf')

fig.show()

#fig.savefig('plots/pol.pdf',bbox_inches='tight')
#fig.savefig("/Users/mbull/Desktop/pol_L46_150d.pdf",bbox_inches='tight')

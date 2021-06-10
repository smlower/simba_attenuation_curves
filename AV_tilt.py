#bump_scattered

#to plot bump strength against the scattered light component of the 2175 flux

import matplotlib
matplotlib.use('Agg')

from pylab import setp 

import numpy as np
import ipdb,pdb
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput

from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants

from attenuation_fit_funcs import *
from scipy.optimize import curve_fit
from astropy.convolution import convolve, Box1DKernel

from glob2 import glob
import itertools

from get_quantities import get_all,get_tauV

from attenuation_fitter_function import *
from get_delta_bump import get_delta_bump


def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx


#===========================================
#MODIFIABLE HEADER

master_file = 'master_attenuation_curves.npz'


COSMOLOGICAL = False
#attenuation_directories = ['/ufrc/narayanan/desika.narayanan/pd_runs/mufasa/m25n512/fh_qr/quick_look_attenuation/snap085/npzfiles']
#for pd stuff
z = 0.0001
minwavbeta = 1000
maxwavbeta = 2200

minwavlir = 8
maxwavlir = 1000

uv_lam = 1670
uv_lam = 1600
b_lam = 4420
v_lam = 5500


#master_file = 'master_attenuation_curves.m25n256_085.npz'

#===========================================
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""



data = np.load(master_file)
master_fit = data['master_fit']
master_delta_calzetti = data['master_delta_calzetti']
master_ebump = data['master_E_bump_kc13']
master_mstar = data['master_mstar']
x = data['x']
wav_af = data['wav_af']*u.micron


'''
#go through and remove all the -1 elements
good_idx = []
bad_idx = []
for i in range(len(master_fit)):
    if master_fit[i][0] == -1:
        bad_idx.append(i)
    else:
        good_idx.append(i)

master_fit = master_fit[good_idx]
master_delta_calzetti = master_delta_calzetti[good_idx]
master_ebump = master_ebump[good_idx]
master_mstar = master_mstar[good_idx]

w09 = np.where( (np.log10(master_mstar) >= 8.75) & (np.log10(master_mstar) <= 9.25))[0].tolist()
w10 = np.where( (np.log10(master_mstar) >= 9.75) & (np.log10(master_mstar) <= 10.25))[0].tolist()
w11 = np.where( (np.log10(master_mstar) >= 10.75) & (np.log10(master_mstar) <= 12.25))[0].tolist()

average_curve_09 = np.zeros(len(x))
for i in range(len(x)):
    average_curve[09] = np.median(master_fit[w09]


average_curve_09 = np.median(master_fit[w09],axis=0)
pdb.set_trace()    
'''


vband_idx = find_nearest(wav_af.to(u.angstrom).value,v_lam)
uvband_idx = find_nearest(wav_af.to(u.angstrom).value,uv_lam)
b_band_idx = find_nearest(wav_af.to(u.angstrom).value,b_lam)

master_tau_v = []
for i in range(master_fit.shape[0]):
    master_tau_v.append(master_fit[i][vband_idx])


#plot AV vs slope deviatino
data = np.loadtxt('datafiles/salim/av_slope.dat',usecols=(0,1,2,3))
salim_av = data[:,0]
salim_delta_average = data[:,1]
salim_delta_average_plus = data[:,2]
salim_delta_average_minus = data[:,3]

#get KC13 stuff

from scipy.io import readsav
kc13_data = readsav('dust.save')
kc13_ebump = kc13_data['b_e_b']
kc13_b_del = kc13_data['b_del']


#get tress data
tressfile = 'datafiles/SHARDS_SFDust.txt'
tress_data = np.loadtxt(tressfile,skiprows=3)
tress_ebump = tress_data[:,5]
tress_b_del = tress_data[:,4]
tress_AV = tress_data[:,6]

fig = plt.figure()
fig.subplots_adjust(hspace=.65)
#fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.scatter(tress_AV,tress_b_del,color='blue',edgecolor='crimson',s=6,alpha=0.55,label='Tress e\
t al. 2018')
ax1.scatter(1.086*np.asarray(master_tau_v),master_delta_calzetti,color='xkcd:sea blue',edgecolor='black',s=15,label='This Study')

ax1.fill_between(salim_av,salim_delta_average_plus,salim_delta_average_minus,alpha=0.5,label='Salim et al. 2018')
ax1.set_xlabel(r'A$_\mathrm{V}$')
ax1.set_ylim([-1.5,1])
ax1.set_xlim([0,2])
ax1.set_ylabel(r'$\delta_\mathrm{cal}$')

salmon_EBV = np.linspace(0,0.6,100)
salmon_delta = 0.62*np.log10(salmon_EBV)+0.26
salmon_delta_low = (0.62-0.05)*np.log10(salmon_EBV)+(0.26+0.05)
salmon_delta_high = (0.62+0.05)*np.log10(salmon_EBV)+(0.26-0.05)
ax1.fill_between(salmon_EBV*3.1,salmon_delta_high,salmon_delta_low,color='coral',alpha=0.5,label='Salmon et al. 2016')

plt.legend(fontsize=6,loc=4)
plt.grid()


#KC13 COMPARISON
ax2 = fig.add_subplot(222)

#tress data
ax2.scatter(tress_b_del,tress_ebump,color='blue',edgecolor='crimson',s=6,alpha=0.55,label='Tress et al. 2018')

w = np.where(np.asarray(master_tau_v)*1.086 >= 0)[0]
ax2.scatter(master_delta_calzetti[w],master_ebump[w],color='xkcd:sea blue',edgecolor='black',label='This Study',s=15)
ax2.set_ylim([0,3])
ax2.set_xlim([-1.5,1])
#linfit this data and plot it
from numpy.linalg import lstsq

A = np.vstack([master_delta_calzetti[w],np.ones(len(master_delta_calzetti[w]))]).T
m,c = np.linalg.lstsq(A,master_ebump[w])[0]
x = np.linspace(-2,2,100)
y = m*x+c
ax2.plot(x,y,color='xkcd:sea blue')
print 'm = ',m
print 'c = ',c
#ax1.scatter(kc13_b_del[1,:],kc13_ebump[1,:],color='red')

ax2.errorbar(kc13_b_del[1,:],kc13_ebump[1,:],yerr=kc13_ebump[1,:]-kc13_ebump[0,:],color='orange',fmt='o')
#plt.errorbar(kc13_b_del[1,:],kc13_ebump[1,:],color='red')
kc13_delta = np.linspace(-2,2,100)
kc_13_ebump = 0.85-(1.9*kc13_delta)
ax2.plot(kc13_delta,kc_13_ebump,color='orange',label='Kriek & Conroy 2013')

ax2.set_xlabel(r'$\delta_\mathrm{cal}$')
ax2.set_ylabel(r'E$_\mathrm{bump}$')
plt.grid()
plt.legend(fontsize=6,loc=2)

'''
#SALIM MASSES COMPARISON
ax3 = fig.add_subplot(223)
data09= np.loadtxt('datafiles/salim/curve_sf_logm09.dat',usecols=(0,1,2))
salim_lam09 = data09[:,0]*u.angstrom
salim_av09 = data09[:,1]
salim_avstd09 = data09[:,2]
salim_x = 1./(salim_lam09.to(u.micron)).value
ax3.plot(salim_x,salim_av09)
'''


fig.savefig('salim1.png',dpi=300)

import numpy as np
from hyperion.model import ModelOutput
from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants
import h5py
from hyperion_input_seds import get_input_seds
from glob import glob
import sys, os

galaxy_idx = int(sys.argv[1])
snap_str = '87'

galaxy_num = galaxy_idx #np.load('/orange/narayanan/s.lower/prospector/simba_runs/simba_galaxy_SFRcut.npz')['ID'][int(galaxy_idx)]

sed_directory = '/blue/narayanan/desika.narayanan/pd_runs/simba/m25n512/NSF_AAG_2020/snap087/'
output_directory = '/orange/narayanan/s.lower/simba/galaxy_properties/attenuation_curves/snap087/'
#========================================================



def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

def sed_downsample(source_lam,lum_source,wav_rest):
    source_lam_downsampled = np.zeros(len(wav_rest))
    lum_source_downsampled = np.zeros(len(wav_rest))
    for i in range(len(wav_rest)):
        idx = find_nearest(source_lam,wav_rest[i].value)
        source_lam_downsampled[i] = source_lam[idx]
        lum_source_downsampled[i] = lum_source[idx]

    return source_lam_downsampled,lum_source_downsampled

#sed_file = glob(sed_directory+'/*galaxy'+str(galaxy_num)+'.rtout.sed')[0]
#print(sed_file)
#stellar_file = glob(sources_directory+'/*.galaxy'+str(galaxy_num)+'.rtout.sed')[0]

print('snap305.galaxy'+str(galaxy_num)+'.rtout.sed')
comp_sed = ModelOutput(sed_directory+'/snap.galaxy'+str(galaxy_num)+'.rtout.sed')
wav_rest_sed,lum_obs_sed = comp_sed.get_sed(inclination=0,aperture=-1)
wav_rest_sed =wav_rest_sed[::-1]* u.micron #wav is in micron                                                                                                 
lum_obs_sed = lum_obs_sed * u.erg/u.s


#stellar_data = np.load(sed_directory+'/stellar_seds_galaxy'+str(galaxy_num)+'.npz')
stellar_data = np.load('/orange/narayanan/s.lower/simba/pd_runs/snap305_mono/snap305/stellar_seds_galaxy'+str(galaxy_num)+'.npz')
stellar_lam = stellar_data['lam'][::-1]
stellar_lum = stellar_data['flam'][::-1]


dwn = sed_downsample(stellar_lam, stellar_lum, wav_rest_sed)

print('downsampled lam: ',dwn[0])
print('lum: ',dwn[1])
print('obs lum: ',lum_obs_sed.to(u.Lsun).value)

extinction = lum_obs_sed.to(u.Lsun).value / (dwn[0] * dwn[1])
tau = -1.0 * np.log(extinction)

vband_st = find_nearest(stellar_lam, 0.55)
vband_ob = find_nearest(wav_rest_sed.value, 0.55)


print('stellar: ',stellar_lam)
print('obs: ',wav_rest_sed)

extinction_v = lum_obs_sed[vband_ob].to(u.Lsun).value / (stellar_lum[vband_st] * stellar_lam[vband_st])
tau_v = -1.0 * np.log(extinction_v)


print('tau_v: ',tau_v)




'''

comp_sed = ModelOutput(sed_file)
wav_rest_sed,dum_lum_obs_sed = comp_sed.get_sed(inclination=0,aperture=-1)
wav_rest_sed =wav_rest_sed* u.micron #wav is in micron 
nu_rest_sed = constants.c.cgs/wav_rest_sed.cgs
lum_obs_sed = dum_lum_obs_sed
lum_obs_sed = lum_obs_sed * u.erg/u.s
nu_rest_sed = constants.c.cgs/(wav_rest_sed.to(u.cm))
fnu_obs_sed = lum_obs_sed.to(u.Lsun)
fnu_obs_sed /= nu_rest_sed.to(u.Hz)
fnu_obs_sed = fnu_obs_sed.to(u.Lsun/u.Hz)

print('size of observed sed: ',np.shape(wav_rest_sed))

#stellar_sed = np.load(stellar_file)
#nu_rest_stellar = stellar_sed['nu'] #Hz
#fnu_rest_stellar = stellar_sed['fnu'] #Lsun/Hz
#fnu_rest_stellar = fnu_rest_stellar * u.Lsun/u.Hz
#nu_rest_stellar = nu_rest_stellar * u.Hz
#lum_rest_stellar = fnu_rest_stellar * nu_rest_stellar
#lam_rest_stellar = constants.c.cgs/(nu_rest_stellar)

#source_lam_downsampled,lum_source_downsampled = sed_downsample(lam_rest_stellar.to(u#.micron),lum_rest_stellar.to(u.Lsun),wav_rest_sed.to(u.micron))
#source_lam_downsampled = source_lam_downsampled * u.micron
#lum_source_downsampled  = lum_source_downsampled * u.Lsun



source_nu,source_fnu = get_input_seds(stellar_file)
source_nu = source_nu * u.Hz
source_lam = constants.c.cgs/source_nu
source_lam = source_lam.to(u.micron)
source_fnu = source_fnu * u.Lsun/u.Hz
lum_source = source_fnu * source_nu


source_lam = source_lam[find_nearest(source_lam, 1.0):find_nearest(source_lam,0.1)+1]
source_lum = lum_source[find_nearest(source_lam, 1.0):find_nearest(source_lam,0.1)+1]
#print('shape of stellar sed: ',np.shape(source_lam))


#source_lam_downsampled,lum_source_downsampled = sed_downsample(source_lam.to(u.micron),lum_source.to(u.Lsun),wav_rest_sed.to(u.micron))
#source_lam_downsampled = source_lam_downsampled * u.micron                                                                 
#lum_source_downsampled  = lum_source_downsampled * u.Lsun



#stellar_sed = ModelOutput(stellar_file)
#wav_rest_stellar,dum_lum_stellar_sed = stellar_sed.get_sed(inclination='all',aperture=-1)
#wav_rest_stellar =wav_rest_stellar* u.micron #wav is in micron                                                     
#nu_rest_stellar = constants.c.cgs/wav_rest_stellar.cgs
#lum_obs_stellar = dum_lum_stellar_sed
#lum_obs_stellar = lum_obs_stellar * u.erg/u.s
#nu_rest_stellar = constants.c.cgs/(wav_rest_stellar.to(u.cm))
#fnu_obs_stellar= lum_obs_stellar.to(u.Lsun)
#fnu_obs_stellar/= nu_rest_stellar.to(u.Hz)
#fnu_obs_stellar= fnu_obs_stellar.to(u.Lsun/u.Hz)


#idx = np.where((lam_rest_stellar.to(u.micron).value > 0.1) & (lam_rest_stellar.to(u.micron).value < 1.0))

#stellar_wav = lam_rest_stellar[idx]
#stellar_flux = lum_obs_stellar[idx]



#find location of vband in lam array
vband_idx = find_nearest(wav_rest_sed.to(u.angstrom).value,3000)
#vband_extinction = np.zeros(fnu_obs_sed.shape[0])
vband_extinction = lum_obs_sed[vband_idx].cgs/lum_source[vband_idx].cgs




#calculate extinction values
#extinction = np.zeros([lum_obs_sed.shape[0],len(wav_rest_sed)])
extinction = lum_obs_sed.cgs/source_lum.cgs'''


#tau is calculated by: e^-tau = I/I_0
#tau = -1.*np.log(extinction)
#tau_v = -1.*np.log(vband_extinction)

inverse_wavelength = 1./wav_rest_sed.to(u.micron)


outfile = output_directory+'attenuation_galaxy'+str(galaxy_num)+'_corrected.npz'

np.savez(outfile,wav_rest=wav_rest_sed[::-1].value,inverse_wavelength = inverse_wavelength.value,tau=tau,tau_v=tau_v)







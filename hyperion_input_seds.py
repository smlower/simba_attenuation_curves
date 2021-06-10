from hyperion.model import ModelOutput,Model
import numpy as np
#import pdb,ipdb
import astropy.units as u
from astropy import constants
#import h5py


def get_input_seds(file, gal_id=None):

    totallum = 0
    m = Model()
    m.use_sources(file, gal_id)
    nsources = len(m.sources)
    #test = h5py.File(file, 'r')
    #nsources = len(test['Input']['Sources'])
    
    for i in range(nsources):
        #source_num = 'source_'+"{:05d}".format(i)
        tempnu = m.sources[i].spectrum["nu"]
        tempfnu = m.sources[i].spectrum["fnu"]
        
        if i == 0: fnu = np.zeros(len(tempnu))
        

        #now we need to scale this because the spectrum is just in
        #terms of an SSP, and we need to scale by the total luminosity
        #that wen t into the model (i.e. by the actual stellar mass
        #used in powderday).



       

        ssp_lum = np.absolute(np.trapz(tempnu,tempfnu))*constants.L_sun.cgs
        lum_scale = np.sum(m.sources[i].luminosity)/ssp_lum #we have to do np.sum in case the sources were in a collection
        tempfnu *= lum_scale.value
        



        for i in range(len(fnu)):
            fnu[i] += tempfnu[i]

    #ipdb.set_trace()
            
    
    return tempnu,fnu

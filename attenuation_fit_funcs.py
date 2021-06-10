import numpy as np
import astropy.units as u
from literature_attenuation import *

def peak_smooth(y):
    for i in range(len(y)-2):
        if (y[i+1] > y[i]) and (y[i+1] > y[i+2]):

            y[i+1] = np.mean([y[i],y[i+2]])
            return y


def noll_fitter(wav,E_bump,gamma,cr,delta,calzetti_norm):

    tau_calzetti = calzetti(wav)

    drude_num = E_bump * (wav**2.)*(gamma**2.)
    drude_denom = (wav**2.-2175.**2.)**2. + (wav**2.*gamma**2.)
    drude = drude_num/drude_denom

    f = calzetti_norm*((tau_calzetti) + drude * ((1.-1.12*cr)/4.05+1.))
    f = f*(wav/5500.)**delta


    return f

    
def noll_fitter_kc13(wav,E_bump,cr,delta,calzetti_norm):
    #same as noll_fitter, but we assume a gamma of 350 Angstrom as in Kriek & Conroy 2013
    tau_calzetti = calzetti(wav)
    
    gamma = 350. #angstrom
    drude_num = E_bump * (wav**2.)*(gamma**2.)
    drude_denom = (wav**2.-2175.**2.)**2. + (wav**2.*gamma**2.)
    drude = drude_num/drude_denom
    
    f = calzetti_norm*((tau_calzetti) + drude * ((1.-1.12*cr)/4.05+1.))
    f = f*(wav/5500.)**delta


    return f
        

def ir_fitter(x,a1,b1,gamma_ir):
    """

    :x is the inverse wavelength (angstrom) defined: 1e4 /
    wav.to(u.angstrom).value

    :a1,b1,gamma_ir take the form 

    f_ir = (a1*x**gamma_ir + b1*x**gamma_ir)/R_v

    """

    a = a1 * x**gamma_ir
    b = b1 * x**gamma_ir
    R_v = 3.1
    return (a+b)/R_v


def opt_fitter(x,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6,b7):
    y = x-1.82

    a= (1 + a1 * y - a2 * y**2 - a3 * y**3 +
        a4 * y**4 + a5 * y**5 - a6 * y**6 +
        a7* y**7)
    
    b = (b1 * y + b2 * y**2 + b3 * y**3 -
          b4 * y**4 - b5 * y**5 + b6 * y**6 -
          b7 * y**7)
    

    R_v = 3.1
    return (a+b)/R_v




def NUV_fitter(x,f_bump):

    R_v = 3.1

    '''
    tmp = (-0.0370 + 0.0469 * f_bump - 0.601 * f_bump / R_v + 0.542 / R_v)
    fa = (3.3 / x)**6. * tmp

    tmp = a3 * f_bump / ((x - a4)**2 + a5)
    a = a1 - a2 * x - tmp + fa
    tmp = b3 * f_bump / ((x - b4)**2 + b5)
    b = b1 + b2 * x + tmp
    '''


    tmp = (-0.0370 + 0.0469 * f_bump - 0.601 * f_bump / R_v + 0.542 / R_v)
    fa = (3.3 / x)**6. * tmp
    tmp = 0.104 * f_bump / ((x - 4.67)**2 + 0.341)
    a = 1.752 - 0.316 * x - tmp + fa
    tmp = 1.206 * f_bump / ((x - 4.62)**2 + 0.263)
    b = -3.09 + 1.825 * x + tmp

    return (a+b)/R_v


def FUV_fitter(x,a1,a2,b1,b2,f_bump,norm):

    fa = -1.*a1 * (x - 5.9)**2.0 - a2 * (x - 5.9)**3
    fb = b1 * (x - 5.9)**2. + b2 * (x - 5.9)**3
    tmp = 0.104 * f_bump / ((x - 4.67)**2 + 0.341)
    a = 1.752 - 0.316 * x - tmp + fa
    tmp = 1.206 * f_bump / ((x - 4.62)**2 + 0.263)
    b = -3.09 + 1.825 * x + tmp + fb

    alam = (a+b)/3.1
    return alam*norm

import numpy as np
import astropy.units as u

from attenuation_fit_funcs import *
from scipy.optimize import curve_fit



# fit parameters that for some reason were in a separate file????
Rv = 3.1
ir1 = 0.3
ir2 = 2.5
opt1 = 2.5
opt2 = 3.2
nuv1 = 3.2
nuv2 = 6.0
fuv1 = 5.9
fuv2 = 10





def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx



def fit_func(npzfile,single_powerlaw=False):

#====================================
#manually load the attenuation curve
#====================================


    data = np.load(npzfile)
    tau = data['tau']
    wav = data['wav_rest']*u.micron
    

#====================================
#load up the reference relations
#====================================
    
    tau_calzetti = calzetti(wav.to(u.angstrom).value)
    tau_conroy = conroy(wav.to(u.angstrom).value)
    
    x = 1e4 / wav.to(u.angstrom).value

 #====================================
 #linearly interpolate over Lya
 #====================================

    lyalpha_max_x = 7.8
    lyalpha_min_x = 8.7
    
    
    min_idx = find_nearest(x,lyalpha_min_x)
    max_idx = find_nearest(x,lyalpha_max_x)
    
    for idx in range(min_idx,max_idx):
        tau[idx] = np.interp(x[idx],[x[max_idx],x[min_idx]],[tau[max_idx],tau[min_idx]])
        
    lyalpha_max_x = 9
    lyalpha_min_x = 10
    
    min_idx = find_nearest(x,lyalpha_min_x)
    max_idx = find_nearest(x,lyalpha_max_x)

    for idx in range(min_idx,max_idx):
        tau[idx] = np.interp(x[idx],[x[max_idx],x[min_idx]],[tau[max_idx],tau[min_idx]])






 #====================================
 #fit the IR using conroy
 #====================================
    ir = (x >= ir1) & (x < ir2)
    ir_popt,ir_pcov = curve_fit(ir_fitter,x[ir],tau[ir])
    a_ir = ir_popt[0]*x[ir]**ir_popt[2]
    b_ir = ir_popt[1]*x[ir]**ir_popt[2]
    ir_fit = (a_ir+b_ir)/Rv

 #====================================
 #fit the optical using conroy
 #====================================

    opt = (x >= opt1) & (x < opt2)
    opt_popt,opt_pcov = curve_fit(opt_fitter,x[opt],tau[opt])
    a1_opt = opt_popt[0]
    a2_opt = opt_popt[1]
    a3_opt = opt_popt[2]
    a4_opt = opt_popt[3]
    a5_opt = opt_popt[4]
    a6_opt = opt_popt[5]
    a7_opt = opt_popt[6]
    b1_opt = opt_popt[7]
    b2_opt = opt_popt[8]
    b3_opt = opt_popt[9]
    b4_opt = opt_popt[10]
    b5_opt = opt_popt[11]
    b6_opt = opt_popt[12]
    b7_opt = opt_popt[13]
    y_opt = x[opt]-1.82
    a_opt= (1 + a1_opt * y_opt - a2_opt * y_opt**2 - a3_opt * y_opt**3 +
            a4_opt * y_opt**4 + a5_opt * y_opt**5 - a6_opt * y_opt**6 +
            a7_opt* y_opt**7)
    
    b_opt = (b1_opt * y_opt + b2_opt * y_opt**2 + b3_opt * y_opt**3 -
             b4_opt * y_opt**4 - b5_opt * y_opt**5 + b6_opt * y_opt**6 -
             b7_opt * y_opt**7)
    opt_fit = (a_opt+b_opt)/Rv
    
    #renormalization
    #opt_fit *= ir_fit[0]/opt_fit[-1]

    
 #====================================
 #fit the NUV and bump using Noll
 #====================================

 #try


    if single_powerlaw == False:
        nuv = (x >= nuv1) & (x < nuv2)
    else:
        nuv = (x >= np.min(x)) & (x <= np.max(x))



    popt,pcov = curve_fit(noll_fitter,wav[nuv].to(u.angstrom).value,tau[nuv])

    
    E_bump = popt[0]
    gamma = popt[1]
    cr = popt[2]
    delta = popt[3]
    calzetti_norm = popt[4]
    
    
    popt_kc13,pcov_kc13 = curve_fit(noll_fitter_kc13,wav[nuv].to(u.angstrom).value,tau[nuv])
    E_bump_kc13 = popt_kc13[0]
    delta_kc13 = popt_kc13[2]


    wav_angstrom = wav.to(u.angstrom).value
    
    drude_num = E_bump * (wav_angstrom[nuv]**2.)*(gamma**2.)
    drude_denom = (wav_angstrom[nuv]**2.-2175.**2.)**2. + (wav_angstrom[nuv]**2.*gamma**2.)
    drude = drude_num/drude_denom

    wav_nuv = wav[nuv].to(u.angstrom)
    idx1 = find_nearest(wav_nuv.to(u.angstrom).value,2175-gamma)
    idx2 = find_nearest(wav_nuv.to(u.angstrom).value,2175+gamma)

    drude_int = np.trapz(drude,wav_angstrom[nuv])
    #drude_int = np.trapz(drude[idx1:idx2],wav_nuv[idx1:idx2].value)


    A_fit = calzetti_norm*((tau_calzetti[nuv]) +drude * ((1.-1.12*cr)/Rv+1.))
    A_fit = A_fit*(wav_angstrom[nuv]/5500.)**delta
    
 #renormalization
    
    #A_fit *= opt_fit[0]/A_fit[-1]


     #====================================
 #fit the fuv using conroy
 #====================================
    fuv = (x>=fuv1) & (x<fuv2)
    fuv_popt,fuv_pcov = curve_fit(FUV_fitter,x[fuv],tau[fuv])
    a1_fuv = fuv_popt[0]
    a2_fuv = fuv_popt[1]
    b1_fuv = fuv_popt[2]
    b2_fuv = fuv_popt[3]
    f_bump = fuv_popt[4]
    norm = fuv_popt[5]

    fa = -1.*a1_fuv * (x[fuv] - fuv1)**2.0 - a2_fuv * (x[fuv] - fuv1)**3
    fb = b1_fuv * (x[fuv] - fuv1)**2. + b2_fuv * (x[fuv] - fuv1)**3
    tmp = 0.104 * f_bump / ((x[fuv] - 4.67)**2 + 0.341)
    a_fuv = 1.752 - 0.316 * x[fuv] - tmp + fa
    tmp = 1.206 * f_bump / ((x[fuv] - 4.62)**2 + 0.263)
    b_fuv = -3.09 + 1.825 * x[fuv] + tmp + fb

    alam = (a_fuv+b_fuv)/Rv
    fuv_fit = alam*norm

    fuv_norm = norm #quanitty to return

    #renornamalization
    #fuv_fit *= A_fit[0]/fuv_fit[-1]



 #except RuntimeError:
 #    print 'RunTimeError: excepting'
 #    A_fit = tau_calzetti
    
    actual_curve = ( (x),(ir),(ir_fit),(opt),(opt_fit),(nuv),(A_fit),(fuv),(fuv_fit))
    return E_bump,delta,gamma,fuv_norm,drude_int,actual_curve,E_bump_kc13,delta_kc13

    
    

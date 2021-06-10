import numpy as np
import tqdm
import glob
import pandas as pd

pd_curves = glob.glob('/orange/narayanan/s.lower/simba/pd_runs/attenuation_curves/snap305/attenuation_curve.305_galaxy*.npz')

pd_av = []
pd_lambda = []
for i in tqdm.tqdm(range(len(pd_curves))):
    data = np.load(pd_curves[i])
    taus = data['tau']
    pd_av.append(taus)
    pd_lambda.append(data['wav_rest'])


data = {'tau_list': pd_av, 'wavelengths': pd_lambda}
s = pd.DataFrame(data)
s.to_pickle('/orange/narayanan/s.lower/simba/pd_runs/attenuation_curves/snap305/combined_data.pkl')

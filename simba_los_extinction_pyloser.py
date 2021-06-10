import caesar, yt
from caesar.pyloser import pyloser
import numpy as np
import sys


#direction = sys.argv[1]

snap_dir = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/'
filt_dir = '/orange/narayanan/s.lower/simba/filtered_snapshots/snap087/'

ds = yt.load(snap_dir+'/snapshot_087.hdf5')
obj = caesar.load(snap_dir+'/Groups/caesar_0087_z5.024.hdf5')
#ds = yt.load(filt_dir+'/galaxy_5.hdf5')
print('initializing pyloser')

for direction in ['x', 'y', 'z']:

    pylose = pyloser.photometry(obj, obj.galaxies[:1000], ds=ds, ssp_table_file='SSP_Kroupa_EL.hdf5', view_dir=direction,nproc=16)


    print('running pyloser in '+str(direction)+' direction')
    test = pylose.run_pyloser()

    print('getting group Av list')
    av_per_group = []
    for i in range(1000):
        av_per_group.append(pylose.groups[i].group_Av)

    print('done. saving')
    np.savez('/orange/narayanan/s.lower/simba/los_AV_caesar_per_galaxy_'+direction+'dir_z5.npz', Av=av_per_group)



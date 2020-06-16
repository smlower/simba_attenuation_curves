import caesar, yt
from caesar.pyloser import pyloser
import numpy as np


snap_dir = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/'
filt_dir = '/orange/narayanan/s.lower/simba/filtered_snapshots/snap305/'

ds = yt.load(snap_dir+'/snapshot_305.hdf5')
obj = caesar.load(snap_dir+'/Groups/caesar_0305_z0.000.hdf5')
#ds = yt.load(filt_dir+'/galaxy_5.hdf5')
print('initializing pyloser')
pylose = pyloser.photometry(obj, obj.galaxies[:2000], ds=ds, ssp_table_file='SSP_Kroupa_EL.hdf5', view_dir='x',nproc=16)


print('running pyloser')
test = pylose.run_pyloser()

print('getting group Av list')
av_per_group = []
for i in range(2000):
    av_per_group.append(pylose.groups[i].group_Av)

print('done. saving')
np.savez('/orange/narayanan/s.lower/simba/los_AV_caesar_per_galaxy_xdir.npz', Av=av_per_group)



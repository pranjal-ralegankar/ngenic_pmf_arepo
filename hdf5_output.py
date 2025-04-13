from ics_to_hdf5 import ReadGadgetMultipleFiles, WriteGadget, ReadGadget
import numpy as np
import os
import sys
import re

if __name__ == "__main__":
    # The first command-line argument (sys.argv[0]) is the name of the script itself
    # The input data will be passed as additional arguments (sys.argv[1:])
    folder_name = " ".join(sys.argv[1:])  # Join all arguments into a single string

folder_sims='ics'
fullpath_snap = os.path.join(folder_sims, folder_name, 'ics')

folder_pylians='hdf5'
fullpath_hdf5 = os.path.join(folder_pylians, folder_name+'.hdf5')

pattern = r'_(\d+)_eta'
match = re.search(pattern, folder_name)
grid=int(match.group(1))

with open('NGENIC-PMFs_new/parameterfiles/'+folder_name+'.param', 'r') as file:
    # Read the entire file contents
    text = file.read()
pattern = r"\s+UnitLength_in_cm\s+([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s+%"
match = re.search(pattern, text)
UnitLength_in_cm = float(match.group(1))

pattern = rf'UnitMass_in_g\s+(\d*\.?\d+(?:[eE][-+]?\d+)?)'
match = re.search(pattern, text)
UnitMass_in_g = float(match.group(1))

pattern = rf'UnitVelocity_in_cm_per_s\s+(\d*\.?\d+(?:[eE][-+]?\d+)?)'
match = re.search(pattern, text)
UnitVelocity_in_cm_per_s = float(match.group(1))

pattern = rf'HubbleParam\s+(\d*\.?\d+(?:[eE][-+]?\d+)?)'
match = re.search(pattern, text)
h = float(match.group(1))

Gauss_to_internal = np.sqrt(UnitLength_in_cm**3)/(np.sqrt(UnitMass_in_g) * UnitVelocity_in_cm_per_s)/h

if grid==64:
    gfile = ReadGadget(fullpath_snap, fformat=0, longids=False, discard=['u'], opt_fields=['u', 'bfld'],  stop_at='', quiet=False, fill_mass=False, fixVelocity=False, verbose=False, swapbytes=False)
else:
    gfile = ReadGadgetMultipleFiles(fullpath_snap, fformat=0, longids=False, discard=['u'], opt_fields=['u', 'bfld'],  stop_at='', quiet=False, fill_mass=False, fixVelocity=False, verbose=False, swapbytes=False, max_num_files=0)

gfile.bfld=gfile.bfld*Gauss_to_internal

WriteGadget(gfile, fullpath_hdf5, fformat=3, longids=False, discard=[], opt_fields=['bfld'], stop_at='', quiet=False, verbose=False)

import h5py
f = h5py.File(fullpath_hdf5, 'r+')

f['Header'].attrs.create('NumPart_Total_HighWord', f['Header'].attrs['NumPart_Total_HW'])
del f['Header'].attrs['NumPart_Total_HW']

f['Header'].attrs.create('Flag_DoublePrecision', 0) 
f.close()
#%%
import netCDF4 as nc4
import numpy as np
#from pymatgen.ext.matproj import MPRester
#from pymatgen.core import Lattice

import json
from os import system, chdir
mpid = "mp-22922"
q_dir  ='q-random-1000'
#%%
def generate_anaddb_input(mpid):
  nqpts = 2000
  inp = dict.fromkeys(['eivec', 'nph2l', 'nph1l', 'qph1l', 'qph2l'])
  inp['eivec'] = 1
  inp['nph1l'] = nqpts+1
  inp['nph2l'] = inp['nph1l']

  with open('LiegeDataset/Repository/phonon/%s.json' %mpid) as f:
    data = json.load(f)
    print(data['metadata'])
    nkptx, nkpty, nkptz = data['metadata']['kpoints_grid']

#  vec_reduced, vec_cartesian = random_q_grid(mpid,nqpts) # Gamma is not a random point, so only nqpts-1 are called
  vec_reduced, vec_cartesian = random_q_grid(mpid,nqpts)

  inp['qph1l'] = 0.0
  inp['qph2l'] = inp['qph1l']

# I did not figure out an easy way to write a dictionary to a text file without going through lots of changes, so for now it just converts everything to strings and writes them down

  s = ''
  system('mkdir -p %s/%s' %(q_dir,mpid))
  with open('%s/%s/anaddb-%s.in' %(q_dir,mpid,mpid), 'w') as f:
    f.write('''eivec 1\nnph1l %i\nnph2l %i\nifcflag 1\nasr 2\nchneut 1\nngqpt %i %i %i\n''' %(nqpts+1, nqpts+1, nkptx, nkpty, nkptz))
    for j,z in zip([1,2],[1.0,0.0]):
      s += '\nqph%il ' %j
      if j == 1:
        vec_s = vec_reduced
      else:
        vec_s = vec_cartesian
      for i in range(len(vec_s)):
         s += ' '.join(map(str, vec_s[i])) + ' %f\n' %z
    f.write(s)
#%%
os.chdir('/Users/frost/Library/CloudStorage/OneDrive-ImperialCollegeLondon/My Stuff/Imperial College/Year 5/Polaron/HighThroughputPolarons/')

generate_anaddb_input(mpid)
# %%
# %%

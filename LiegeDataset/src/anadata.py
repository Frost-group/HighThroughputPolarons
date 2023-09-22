import numpy as np
#from pymatgen.ext.matproj import MPRester
#from pymatgen.core import Lattice


import functions as fn
import alps_units as un
import json
from os import system, chdir
import netCDF4 as nc4


q_norm = 1.e-4
q_dir  ='q-random-1000'

def lattice(mpid):
# cell = []

# with open('phonon/%s.json' %mpid) as f:
#   data = json.load(f)
#   for i in ['_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma']:
#     j = data['metadata']['structure'].split().index(i)
#     cell.append(float(data['metadata']['structure'].split()[j+1]))

# print ("cell ", cell)
# lattice = Lattice.from_parameters(a = cell[0], b = cell[1], c = cell[2], alpha = cell[3], beta = cell[4], gamma = cell[5])
# direct_lattice = lattice.__dict__['_matrix']
# print ("direct_lattice ", direct_lattice)

  f = nc4.Dataset('q-random-1000/%s/anaddb-%s.out_PHBST.nc' %(mpid,mpid),'r')
  direct_lattice = f.variables['primitive_vectors']   #phonon frequencies are in eV, and I need them in Ha

  reciprocal_lattice = np.linalg.inv(direct_lattice)
# print ("reciprocal_lattice ", reciprocal_lattice)

  return direct_lattice, reciprocal_lattice

def random_q_grid(mpid,nqpts):
  
  vec_cartesian = np.zeros((nqpts + 1,3))

  for i in range(1,nqpts+1):
    x = np.random.normal(0,1)
    y = np.random.normal(0,1)
    z = np.random.normal(0,1)
    r = np.sqrt(x**2 + y**2 + z**2)
    vec_cartesian[i] = q_norm*np.array([x,y,z])/r
  
  direct_lattice, reciprocal_lattice = lattice(mpid)

  vec_reduced = np.array([direct_lattice@q for q in vec_cartesian])

  return vec_reduced, vec_cartesian

def regular_q_grid(mpid,nqpts):

  vec_cartesian = []#np.zeros((nqpts + 1,3))

  dtheta = np.pi/100.0
  dphi = 2*np.pi/100.0

  for i in range(-50,51):
    for j in range(0,100):
      x = q_norm*np.cos(i*dtheta)*np.cos(j*dphi)
      y = q_norm*np.cos(i*dtheta)*np.cos(j*dphi)
      z = q_norm*np.sin(i*dtheta)
      vec_cartesian.append([x,y,z])
  
  vec_cartesian = np.array(vec_cartesian)
  direct_lattice, reciprocal_lattice = lattice(mpid)

  vec_reduced = np.array([direct_lattice@q for q in vec_cartesian])

  return vec_reduced, vec_cartesian

def generate_q_grid(mpid):
  a = np.sqrt(3)/3.
  vec_cartesian = q_norm*np.array([[0.0,0.0,0.0],
                             [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0], [-1.0,0.0,0.0], [0.0,-1.0,0.0], [0.0,0.0,-1.0],
                             [a,a,a], [-a,a,a], [a,-a,a], [a,a,-a], [-a,-a,a], [-a,a,-a], [a,-a,-a], [-a,-a,-a]])

  direct_lattice, reciprocal_lattice = lattice(mpid)

  vec_reduced = np.array([direct_lattice@q for q in vec_cartesian])

  vec_reduced = vec_reduced.tolist()
  vec_cartesian = vec_cartesian.tolist()

  return vec_reduced, vec_cartesian

def generate_anaddb_input(mpid):
  nqpts = 2000
  inp = dict.fromkeys(['eivec', 'nph2l', 'nph1l', 'qph1l', 'qph2l'])
  inp['eivec'] = 1
  inp['nph1l'] = nqpts+1
  inp['nph2l'] = inp['nph1l']

  with open('phonon/%s.json' %mpid) as f:
    data = json.load(f)
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

def run_anaddb(mpid):
  chdir('%s/%s' %(q_dir,mpid))
  with open('run_files.files', 'w') as f:
    f.write("""anaddb-%s.in
anaddb-%s.out
/home/pmelo/codes/alps/ddbs/%s_DDB
anaddb-%s-0
anaddb-%s-1
anaddb-%s-2
anaddb-%s-3
""" %(mpid,mpid,mpid,mpid,mpid,mpid,mpid))
  system("anaddb < run_files.files > anaddb-%s.log" %mpid)
  chdir('../../')

def read_eigendisplacements(mpid):
  # 1 Bohr, in Angstrom
  Bohr_Ang = 0.52917720859
  # 1 Angstrom in Bohr
  Ang_Bohr = 1.0/Bohr_Ang

  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  eigen = f.variables['phdispl_cart'][:]*Ang_Bohr       #displacements are in Angstrom and I need them in Bohr
  nqpts = f.dimensions['number_of_qpoints'].size
  modes = f.dimensions['number_of_phonon_modes'].size
  atoms = f.dimensions['number_of_atoms'].size

# print(eigen)

  store = np.zeros((modes,atoms,3,nqpts))

  store_img = np.zeros((modes,atoms,3,nqpts))

  for i in range(modes):
    for j in range(atoms):
      for k in range(3):
        for v in range(nqpts):
          store[i][j][k][v] = eigen[v][i][3*j+k][0]
          store_img[i][j][k][v] = eigen[v][i][3*j+k][1]

# print(store)

  grid_reduced, grid_cart = generate_q_grid(mpid)

  grid_cart = np.array(grid_cart)

  mat_to_inv = []
  vec_source = np.zeros((modes,atoms,3,10))
  vec_img = np.zeros((modes,atoms,3,10))

  coeff_mat = []
  coeff_mat_img = []

  for i in range(1,11):
    l = [1.0, round(grid_cart[i][0],10), round(grid_cart[i][1],10), round(grid_cart[i][2],10), round(grid_cart[i][0]*grid_cart[i][0],10), round(grid_cart[i][1]*grid_cart[i][1],10), round(grid_cart[i][2]*grid_cart[i][2],10), round(grid_cart[i][0]*grid_cart[i][1],10), round(grid_cart[i][0]*grid_cart[i][2],10), round(grid_cart[i][1]*grid_cart[i][2],10)]
    mat_to_inv.append(l)
  
  mat_to_inv = np.array(mat_to_inv)
# print(mat_to_inv)
 
  error = []
  error_img = []
 
  for i in range(modes):
    for j in range(atoms):
      for k in range(3):
        vec_source[i][j][k] = np.round(store[i][j][k][1:11],6)
        coefficients = np.round(np.linalg.solve(mat_to_inv, vec_source[i][j][k]),6)
        coeff_mat.append(coefficients.tolist())
        a = mat_to_inv@coefficients
        error.append(np.sum([np.dot(x,x) for x in a - vec_source[i][j][k]]))

        vec_img[i][j][k] = np.round(store_img[i][j][k][1:11],6)
        coeff_img = np.round(np.linalg.solve(mat_to_inv, vec_img[i][j][k]),6)
        coeff_mat_img.append(coeff_img)

        b = mat_to_inv@coeff_img
        error_img.append(np.sum([np.dot(x,x) for x in b - vec_img[i][j][k]]))
        

  error = np.sqrt(np.sum(error)/atoms/3.0/modes/10)
  error_img = np.sqrt(np.sum(error_img)/atoms/3.0/modes/10)

  coeff_mat = np.reshape(np.array(coeff_mat),(modes,atoms,3,10))
  coeff_mat_img = np.reshape(np.array(coeff_mat_img),(modes,atoms,3,10))
 
  return coeff_mat, coeff_mat_img, error, error_img, modes, nqpts, atoms

def read_omega_q(mpid):
  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  eigen = f.variables['phfreqs'][:]   #phonon frequencies are in eV, and I need them in Ha
  modes = f.dimensions['number_of_phonon_modes'].size

  Ha_eV = 27.21138282600827
  eigen = eigen.T/Ha_eV
  print('phonon_energies', eigen)

  grid_reduced, grid_cart = generate_q_grid(mpid)

  grid_cart = np.array(grid_cart)

  mat_to_inv = []
  vec_source = np.zeros((modes,10))

  coeff_mat = []

  for i in range(1,11):
    l = [1.0, round(grid_cart[i][0],10), round(grid_cart[i][1],10), round(grid_cart[i][2],10), round(grid_cart[i][0]*grid_cart[i][0],10), round(grid_cart[i][1]*grid_cart[i][1],10), round(grid_cart[i][2]*grid_cart[i][2],10), round(grid_cart[i][0]*grid_cart[i][1],10), round(grid_cart[i][0]*grid_cart[i][2],10), round(grid_cart[i][1]*grid_cart[i][2],10)]
    mat_to_inv.append(l)

  mat_to_inv = np.array(mat_to_inv)
  error = []
  
  for n in range(modes):
    coefficients = np.round(np.linalg.solve(mat_to_inv, eigen[n][1:11]),8)
    coeff_mat.append(coefficients.tolist())
    a = mat_to_inv@coefficients
    error.append(np.sum([np.dot(x,x) for x in a-eigen[n][1:11]]))

  error = np.sqrt(np.sum(error)/10/modes)
# print(np.sqrt(np.sum(error)/10/modes))

  return coeff_mat, error

def get_qpt_grid_cart(mpid):
  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  qpts  = f.variables['qpoints'][1:]

  direct_lattice = f.variables['primitive_vectors']
  reciprocal_lattice = np.linalg.inv(direct_lattice)

  qpts_cart = np.array([reciprocal_lattice@q for q in qpts])

  return qpts_cart


def test_orthogonality(mpid):
  # 1 Bohr, in Angstrom
  Bohr_Ang = 0.52917720859
  # 1 Angstrom in Bohr
  Ang_Bohr = 1.0/Bohr_Ang

  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  displ = f.variables['phdispl_cart'][1:]*Ang_Bohr       #displacements are in Angstrom and I need them in Bohr
  nqpts = f.dimensions['number_of_qpoints'].size - 1
  modes = f.dimensions['number_of_phonon_modes'].size
  atoms = f.dimensions['number_of_atoms'].size

  becs = fn.get_becs(mpid)
 
  qpts = get_qpt_grid_cart(mpid) 
  for iq in range(0,10):
    for nmode in range(modes):
      u = 0.0
      pol = np.zeros(3)
      for n in range(atoms):
        for i in range(3):
          for j in range(3):
            pol[j] += becs[n][i][j] * displ[iq][nmode][3*n+i][0]
            #u += np.round(np.round(becs[n][i][j],7)*np.round(displ[iq][nmode][3*n+i],7)*np.round(qpts[iq][j],7),7)
            u += becs[n][i][j] * displ[iq][nmode][3*n+i][0] * qpts[iq][j]
      print(iq, nmode,  u)


def get_q_weights(mpid):
  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  nqpts = f.dimensions['number_of_qpoints'].size - 1

  weight = np.zeros(nqpts) #The q-point list in the databases will contain Gamma, and we do not want it

  qpts = get_qpt_grid_cart(mpid)

  for iq in range(nqpts):
    weight[iq] = np.sqrt(qpts[iq][0]**2 + qpts[iq][1]**2)/q_norm
  return weight, nqpts


def get_q_dot_p_squared(mpid):
  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  eigen = f.variables['phdispl_cart'][1:]*un.Ang_Bohr       # Displacements are in Angstrom and I need them in Bohr
                                                            # Gamma point is at eigen[0], and it cannot be in the final sum
  nqpts = f.dimensions['number_of_qpoints'].size - 1
  modes = f.dimensions['number_of_phonon_modes'].size
  atoms = f.dimensions['number_of_atoms'].size

  becs = fn.get_becs(mpid)

  qpts = get_qpt_grid_cart(mpid)

  store = np.zeros((modes,nqpts))

  for q in range(0,nqpts):
    for nmode in range(modes):
      tmp = np.zeros((2))
      for j in range(3):
        for n in range(atoms):
          for i in range(3):
            tmp += becs[n][i][j]*eigen[q][nmode][3*n+i]*qpts[q][j]
      store[nmode][q] = (tmp[0]**2 + tmp[1]**2)/q_norm**2

  return store

def get_omega_q(mpid):
  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  eigen = f.variables['phfreqs'][1:]   #phonon frequencies are in eV, and I need them in Ha
  modes = f.dimensions['number_of_phonon_modes'].size

  Ha_eV = 27.21138282600827
  eigen = eigen/Ha_eV

  return eigen, modes

def get_masses(mpid):
  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  masses = f.variables['atomic_mass_units'][:] 

  return masses

def get_typat(mpid):
  f = nc4.Dataset('%s/%s/anaddb-%s.out_PHBST.nc' %(q_dir,mpid,mpid),'r')
  typat = f.variables['atom_species'][:] 

  return typat



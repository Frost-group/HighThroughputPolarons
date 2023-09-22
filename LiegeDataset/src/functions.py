#!/usr/bin/env python3
#  -*- coding: utf-8 -*-
"""
Routine: Functions to call when computing dielectric tensors

Author: Pedro Melo
Email: pmmelo@uliege.be

Final Purpose: provide functions to main routine - input file generation, job submission, etc 
"""
###############################################################################
# Import modules that already exist in python

from __future__ import print_function, division, unicode_literals

import sys
import json
from numpy import cos, sin, sqrt, pi
import numpy as np
import alps_units as un

import anadata as an


#from pymatgen.ext.matproj import MPRester



import scipy.integrate as integrate
import scipy.special as special

import matplotlib.pyplot as plt

#
# Check fucntions for epslion (electronic and total), effective mass, and phonon LO frequency. Note that they may or may not exist in different databases,
# so different checks are needed. Right now it works only with mp-ids
#

#
# right now check_epsilon prefers data from materials project. this is because I cannot understand the units from phonon database
# and how to convert them.
#

band_type = ['valence','conduction']

def get_formula(mpid):
  data  ="" # m.get_data(mpid)
  try:
    formula = data[0]['pretty_formula']
  except:
    with open('phonon/%s.json' %mpid, 'r') as g:
      guido_data = json.load(g)
      formula = guido_data['metadata']['formula']
  return formula

def get_atomic_info():
    with open('atoms/atomic_data.json', 'r') as g:
        data = json.load(g)
    return data


def check_epsilon(mpid):
# data  = m.get_data(mpid)
# try:
#   epsilon_t = data[0]['diel']['e_total']
#   epsilon_e = data[0]['diel']['e_electronic']
#   epsilon_s = 'mp-data'
# except:
  try:
    with open('phonon/%s.json' %mpid, 'r') as f:
      guido_data = json.load(f)
      epsilon_t = guido_data['dielectric']['eps_total']
      epsilon_e = guido_data['dielectric']['eps_electronic']
      epsilon_s = 'guido'
      if epsilon_t is None or epsilon_e is None:
          try:
              data_mp = m.get_data(mpid, prop="diel")[0]
              epsilon_e = data_mp["diel"]["e_electronic"]
              epsilon_t = data_mp["diel"]["e_total"]
              epsilon_s ="mp"
          except:
              print("%s not in any dielectric database" %mpid)
              epsilon_t = -1.0
              epsilon_e = -1.0
              epsilon_s = False
  except:
    print("%s not in any dielectric database" %mpid)
    epsilon_t = -1.0#np.array([1.0,1.0,1.0])
    epsilon_e = -1.0#np.array([1.0,1.0,1.0])
    epsilon_s = False
  return epsilon_t, epsilon_e, epsilon_s

def get_becs(mpid):
  try:
    with open('phonon/%s.json' %mpid,'r') as f:
      guido_data = json.load(f)
      becs = np.array(guido_data['dielectric']['becs'])
  except:
    becs = -1.0
    print('No Born effective charges for %s' %mpid)
  return becs

def ang_int(beta):
  result = integrate.quad(lambda x: np.cos(x)/sqrt((np.sin(x))**2 + (beta*np.cos(x))**2), -np.pi*0.5, np.pi*0.5)
  return result[0]*0.5

def check_mstar(mpid):
# try:
#   project = 'carrier_transport'
#   client = Client('zSHgsQzpW5sDel8WBhcX50OS2ziufHJM')

#   m_xx = '\u03b5\u2081'
#   m_yy = '\u03b5\u2082'
#   m_zz = '\u03b5\u2083'
#   m_econd = 'm\u2091\u1d9c\u1d52\u207f\u1d48'

#   u = client.contributions.get_entries(project=project,identifier=mpid,_fields=["data"]).result()

# conduction band
#   m_star_x = u['data'][0]['data'][m_econd]['n'][m_xx]['value']
#   m_star_y = u['data'][0]['data'][m_econd]['n'][m_yy]['value']
#   m_star_z = u['data'][0]['data'][m_econd]['n'][m_zz]['value']
# valence bands
#   m_star_x = u['data'][0]['data'][m_econd]['p'][m_xx]['value']
#   m_star_y = u['data'][0]['data'][m_econd]['p'][m_yy]['value']
#   m_star_z = u['data'][0]['data'][m_econd]['p'][m_zz]['value']

#   beta = m_star_x/m_star_z
#   m_star = m_star_x*ang_int(beta)
#   m_star_s = 'guido'
# except:
  try:
    with open('eff_masses/%s.json' %mpid,'r') as f:
      data = json.load(f)
      m_star_x, m_star_y, m_star_z = data['cond_eff_mass'][0]['p']['300']['1e+18']
      m_star = m_star_x*ang_int(m_star_x/m_star_z)
      m_star_s = 'joao'
  except:
    m_star = -1.0
    m_star_s = False
    print("%s not in any effective mass database" % mpid)
  return m_star, m_star_s


def read_abinit_effmasses(filename):
  effmass_flag = False
  val_idx = -1
  cond_idx = -1
  flag_occ = True
  val_tensor_deg =  np.zeros((3,3))
  val_eigen_deg =  np.zeros((3))
  cond_tensor_deg =  np.zeros((3,3))
  cond_eigen_deg =  np.zeros((3))
  with open(filename,"r") as f:
    for lines in f:
      if "  occ  " in lines and flag_occ:
        is_number = True
        finished = False
        l_split = lines.split()[1:]
        val_idx = 0
        while is_number and not finished:
          for ls in l_split:
            if float(ls) > 0.0:
              val_idx+=1
            else:
              finished = True
              break
          lines = next(f)
          l_split = lines.split()
          try:
            is_number = bool(float(l_split[0]))
          except:
            is_number = False
        cond_idx=val_idx+1
      if "CALCULATION OF EFFECTIVE MASSES" in lines:
        effmass_flag = True

      if "END OF EFFECTIVE MASSES SECTION" in lines:
        effmass_flag = False

      if effmass_flag:
        if "At k-point" in lines:
          l_split = lines.split()
          deg_bands = np.arange(int(l_split[-3]), int(l_split[-1])+1)
          if val_idx in deg_bands:
            val_tensor = np.zeros((len(deg_bands),3,3))
            val_eigen = np.zeros((len(deg_bands),3))
            val_bands = deg_bands
            cond_bands = []
          elif cond_idx in deg_bands:
            cond_tensor = np.zeros((len(deg_bands),3,3))
            cond_eigen = np.zeros((len(deg_bands),3))
            cond_bands = deg_bands

        if "K-point" in lines:
          l_split = lines.split()
          band = int(l_split[-1])

          if val_idx == band and "val_tensor" not in locals():
            val_tensor = np.zeros((1,3,3))
            val_eigen = np.zeros((1,3))
            val_bands = np.array([band])
            cond_bands = []
          if cond_idx == band and "cond_tensor" not in locals():
            cond_tensor = np.zeros((1,3,3))
            cond_eigen = np.zeros((1,3))
            cond_bands = np.array([band])

        if "mass tensor" in lines and ":" in lines: 
          if "eigenvalues" in lines:
            mdir = 1
          else:
            mdir = 3
          for idir in range(mdir):
            lines = next(f)
            l_split = lines.split()
            if band in val_bands:
              idx = np.where(val_bands == band)[0][0]
              if "SADDLE" in lines:
                val_tensor[idx,:,:] = [None, None, None]
                val_eigen[idx,:] = [None, None, None]
                break
              if mdir == 3:
                  val_tensor[idx,idir,:] = np.array([float(l_split[0]), float(l_split[1]), float(l_split[2])])
              elif mdir == 1:
                  val_eigen[idx,:] = np.array([float(l_split[0]), float(l_split[1]), float(l_split[2])])
            elif band in cond_bands:
              idx = np.where(cond_bands == band)[0][0]
              if "SADDLE" in lines:
                cond_tensor[idx,:,:] = [None, None, None]
                cond_eigen[idx,:] = [None, None, None]
                break
              if mdir == 3:
                  cond_tensor[idx,idir,:] = np.array([float(l_split[0]), float(l_split[1]), float(l_split[2])])
              elif mdir == 1:
                  cond_eigen[idx,:] = np.array([float(l_split[0]), float(l_split[1]), float(l_split[2])])

# if len(val_bands) == 1:
#   val_tensor_deg[:,:] = val_tensor[0,:,:]
#   val_eigen_deg[:] = val_eigen[0,:]
# if len(cond_bands) == 1:
#   cond_tensor_deg[:,:] = cond_tensor[0,:,:]
#   cond_eigen_deg[:] = cond_eigen[0,:]



# for i in range(3):
#   if len(val_bands) > 1:
#     val_eigen_deg[i] = np.sum(np.abs(val_eigen[:,i])**0.5) / len(val_bands)
#   if len(cond_bands) > 1:
#     cond_eigen_deg[i] = np.sum(np.abs(cond_eigen[:,i])**0.5) / len(cond_bands)
#   for j in range(3):
#     if len(val_bands) > 1:
#       val_tensor_deg[i,j] = np.sum(np.abs(val_tensor[:,i,j])**0.5) / len(val_bands)

#     if len(cond_bands) > 1:
#       cond_tensor_deg[i,j] = np.sum(np.abs(cond_tensor[:,i,j])**0.5) / len(cond_bands)
  
  return [ [ -1*val_tensor, np.abs(val_eigen), val_bands ], [ cond_tensor, cond_eigen, cond_bands ] ]
#  return [[np.abs(val_eigen), val_bands], [np.abs(cond_eigen), cond_bands]]


def get_mstar(mpid, opt, what):
  m_tensor = 0.0*np.ones((3,3))
  m_eigen = [ [-1.0,-1.0,-1.0] ]
  m_deg = [-1.0]
  m_star_s = False
  NotWorking = True

  if opt["eff_masses"] == "mp" or opt["eff_masses"] == "all":
    try:
      with open('eff_masses/%s.json' %mpid,'r') as f:
        data = json.load(f)
        if what == 1:
            m_star_x, m_star_y, m_star_z = data['cond_eff_mass'][0]['n']['300']['1e+18']
        else:
            m_star_x, m_star_y, m_star_z = data['cond_eff_mass'][0]['p']['300']['1e+18']
        m_eigen = np.abs([[m_star_x, m_star_y, m_star_z]])
        m_tensor = np.array([ m_tensor] )
        m_deg = [1]
        m_star_s = 'mp'
        NotWorking = False
    except:
      pass

  if opt["eff_masses"] == "abinit" or (opt["eff_masses"] == "all" and NotWorking):
    try:
      
      #tensor = read_abinit_effmasses("abinit_effmasses/%s_output" %mpid)[1][0]
      # While tensor is not applicable in the integration of the q-grid get eigenvalues
      #m_star_x, m_star_y, m_star_z = read_abinit_effmasses("abinit_effmasses/%s_output" %mpid)[0][0][0]
      m_tensor, m_eigen, m_deg = read_abinit_effmasses("abinit_effmasses/%s_output" %mpid)[what]
      #m_star_y  =1.0
      #m_star_z = 1.0
      m_star_s = "abinit"
      NotWorking = False
    except:
      pass

  if NotWorking:
    print("%s not in any effective mass database" % mpid)  
    return -1
  else:
    print("%s effective mass calculated from %s database" % (mpid, m_star_s))
    return m_tensor, m_eigen, m_deg, m_star_s

#
# COMMENT: WHILE MATHEMATICALLY EQUIVALENT, THE 1/SQRT(1-X^2) MAKES THIS INTEGRAL UNSTABLE
#          IT IS PREFERABLE TO USE THE ANGULAR FORM
#def angular_int(beta):
#  result = integrate.quad(lambda x: x/(sqrt(1.0 - x**2)*sqrt(x**2 + beta**2*(1.0-x**2))), 0.0, 0.9999999999999999)
#  return result[0]

#def angular_m(mx,my,mz,y,x):
#  b = mx/mz
#  a = mx/my
#  return np.sqrt(mx)/np.sqrt((np.sin(y)*np.cos(x))**2 + a*(np.sin(y)*np.sin(x))**2 + b*(np.cos(y))**2)*cos(y)

def angular_m(m_tensor,y,x):
  m1 = np.linalg.inv(m_tensor)
  #q = np.array([cos(x)*sin(y), sin(x)*sin(y),cos(y)])
  q = np.array([cos(x)*cos(y), sin(x)*cos(y),sin(y)])
  return 1.0/np.sqrt(q@m1@q)*cos(y)


def m2integrate(m_tensor):
  #print(mx,my,mz)
  #f = lambda y,x: angular_m(mx,my,mz,y,x)
  f = lambda y,x: angular_m(m_tensor,y,x)
  #print(f(0,2*np.pi))
  m_star_d, error = integrate.dblquad(f, 0, 2*np.pi, lambda x: -np.pi/2.0, lambda x: np.pi/2.0,epsabs=1e-10)
  #m_star_d, error = integrate.dblquad(f, 0, 2*np.pi, lambda x: 0.0, lambda x: np.pi,epsabs=1e-10)
  print("Any error during angular integration? ", error)
  return m_star_d

    
crystal_groups = [{ 'triclinic'    : {'triclinic', '1', '-1'},
                    'monoclinic'   : {'monoclinic', '2', 'm', '2/m'},
                    'orthorhombic' : {'orthorhombic', '222', 'mm2', 'mmm'},
                    'tetragonal'   : {'tetragonal', '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm'},
                    'trigonal'     : {'trigonal', '3', '-3', '32', '3m', '-3m'},
                    'hexagonal'    : {'hexagonal', '6', '-6', '6/m', '622', '6mm', '-6m2', '6/mmm'},
                    'cubic'        : {'23','m-3','432','-43m','m-3m','cubic'}}]

# OLD FUNCTION WHICH READS THE EIGENDISPACEMENTS FROM A DATABASE MADE FROM PARSING ANADDB OUTPUTS - TO BE EDITED TO READ THEM FROM THE NEW DATABASES
def get_eigendisplacements(mpid):
  natoms = get_num_atoms(mpid)
  #with open('phonon/%s.json' %mpid) as f:
  #  data = json.load(f)
  #  natoms = data['metadata']['nsites']
  nmodes = 3*natoms
  block = 1+natoms*2

  eigendisplacements = []

  with open('eigen-displacements/eigendisplacements-%s_DDB' %mpid) as f:
    data_raw = f.readlines()
    for i in range(1,nmodes*(1+2*natoms)+1,1+2*natoms):
      for j in range(0,natoms*2,2):
        z = [complex(float(x),float(y)) for x in data_raw[i+j].split()[1:] for y in data_raw[i+1+j].split()]
        eigendisplacements.append([float(x) for x in data_raw[i+j].split()[1:]])

  u = np.reshape(eigendisplacements,(nmodes,natoms,3))
  return u


def get_omega_gamma(mpid):
  try:
    with open('phonon/%s.json' %mpid) as f:
      data = json.load(f)
      omega_gamma = [x*un.cmm12meV*0.001*un.eV2Ha for x in data['phonon']['ph_bandstructure'][0]]      #phonon frequencies are in cm^-1 in the Phonon database and we need them in a.u., i.e, Ha
      omega_s = True
  except:
    omega_gamma = -1.0
    omega_s = False
  return omega_gamma, omega_s

def get_num_atoms(mpid):
  with open('phonon/%s.json' %mpid, 'r') as f:
    data = json.load(f)
    return data['metadata']['nsites']

def get_cell_volume(mpid):
  with open('phonon/%s.json' %mpid) as f:
    data = json.load(f)
    where = data['metadata']['structure'].split().index('_cell_volume')
    return float(data['metadata']['structure'].split()[where+1])*(un.Ang_Bohr)**3  # convert from cubic angstrom to cubic Bohr

def get_polarity_vector_gamma(mpid):
  try:
    natoms = get_num_atoms(mpid)
    #with open('phonon/%s.json' %mpid) as f:
    #  data = json.load(f)
    #
    #  natoms = data['metadata']['nsites']
    nmodes = 3*natoms
    block = 1+natoms*2
  
    becs = get_becs(mpid)
    eigendisplacements = get_eigendisplacements(mpid)
        
    u = []

    for nmode in range(nmodes):
      z = 0
      for natom in range(natoms):
        z += np.matmul(np.array(becs[natom]),np.array(eigendisplacements[nmode][natom].T))
      u.append(z)
    polarity_vec = u
    polarity_s = True
  except:
    polarity_vec = []
    polarity_s = False
  return polarity_vec, polarity_s

def get_polarity_vector_q(mpid,coeff_mat,coeff_mat_img,natoms,q_norm,y,x):
  becs = get_becs(mpid)
  
  q1 = np.array([q_norm*cos(x)*sin(y), q_norm*sin(x)*sin(y), q_norm*cos(y)])
  q2 = np.array([1.0, q1[0], q1[1], q1[2], q1[0]*q1[0], q1[1]*q1[1], q1[2]*q1[2], q1[0]*q1[1], q1[0]*q1[2], q1[1]*q1[2]]) 

  u = 0
  v = 0

  for natom in range(natoms):
   for i in range(0,3):
     for j in range(0,3):
       for n in range(10):
         u += becs[natom][j][i]*coeff_mat[natom][j][n]*q2[n]*q1[i]/q_norm
         v += becs[natom][j][i]*coeff_mat_img[natom][j][n]*q2[n]*q1[i]/q_norm

  return u**2 + v**2

def interpol(q_grid,quantity,point):
  return griddata(q_grid,quantity,point)
  
def angular_epsilon(epsilon,y,x):
  q = np.array([cos(x)*sin(y), sin(x)*sin(y),cos(y)])
  return q@epsilon@q.T

def interpol_omega(omega_coeff,q_norm,y,x):

  qx = q_norm*cos(x)*sin(y)
  qy = q_norm*sin(x)*sin(y)
  qz = q_norm*cos(y)
  q = np.array([1.0, qx, qy, qz, qx*qx, qy*qy, qz*qz, qx*qy, qx*qz, qy*qz]) 

  w = 0

  for i in range(10):
    w += q[i]*omega_coeff[i]

  return w

def mass_q(m_tensor,q,q_norm):
#  a = mx/my
#  b = mx/mz
  m_1 = np.linalg.inv(m_tensor)
  return  q_norm/np.sqrt(q@m_1@q)
# return  np.sqrt(mx)*q_norm/np.sqrt(q[0]**2 + a*(q[1])**2 + b*(q[2])**2)

def epsilon_q(epsilon,q,q_norm):
  return q@epsilon@q/q_norm**2

def generalised_zpr(mpid, opt, what):
  volume = get_cell_volume(mpid)
  #print("vol: ",volume)
  omega_gamma, os = get_omega_gamma(mpid)
  #print("omega_gamma", omega_gamma)
  try:
    #m_eigen, m_deg = read_abinit_effmasses("abinit_effmasses/%s_output" %mpid)[0]
    m_tensor, m_eigen, m_deg, m_s = get_mstar(mpid, opt, what)
    #print(mpid,m_deg)
    ms = True
  except:
    ms = False
  #print("mx,my,mz: ",mx,my,mz)
  epsilon_t, epsilon_e, epsilon_s = check_epsilon(mpid) 
  #print(epsilon_t,epsilon_e,epsilon_s)

  zpr = []
  average_omega = []
  average_epsm1 = []
  average_mstar = 0.0
  average_alpha = []

  if m_s and epsilon_s and os:
    if m_s == "mp":
        for i in range(3):
            m_tensor[0,i,i] = m_eigen[0,i]
    polarity = an.get_q_dot_p_squared(mpid)

    ph_energies, nmodes = an.get_omega_q(mpid)
 
    q_grid = an.get_qpt_grid_cart(mpid)

    weight, nqpts = an.get_q_weights(mpid)

    for d in range(len(m_deg)):
       mx, my, mz = m_eigen[d]
#      print('print mass tensor, degeneracies, eigenvalues', m_eigen, m_deg, mx, my, mz)
       #average_mstar += m2integrate(mx,my,mz)/4./np.pi/len(m_deg)
       average_mstar += m2integrate(m_tensor[d])/4./np.pi/len(m_deg)

    for nmode in range(0,nmodes):
      q_norm = 1.e-4 #norm of the small q-vector

      tmp_zpr = 0.0
      tmp_eps = 0.0
      tmp_alp = 0.0
      tmp_omg = 0.0
      check_sum = 0.0
      mass_deg = 0.0
      if omega_gamma[nmode] <= 0.0:
        pass 
      else:
        for iq in range(nqpts):
          tmp_mass = 0.0
          tmp = polarity[nmode][iq]/(ph_energies[iq][nmode]**(3./2.)*epsilon_q(epsilon_e,q_grid[iq],q_norm)**2)*weight[iq]*2*pi**2

#I tried to see if making the degeneracy=1 explicit would help, but it did not 
#         if len(m_deg) == 1:
#           mx, my, mz = m_eigen[0]
#           tmp_zpr += tmp*mass_q(mx,my,mz,q_grid[iq],q_norm)
#           tmp_alp += tmp/ph_energies[iq][nmode]*mass_q(mx,my,mz,q_grid[iq],q_norm)
#         else:
          for d in range(len(m_deg)):
            mx, my, mz = m_eigen[d]
  
            tmp_zpr += tmp*mass_q(m_tensor[d],q_grid[iq],q_norm)/len(m_deg)
            tmp_alp += tmp/ph_energies[iq][nmode]*mass_q(m_tensor[d],q_grid[iq],q_norm)/len(m_deg)
            tmp_mass += mass_q(m_tensor[d],q_grid[iq],q_norm)

#           print('intermediate quantities',mass_deg, mass_q(mx,my,mz,q_grid[iq],q_norm),d)
          mass_deg += tmp_mass/len(m_deg)*weight[iq]*2*pi**2
#         print('after degeneracy loop', mass_deg, len(m_deg), weight[iq], 2*pi**2)
          #print("mass_deg:",mass_deg)
          tmp_alp2 = tmp/ph_energies[iq][nmode]*mass_deg

          tmp_eps += tmp/np.sqrt(ph_energies[iq][nmode])
          tmp_omg += weight[iq]/np.sqrt(ph_energies[iq][nmode])
          check_sum += weight[iq]*2*pi**2

        mass_deg = mass_deg/(4.*np.pi*nqpts*len(m_deg))
      #print('mass_deg',mass_deg, 'check_sum', check_sum/nqpts/4./np.pi, 'average_mstar',average_mstar, flush=True)

      zpr.append(-tmp_zpr/(np.sqrt(2)*volume*un.eV2Ha)*1000/nqpts)
      average_epsm1.append(tmp_eps/(volume*nqpts))
      average_alpha.append(tmp_alp/(volume*sqrt(2)*nqpts))
      average_omega.append(tmp_omg/(4.*np.pi)*sqrt(0.001*un.eV2Ha)/nqpts)

      zpr_s = True

    with open('%s-data-per-mode-%s.dat' %(mpid,band_type[what]), 'w') as fp:
     fp.write("#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \%-20s\n" %('PHONON MODE', 'EPSM1', 'OMEGAM05[meV]', 'ALPHA', 'ZPR[meV]', 'AREA/4pi'))
     for i in range(nmodes):
       fp.write("%-20i \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f\n" %(i, average_epsm1[i], average_omega[i], average_alpha[i], zpr[i], check_sum/nqpts/4/pi))
    #print(sum(average_epsm1), average_mstar, sum(average_omega), sum(average_alpha), sum(zpr), zpr_s,check_sum/nqpts/4.0/pi, len(m_deg))
    return sum(average_epsm1), average_mstar, sum(average_omega), sum(average_alpha), sum(zpr), zpr_s, check_sum/nqpts/4.0/pi, len(m_deg)

  else:
    zpr_s = False
    zpr = 0.0
    average_omega =-1.0
    average_epsm1 =-1.0
    average_mstar =-1.0
    average_alpha =-1.0
    check_sum = -1.0
    return average_epsm1, average_mstar, average_omega, average_alpha, zpr, zpr_s, check_sum, len(m_deg)


def compare_epsilon(mpid):
  comp_file = open("epsilon_comparison-3.dat","a")
# data = m.get_data(mpid)
# try:
#   eps_t_mp = data[0]['diel']['e_total'][0][0]
#   eps_e_mp = data[0]['diel']['e_electronic'][0][0]
# except:
#   eps_t_mp = -1000
#   eps_e_mp = -1000
# with open('phonon/%s.json' %mpid, 'r') as f:
#   guido_data = json.load(f)
#   try:
#     eps_t_g = guido_data['dielectric']['eps_total'][0][0]
#     eps_e_g = guido_data['dielectric']['eps_electronic'][0][0]
#   except:
#     eps_t_g = -1000
#     eps_e_g = -1000 
# if (abs(eps_t_mp-eps_t_g) > 0.15*eps_t_mp or abs(eps_e_mp-eps_e_g) > 0.15*eps_e_mp) and (eps_t_mp != -1000 and eps_t_g != -1000) and (eps_e_mp != -1000 and eps_e_g != -1000):
#   bs = m.get_bandstructure_by_material_id(mpid)
#   try:
#     dir_band_gap = bs.get_direct_band_gap()
#   except:
#     dir_band_gap = -1.0
#   formula = get_formula(mpid)
#   m_star, m_s = check_mstar(mpid)
#   comp_file.write("%-15f \t %-15f \t %-15f \t %-15f \t %-15s \t %-15s \t %-15f \t %-15f\n" %(eps_t_mp, eps_t_g, eps_e_mp, eps_e_g, formula, mpid, dir_band_gap, m_star))
  comp_file.close()


def check_wlo(mpid):
  try:
    with open('phonon/%s.json' %mpid, 'r') as f:
      guido_data = json.load(f)
      w_lo = guido_data['phonon']['ph_bandstructure'][0][-1]*un.cmm12meV
      w_lo_s = 'guido'
  except:
    print("%s not in any phonon database" %mpid)
    w_lo = -1.0
    w_lo_s = False
  return w_lo, w_lo_s

def check_wlo_hellwarth(mpid, eps_star, eps_inf, atominfo):
  #try:
  if True:
    #with open('phonon/%s.json' %mpid, 'r') as f:
# from json file Guido made
    #guido_data = json.load(f)

    volume = get_cell_volume(mpid)

    natoms = get_num_atoms(mpid)
      #guido_data['metadata']['nsites']

    becs = get_becs(mpid)

# from netcdf file anaddb produced
    typat = np.array(atominfo["atom_species"])

    masses = np.array(atominfo["atomic_mass_units"]) * un.amu2me # Convert from amu to electron mass 



     
# precalculate sum of born charges squared, divided by mass
    wlo_hellwarth = 0.0
    for iatom in range(natoms):
      wlo_hellwarth += sum(sum(becs[iatom][:][:]**2)) / masses[typat[iatom]-1]


# w_hell**2 =  eps* / eps_inf**2 * 4 pi / 3 Volume * 
#              [sum_kappa  1/M_kappa * sum_alpha,beta Z*_kappa,alpha,beta**2 ] #from previous loop
    #wlo_hellwarth = eps_star / eps_inf**2 * .75 * np.pi  / volume * wlo_hellwarth 
    #PART1 = eps_star / eps_inf**2 * 4 * np.pi  / volume/3
    #PART2 = wlo_hellwarth
    #print("PART1: ", eps_star / eps_inf**2 * 4 * np.pi  / volume/3)
    #print("PART2: ", wlo_hellwarth)
    #print("PART1/PART2: ",PART1/PART2)
    wlo_hellwarth = eps_star / eps_inf**2 * 4 * np.pi  / volume * wlo_hellwarth  /3
    wlo_hellwarth = wlo_hellwarth**0.5 * un.Ha2eV * 1000 # Hartree to eV and eV to meV 
    
    wlo_hellwarth_s = "guido"
    

  #except:
  #  print("%s not in any phonon database" %mpid)
  #  wlo_hellwarth = -1.0
  #  wlo_hellwarth_s = False
  return wlo_hellwarth, wlo_hellwarth_s #, PART1, PART2

#COMPUTES THE STANDARD FROELICH MODEL ALPHA AND ZPR 
def standard_zpr(mpid,opt,what):
  e_0, e_inf, e_s = check_epsilon(mpid)
  # TODO calculate and use the 1/epsilon average as it is done for Hellwarth
  #print('e_0, e_inf, e_s', e_0, e_inf, e_s)
  wlo, wlo_s = check_wlo(mpid)
  m_tensor, m_eigen, m_deg, m_s = get_mstar(mpid, opt, what)
  #m_eigen, m_deg = read_abinit_effmasses("abinit_effmasses/%s_output" %mpid)[0]
  m_star = 0.0

  if m_s:
    if m_s == "mp":
        for i in range(3):
            m_tensor[0,i,i] = m_eigen[0,i]
    for d in range(len(m_deg)):
#     print(m_eigen[d])
      mx, my, mz = m_eigen[d]
#      print(mx,my,mz)
      #print(m_tensor[d])
      #print(m_deg)
      m_star += m2integrate(m_tensor[d])/4/np.pi/len(m_deg)
    #m_star += m2integrate(mx,my,mz)/4/np.pi/len(m_deg)
#     print("m_star: ",m_star,m_deg)

  if e_s and wlo_s and e_0 and e_inf and wlo and m_star:
    #e_0 = np.trace(e_0)/3.
    #e_inf = np.trace(e_inf)/3.
    e_inf = 3./ np.trace(np.linalg.inv(e_inf)) #/ 3.
    e_0 = 3./np.trace(np.linalg.inv(e_0)) #/ 3.
    e_star = np.trace(np.linalg.inv(e_inf) - np.linalg.inv(e_0)) / 3.
    eps_star = 1./e_star

    alpha = (1.0/e_inf - 1.0/e_0)*np.sqrt(0.5/(wlo*0.001*0.0367493))*m_star
    zpr = -alpha*wlo
    zpr_s = True
  else:
    e_0 = -1.0
    e_inf = -1.0
#   wlo = -1.0
    m_star = -1.0
    alpha = -1.0 
    zpr = -1.0
    zpr_s = False
  return round(eps_star,4), round(wlo,4), round(m_star,4), round(alpha,4), round(zpr,4), zpr_s

#COMPUTES THE Hellwarth FROELICH MODEL ALPHA AND ZPR 
def hellwarth_zpr(mpid,opt,what, atominfo):

# extract dielectric tensor from mpid in database
  e_0, e_inf, e_s = check_epsilon(mpid)
  #print('e_0, e_inf, e_s', e_0, e_inf, e_s)
  if e_s and e_0 and e_inf:
    eps_inf = 3./ np.trace(np.linalg.inv(e_inf)) #/ 3.
    eps_0 = 3./np.trace(np.linalg.inv(e_0)) #/ 3.
    e_star = np.trace(np.linalg.inv(e_inf) - np.linalg.inv(e_0)) / 3.
    eps_star = 1./e_star
    #wlo, wlo_s, PART1, PART2 = check_wlo_hellwarth(mpid, eps_star, eps_inf)
    wlo, wlo_s = check_wlo_hellwarth(mpid, eps_star, eps_inf, atominfo)
  else:
    print("Problems in calculating epsilon in %s" %mpid) 

# extract average Hellwarth phonon frequency from mpid number

# effective mass from database?
  m_tensor, m_eigen, m_deg, m_s = get_mstar(mpid, opt, what)

  #m_eigen, m_deg = read_abinit_effmasses("abinit_effmasses/%s_output" %mpid)[0]
  m_star = 0.0

# average effective mass if needed
  if m_s:
    if m_s == "mp":
        for i in range(3):
            m_tensor[0,i,i] = m_eigen[0,i]
    for d in range(len(m_deg)):
      #print(m_eigen[d])
      
      mx, my, mz = m_eigen[d]
#     print(m_tensor)
      m_star += m2integrate(m_tensor[d])/4/np.pi/len(m_deg)
 # print(m_tensor)
 # print(m_eigen)
  #wlo = wlo * un.Ha2eV
    #m_star += m2integrate(mx,my,mz)/4/np.pi/len(m_deg)
#     print("m_star: ",m_star,m_deg)

# if we have all we need
  if e_star and wlo_s and eps_0 and eps_inf and wlo and m_star:
# average inverse of dielec constants, trace * 1/3

    alpha = e_star * np.sqrt(0.5/(wlo * 0.001 * 0.0367493)) * m_star
    zpr = -alpha*wlo
    zpr_s = True


  else:
    eps_0 = -1.0
    eps_inf = -1.0

#   wlo = -1.0
    m_star = -1.0
    alpha = -1.0 
    zpr = -1.0
    zpr_s = False
  #print(eps_0, eps_inf, wlo, m_star, alpha, zpr, zpr_s)

  #return round(eps_0,4), round(eps_inf,4), round(wlo,4), round(m_star,4), round(alpha,4), round(zpr,4), zpr_s#, round(PART1,4), round(PART2,4)
  return round(eps_star,4), round(wlo,4), round(m_star,4), round(alpha,4), round(zpr,4), zpr_s#, round(PART1,4), round(PART2,4)


#COMPUTES THE FEYNMAN PATH INTEGRAL FOR STANDARD FROELICH MODEL ALPHA AND ZPR 
def feynman_zpr(mpid,opt,what,atominfo):

# extract dielectric tensor from mpid in database
  e_0, e_inf, e_s = check_epsilon(mpid)
  #print('e_0, e_inf, e_s', e_0, e_inf, e_s)
  if e_s and e_0 and e_inf:
    eps_inf = 3./ np.trace(np.linalg.inv(e_inf)) #/ 3.
    eps_0 = 3./np.trace(np.linalg.inv(e_0)) #/ 3.
    e_star = np.trace(np.linalg.inv(e_inf) - np.linalg.inv(e_0)) / 3.
    eps_star = 1./e_star
    #wlo, wlo_s, PART1, PART2 = check_wlo_hellwarth(mpid, eps_star, eps_inf)
    wlo, wlo_s  = check_wlo_hellwarth(mpid, eps_star, eps_inf,atominfo)
  else:
    print("Problems in calculating epsilon in %s" %mpid) 

# extract average Hellwarth phonon frequency from mpid number

# effective mass from database?
  m_tensor, m_eigen, m_deg, m_s = get_mstar(mpid, opt, what)

  #m_eigen, m_deg = read_abinit_effmasses("abinit_effmasses/%s_output" %mpid)[0]
  m_star = 0.0

# average effective mass if needed
  if m_s:
    if m_s == "mp":
        for i in range(3):
            m_tensor[0,i,i] = m_eigen[0,i]

    for d in range(len(m_deg)):
      #print(m_eigen[d])
      
      mx, my, mz = m_eigen[d]
#     print(m_tensor)
      m_star += m2integrate(m_tensor[d])/4/np.pi/len(m_deg)
 # print(m_tensor)
 # print(m_eigen)
  #wlo = wlo * un.Ha2eV
    #m_star += m2integrate(mx,my,mz)/4/np.pi/len(m_deg)
#     print("m_star: ",m_star,m_deg)

# if we have all we need
  if e_star and wlo_s and eps_0 and eps_inf and wlo and m_star:
# average inverse of dielec constants, trace * 1/3

    
    #print("e star: ",e_star)
    #print("sqrt(0.5/wlo)",np.sqrt(0.5/(wlo * 0.001 * 0.0367493)))
    #print("m_star", m_star)
    alpha = e_star * np.sqrt(0.5/(wlo * 0.001 * 0.0367493)) * m_star
    #print("ALPHA,WLO")
    #print(alpha,wlo)
    if alpha <= 5:
        zpr = -wlo*(alpha + 0.98*(alpha/10)**2 + 0.60*(alpha/10)**3 + 0.14*(alpha/10)**4)
    else:
        #print("ALPHA**2", alpha**2, "WLO*ALPHA**2", wlo*alpha**2)
        zpr = -wlo*(0.106*alpha**2 + 2.83)
    zpr_s = True


  else:
    eps_0 = -1.0
    eps_inf = -1.0

#   wlo = -1.0
    m_star = -1.0
    alpha = -1.0 
    zpr = -1.0
    zpr_s = False
  #print(eps_0, eps_inf, wlo, m_star, alpha, zpr, zpr_s)

  return round(eps_star,4), round(wlo,4), round(m_star,4), round(alpha,4), round(zpr,4), zpr_s




def plot_alpha_zpr(x,y,names,c=[],size=[],color_ticks_values=[],color_ticks_labels=[],title='',x_label='',y_label='',hardcopy=False,xyline=False,yline=False,xmax=None,xmin=None,ymax=None,ymin=None,labels='',logx=False,logy=False):

  fig,ax = plt.subplots()
# sc = plt.scatter(x,y,c=c, s=100, cmap=cmap, norm=norm)
  
#  if c.any():
#    norm = plt.Normalize(min(c),max(c))
#    cmap = plt.cm.rainbow
#    sc = plt.scatter(x,y,c=c, s=100, cmap=cmap, norm=norm)
#  else:
  if any(size):
    sc = plt.scatter(x,y,c=c, s=size)
  else:
    sc = plt.scatter(x,y,c=c, s=100)
 
  if yline:
    plt.axhline(y=0, color='black', linestyle='--')

  if xmax and xmin:
    plt.xlim=((xmin,xmax))

  if logx:
    plt.xscale('log')

  if logy:
    plt.yscale('log')

  plt.xlim=((0.0,60.0))
  if ymax and ymin:
    plt.ylim=((ymin,ymax))

  ax1 = plt.gca()
  ax1.set_xlim([xmin, xmax])
  ax1.set_ylim([ymin, ymax])

  if xyline:
    if title:
      if 'ZPR' in title:
        linestart = np.min((min(x), min(y)))
        lineend = np.max((0.0, 0.0))
      if 'alpha' in title:
        linestart = np.min((0.00, 0.00))
        lineend = np.max((max(x),max(y)))

    extrema = np.linspace(linestart,lineend,2)
    ax.plot(extrema,extrema,color='black',ls='-')
    ax.plot(extrema,-extrema,color='black',ls='-')

#   ax.plot(extrema,1.5*extrema,color='black',ls='--')
#   ax.plot(extrema,0.5*extrema,color='black',ls='--')

#   ax.plot(extrema,-1.5*extrema,color='black',ls='--')
#   ax.plot(extrema,-0.5*extrema,color='black',ls='--')

#  if list(color_ticks_values):# and color_ticks_labels.any():
#    bar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),ticks=color_ticks_values)
#    bar.ax.set_yticklabels(color_ticks_labels)

  annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->"))
  annot.set_visible(False)
  
  if x_label:
    plt.xlabel(x_label)
  if y_label:
    plt.ylabel(y_label)
 
  labels_list = {'RbN3' : r'RbN$_3$',
         '         KN3' : r'KN$_3$', 
                 'NaN3' : r'NaN$_3$', 
                 'LiN3' : r'LiN$_3$', 
               'K2TiF6' : r'K$_2$TiF$_6$' , 
            'Cs2NaScF6' : r'Cs$_2$NaScF$_6$', 
                'CsNO2' : r'CSNO$_2$'}
#Missing Li2CaHfF6?

  labels_coord = []
  label2_coord = []

  if labels:
    for l in labels_list:
      try:
        xl = x[names.index(l)]
        yl = y[names.index(l)]
        labels_coord.append([xl,yl]) 
        xl = x[names.index(l,1035)]
        yl = y[names.index(l,1035)]
        label2_coord.append([xl,yl])

      except:
        pass

    if 'alpha' in title:
      plt.text(labels_coord[-1][0] + 5, labels_coord[-1][1] + 5, labels_list[list(labels_list)[-1]])
      plt.arrow(labels_coord[-1][0],labels_coord[-1][1],5,5)
    
      plt.text(labels_coord[-2][0] + 15, labels_coord[-2][1] + 5, labels_list[list(labels_list)[-2]])
      plt.arrow(labels_coord[-2][0],labels_coord[-2][1],15,5)

      plt.text(labels_coord[-3][0] + 20, labels_coord[-3][1] + 5, labels_list[list(labels_list)[-3]])
      plt.arrow(labels_coord[-3][0],labels_coord[-3][1],20,5)

      plt.text(labels_coord[-4][0] + 20, labels_coord[-4][1] - 2, labels_list[list(labels_list)[-4]])
      plt.arrow(labels_coord[-4][0],labels_coord[-4][1],20,-2)

      plt.text(labels_coord[-5][0] + 20, labels_coord[-5][1]+5, labels_list[list(labels_list)[-5]])
      plt.arrow(labels_coord[-5][0],labels_coord[-5][1],20,+5)

      plt.text(labels_coord[-6][0] + 5, labels_coord[-6][1]+40, labels_list[list(labels_list)[-6]])
      plt.arrow(labels_coord[-6][0],labels_coord[-6][1],5,40)

      plt.text(labels_coord[0][0], labels_coord[0][1]+40, labels_list[list(labels_list)[-7]])
      plt.arrow(labels_coord[0][0],labels_coord[0][1],0,40)

      for l in labels_list:
        try:
          xl = x[names.index(l,1035)]
          yl = y[names.index(l,1035)]
          label2_coord.append([xl,yl])
        except:
          pass

      plt.arrow(label2_coord[-1][0],label2_coord[-1][1],labels_coord[-1][0] + 5 - label2_coord[-1][0], labels_coord[-1][1] + 5 - label2_coord[-1][1],ls=':')

      plt.arrow(label2_coord[-2][0],label2_coord[-2][1],labels_coord[-2][0] + 15 - label2_coord[-2][0],labels_coord[-2][1] + 5 - label2_coord[-2][1],ls=':')

      plt.arrow(label2_coord[-3][0],label2_coord[-3][1],labels_coord[-3][0] + 20 - label2_coord[-3][0],labels_coord[-3][1] + 5 - label2_coord[-3][1],ls=':')

      plt.arrow(label2_coord[-4][0],label2_coord[-4][1],labels_coord[-4][0] + 20 - label2_coord[-4][0], labels_coord[-4][1] - 2 - label2_coord[-4][1],ls=':')

      plt.arrow(label2_coord[-5][0],label2_coord[-5][1],labels_coord[-5][0] + 20 - label2_coord[-5][0], labels_coord[-5][1] + 5 - label2_coord[-5][1],ls=':')

      plt.arrow(label2_coord[-6][0],label2_coord[-6][1],labels_coord[-6][0] + 5 - label2_coord[-6][0], labels_coord[-6][1] + 40 - label2_coord[-6][1],ls=':')

      plt.arrow(label2_coord[0][0],label2_coord[0][1],labels_coord[0][0] - label2_coord[0][0] , labels_coord[0][1]+40-label2_coord[0][1],ls=':')

    if 'ZPR' in title:
      plt.text(labels_coord[-1][0] - 200, labels_coord[-1][1] -1700, labels_list[list(labels_list)[-1]])
      plt.arrow(labels_coord[-1][0],labels_coord[-1][1],-200,-1700)
    
      plt.text(labels_coord[-2][0] + 15, labels_coord[-2][1] -1200, labels_list[list(labels_list)[-2]])
      plt.arrow(labels_coord[-2][0],labels_coord[-2][1],15,-1200)

      plt.text(labels_coord[-3][0] -500, labels_coord[-3][1] -1500, labels_list[list(labels_list)[-3]])
      plt.arrow(labels_coord[-3][0],labels_coord[-3][1],-500,-1500)

      plt.text(labels_coord[-4][0] + 20, labels_coord[-4][1] -1000, labels_list[list(labels_list)[-4]])
      plt.arrow(labels_coord[-4][0],labels_coord[-4][1],20,-1000)

      plt.text(labels_coord[-5][0] + 400, labels_coord[-5][1] -1200, labels_list[list(labels_list)[-5]])
      plt.arrow(labels_coord[-5][0],labels_coord[-5][1],400,-1200)

      plt.text(labels_coord[-6][0] + 600, labels_coord[-6][1] -1000, labels_list[list(labels_list)[-6]])
      plt.arrow(labels_coord[-6][0],labels_coord[-6][1],600,-1000)

      plt.text(labels_coord[0][0]+600, labels_coord[0][1] -600, labels_list[list(labels_list)[-7]])
      plt.arrow(labels_coord[0][0],labels_coord[0][1],600,-600)

      for l in labels_list:
        try:
          xl = x[names.index(l,1035)]
          yl = y[names.index(l,1035)]
          label2_coord.append([xl,yl])
        except:
          pass

      plt.arrow(label2_coord[-1][0],label2_coord[-1][1],labels_coord[-1][0] -200 - label2_coord[-1][0], labels_coord[-1][1] - 1700 - label2_coord[-1][1],ls=':')

      plt.arrow(label2_coord[-2][0],label2_coord[-2][1],labels_coord[-2][0] + 15 - label2_coord[-2][0],labels_coord[-2][1] -1200 - label2_coord[-2][1],ls=':')

      plt.arrow(label2_coord[-3][0],label2_coord[-3][1],labels_coord[-3][0] -500 - label2_coord[-3][0],labels_coord[-3][1] -1500 - label2_coord[-3][1],ls=':')

      plt.arrow(label2_coord[-4][0],label2_coord[-4][1],labels_coord[-4][0] + 20 - label2_coord[-4][0], labels_coord[-4][1] -1000 - label2_coord[-4][1],ls=':')

      plt.arrow(label2_coord[-5][0],label2_coord[-5][1],labels_coord[-5][0] + 400 - label2_coord[-5][0], labels_coord[-5][1] -1200 - label2_coord[-5][1],ls=':')

      plt.arrow(label2_coord[-6][0],label2_coord[-6][1],labels_coord[-6][0] + 600 - label2_coord[-6][0], labels_coord[-6][1] -1000 - label2_coord[-6][1],ls=':')

      plt.arrow(label2_coord[0][0],label2_coord[0][1],labels_coord[0][0] +600 - label2_coord[0][0] , labels_coord[0][1] -600 - label2_coord[0][1],ls=':')

# if title:
#   plt.title(title)
   
  def update_annot(ind):
    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
#    annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
    annot.get_bbox_patch().set_alpha(0.4)


  def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
      cont, ind = sc.contains(event)
      if cont:
        update_annot(ind)
        annot.set_visible(True)
        fig.canvas.draw_idle()
      elif vis:
        annot.set_visible(False)
        fig.canvas.draw_idle()

  fig.canvas.mpl_connect("motion_notify_event", hover)
  if hardcopy:
    name = title.replace(' ','_')
    name = name.replace('$','')
    name = name.replace('{','')
    name = name.replace('}','')
    name = name.replace('\\','')
    name = name.replace('^','')
    name = name.replace('*','')
    if not title:
      name = 'ZPR_vs_alpha_standard_all_mats'

    plt.savefig(name+'.png', bbox_inches="tight")
  plt.show()


def gs_input(in_structure,psp_list,psp_dir,ecut=20,ngkpt=(4,4,4)):

#  from collections import Counter

  # Build AbinitInput from structure and pseudos from pseudo-dojo
  inp = abilab.AbinitInput(structure=in_structure, pseudos=psp_list,pseudo_dir=psp_dir)
  
  # Set value for other variables
  inp.set_vars(
    ecut=ecut,
    ngkpt=ngkpt,
    nshiftk=4,
    shiftk = [0.0, 0.0, 0.5,   # Shift for the FCC MK grid - the set_kmesh function should be 
             0.0, 0.5, 0.0,   # used here to set the shift in accordance to the cell type
             0.5, 0.0, 0.0,
             0.5, 0.5, 0.5],
    ixc=1,
    nstep=25,
    diemac=12.0,              # I know that this is a model dielectric function, but is there
                              # a better way to set it up?
    tolvrs=1.0e-10
  )
  # Define k-point sampling
#  inp.set_kmesh(ngkpt=(6, 6, 6), shiftk=(0, 0, 0))
  
  return inp

def gamma0_dfpt(in_structure,psp_list,psp_dir):
  scf_input = gs_input(in_structure,psp_list,psp_dir)
# print(scf_input)
# return flowtk.phonon_conv_flow("flow_alas_ecut_conv", scf_input, qpoints=(0, 0, 0),
#                                  params=["ecut", [20, 25, 30]])
  

# Energy cutoffs, kpoint sampling, and number of bands should be decided automatically 
# depending on the system. I have not found a way to obtain these from materialsproject
# but I guess there is a way to query it

# If there is none, k-points should be derived from symmetries(?) and ecut should be converged


###############################################################################


# Old stuff that I used at some point

#   from collections import Counter

# # Build structure from dictionary with input variables.
# a, b, c            = in_structure.lattice.abc
# alpha, beta, gamma = in_structure.lattice.angles

# t = []
# z = list(Counter(in_structure.species).values())
# for i in range(0,len(z)):
#   t += z[i]*[i+1]
# atom_types = t
# 
# coords = np.reshape(in_structure.frac_coords,(3*len(in_structure.frac_coords)),order='C')

# structure = abilab.Structure.from_abivars(
#      ntypat=in_structure.ntypesp,        # Number of atom types.
#      znucl=in_structure.atomic_numbers,  # Atomic numbers of the type(s) of atom.
#      natom=in_structure.num_sites,       # Number of atoms.
#      typat=atom_types,                   # List of atom types.
#      xred=coords,                        # Reduced coordinates
#      acell=[a,b,c],                      # Lengths of the primitive vectors (in Bohr).
#      angdeg=[alpha,beta,gamma]           # Angles between primitive vectors 
# )
# print('in # types of atoms',  in_structure.ntypesp,        'out # types of atoms',  structure.ntypesp)
# print('in # atomic numbers',  in_structure.atomic_numbers, 'out # atomic numbers',  structure.atomic_numbers)
# print('in # number of atoms', in_structure.num_sites,      'out # number of atoms', structure.num_sites)
#
#project = 'carrier_transport'
#client = load_client('zSHgsQzpW5sDel8WBhcX50OS2ziufHJM')
#
#print(dir(client))
#print(dir(client.projects))
#print(dir(client.redox_thermo_csp))
#
#r = client.projects.get_entry(pk=project, _fields=['title', 'authors', 'description', 'urls']).result()
#print(r)
#
#u = client.contributions.get_entries(project=project,identifier='mp-1138',_fields=['data']).result()
#
#print(u['data'][0]['data']['m\u2091\u1d9c\u1d52\u207f\u1d48|p']['\u03b5\u0304']['value'])
#
#
#identifier = 'mp-1002'
#cid = client.contributions.create_entry(contribution={
#    'project': project, 'identifier': identifier,
#    'data': {'E': '3.33 eV', 'E|V': {'a': 1, 'b': '3 cm'}}
#}).result()['id']
#client.contributions.get_entry(pk=cid, _fields=['id', 'identifier', 'data']).result()

#exit()
###############################################################################

#def __init__(self, mp_id, *kwargs):
#  if mp_id is None: 
#    raise ValueError('No structure given, cannot start input file.')
#  else:
#    self.mp_id = mp_id
#  self.initialize()

#
# Guido's database
#
#with open('phonon_dielectric_mp.json', 'r') as f:
#  guido_data = json.load(f)
#
#
# MPRester key needed for materials project - add it to .pmgrc.yaml to ease usage
#
#mpkey    = 'rMTANhpU06bjiedXX8Vg' 
#m = mg.MPRester(mpkey)

#m_e = 0.51099895000*10**(9)
#
#
#def check_epsilon(mpid):
# try:
#   with open('phonon/%s.json' %mpid, 'r') as f:
#     guido_data = json.load(f)
#     epsilon_t = guido_data['dielectric']['eps_total'][0][0]
#     epsilon_e = guido_data['dielectric']['eps_electronic'][0][0]
#     epsilon_s = 'guido'
#   return epsilon_t, epsilon_e, epsilon_s
# except:
#   try:
#     mpkey = 'rMTANhpU06bjiedXX8Vg' 
#     m     = mg.MPRester(mpkey)
#     data  = m.get_data(mp_id)
#     epsilon_t = data[0]['diel']['e_total']
#     epsilon_0 = data[0]['diel']['e_electronic']
#     epsilon_s = 'mp-data'
#     return epsilon_t, epsilon_e, epsilon_s
#   except: 
#     print("%s not in dielectric database" %mpid)
#     return 1.0, 1.0, False

""" ! """#!/usr/bin/env python3
#  -*- coding: utf-8 -*-
"""
Routine: 

Author: Pedro Melo
Email: pmmelo@uliege.be

Final Purpose: Check existing databases for parameters with which to compute alpha
               alpha = left(frac{1}{epsilon_infty} - frac{1}{epsilon_0}right)sqrt{frac{m*}{2omega_{LO}}}
               If no database or parameter exists, query MaterialsProject (MP) for the info. If info gathered from MP is incomplete, compute missing data with abinit
"""
import numpy as np

import argparse
import functions as fn
import anadata as an
import json

from os import path

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Check databases for variables needed to compute alpha. Check status of calculations for missing data. Updates master database')
parser.add_argument('-rf' ,'--ReadFile',  type=str, help="Read information from file",              default='')
parser.add_argument('-m'  ,'--MatProjID', type=str, help='User provided in-line MP-ID list',nargs='*', metavar="mp-ID",   default=[])
#parser.add_argument('-f'  ,'--Formula',   type=str, help='User provided in-line chemical formula', default='')
parser.add_argument('-gF' ,'--generalisedFroehlich', action='store_true', help='Run calculations for the Generalised Frohlich model')
parser.add_argument('-osF','--oldstandardFroelich',    action='store_true', help='(Outdated) Run calculations for the Standard Frohlich model highest phonon frequency')
parser.add_argument('-sF' ,'--standardFroelich',    action='store_true', help='Run calculations for the Standard Frohlich model')
parser.add_argument('-sfF','--standardFeynman',    action='store_true', help='Run calculations for the Standard Frohlich model based on Feynman path integral')
parser.add_argument('-ai' ,'--createAnaddbInp',     action='store_true', help='Create anaddb inputs to compute phonon eigendisplacements')
parser.add_argument('-ar' ,'--runAnaddb',           action='store_true', help='Run anaddb to compute phonon eigendisplacements')
parser.add_argument('-ce' ,'--compareEpsilon',      action='store_true', help='Compare dielectric function on Guido and MP database')
parser.add_argument('-fi' ,'--freqInterpol',   action='store_true', help='Compare interpolated and calculated phonon frequencies')
parser.add_argument('-di' ,'--dispInterpol',   action='store_true', help='Compare interpolated and calculated phonon eigendisplacements')
parser.add_argument('-to' ,'--testOrtho',       action='store_true', help='Test mode-polarity vector orthogonality')
parser.add_argument('-sd' ,'--showDisp',       action='store_true', help='Show eigendisplacements and q-point list')
parser.add_argument('-emdb','--effmassDB', type=str,  help=" Select the database from calculated effective masses: mp, abinit. If not explicit, it checks if exists by hierachy ( left to right in option)", default="all")
parser.add_argument('-gAEM' ,'--genAbinitEffMasses', action='store_true', help='Generate Input Files to calculate effective masses using abipy. Copy folder calc_effmasses and files inside extra_files/abinit_effmasses to the folder where effective masses runs will occur and start calculation using flow_effmass.py, abipy commands and in the end use join_files.py')
parser.add_argument('-c'  ,'--conduction', action='store_true', help='Conduction ZPR and alpha')
parser.add_argument('-v'  ,'--valence'   , action='store_true', help='Valence ZPR and alpha')
parser.add_argument('-dir' ,'--directory', type=str, help="Save the results in a target directory", default=".")

band = 1 # conduction as default
args = parser.parse_args()

argsdict = args.__dict__

opt = {}

type_list = ['valence','conduction']

def read_file(in_file):
  sys_list = np.loadtxt(in_file,dtype='str',ndmin=1, usecols=(0))
  if sys_list.size == 0:
      raise EmptyError("No material inside file")

  return sys_list

def run_gFr(opt,band):
  
#  f = open(target_dir + 'generalised_froelich_data','a')
#  f.write('#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \n' %('MPID', 'FORMULA', 'AVERAGE_EPSM1', 'AVERAGE_OMEGAM1', 'AVERAGE_MSTAR', 'AVERAGE_ALPHA', 'ZPR[meV]'))

  for mpid in sys_list:
    if path.exists("%s-data-per-mode-%s.dat" %(mpid,type_list[band])):
      print("Calculation for %s %s already finished" %(type_list[band],mpid))
    else:
      formula = fn.get_formula(mpid)
      average_epsm1, average_mstar, average_omega, average_alpha, zpr , zpr_s, check_sum, deg = fn.generalised_zpr(mpid, opt, band)
      #print(mpid, formula, average_epsm1, average_mstar, average_omega, average_alpha, zpr, zpr_s)
      if zpr_s:
        f = open(target_dir + '/gFr-%s-%s.dat' %(mpid,type_list[band]),'w')
        f.write('#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s\n' %('MPID', 'FORMULA', 'AVERAGE_EPSM1', 'AVERAGE_OMEGAM05', 'AVERAGE_MSTAR', 'AVERAGE_ALPHA', 'ZPR[meV]', 'AREA/4pi', 'DEGENERACIES'))
        f.write('%-20s \t %-20s \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f\n' %(mpid, formula, average_epsm1, average_omega, average_mstar, average_alpha, zpr, check_sum, deg))
        f.close()

def run_oldsFr(opt,band):
  f = open( target_dir + '/old_standard_frohlich_%s' %type_list[band],'w')
  #f.write('#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s\n' %('MPID', 'FORMULA', 'EPS_static', 'EPS_electronic', 'OMEGA_LO[meV]', 'MSTAR', 'ALPHA', 'ZPR[meV]'))
  f.write('#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s\n' %('MPID', 'FORMULA', 'EPS', 'OMEGA_LO[meV]', 'MSTAR', 'ALPHA', 'ZPR[meV]'))

  for mpid in sys_list:
    formula = fn.get_formula(mpid)
    #print(formula)
    #e_static, e_electronic, omega_lo, mstar, alpha, zpr , zpr_s = fn.standard_zpr(mpid, opt, band)
    eps, omega_lo, mstar, alpha, zpr , zpr_s = fn.standard_zpr(mpid, opt, band)
    #print(mpid, formula, e_static, e_electronic, omega_lo, mstar, alpha, zpr, zpr_s)
    if zpr_s:
      #f.write('%-20s \t %-20s \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f \t %-20s\n' %(mpid, formula, e_static, e_electronic, omega_lo, mstar, alpha, zpr))
      f.write('%-20s \t %-20s \t %-20f \t %-20f \t %-20f \t %-20f \t %-20s\n' %(mpid, formula, eps, omega_lo, mstar, alpha, zpr))
  f.close()

def run_sFr(opt,band):
  f = open(target_dir + '/standard_frohlich_data_%s' %type_list[band],'w')
  #f.write('#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s\n' %('MPID', 'FORMULA', 'EPS_static', 'EPS_electronic', 'OMEGA_LO[meV]', 'MSTAR', 'ALPHA', 'ZPR[meV]'))
  f.write('#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s\n' %('MPID', 'FORMULA', 'EPS', 'OMEGA_LO[meV]', 'MSTAR', 'ALPHA', 'ZPR[meV]'))

  atomic_info = fn.get_atomic_info()

  for mpid in sys_list:
    formula = fn.get_formula(mpid)
    #print(formula)
    #e_static, e_electronic, omega_lo, mstar, alpha, zpr , zpr_s, = fn.hellwarth_zpr(mpid, opt, band)
    atominfo = atomic_info[mpid]
    eps, omega_lo, mstar, alpha, zpr , zpr_s, = fn.hellwarth_zpr(mpid, opt, band, atominfo)
    #print(mpid, formula, e_static, e_electronic, omega_lo, mstar, alpha, zpr, zpr_s)
    if zpr_s:
      #f.write('%-20s \t %-20s \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f\n' %(mpid, formula, e_static, e_electronic, omega_lo, mstar, alpha, zpr))
      f.write('%-20s \t %-20s \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f\n' %(mpid, formula, eps, omega_lo, mstar, alpha, zpr))
  f.close()

def run_sfFr(opt,band):
  f = open( target_dir + '/standard_feynman_data_%s' %type_list[band],'w')
  #f.write('#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s\n' %('MPID', 'FORMULA', 'EPS_total', 'EPS_electronic', 'OMEGA_LO[meV]', 'MSTAR', 'ALPHA', 'ZPR[meV]'))
  f.write('#%-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s \t %-20s\n' %('MPID', 'FORMULA', 'EPS', 'OMEGA_LO[meV]', 'MSTAR', 'ALPHA', 'ZPR[meV]'))

  atomic_info = fn.get_atomic_info()
  for mpid in sys_list:
    formula = fn.get_formula(mpid)
    #print(formula)
    #e_total, e_electronic, omega_lo, mstar, alpha, zpr , zpr_s = fn.feynman_zpr(mpid, opt, band)
    atominfo = atomic_info[mpid]
    eps, omega_lo, mstar, alpha, zpr , zpr_s = fn.feynman_zpr(mpid, opt, band, atominfo)
    #print(mpid, formula, e_total, e_electronic, omega_lo, mstar, alpha, zpr, zpr_s)
    if zpr_s:
      #f.write('%-20s \t %-20s \t %-20f \t %-20f \t %-20f \t %-20f \t %-20f \t %-20s\n' %(mpid, formula, e_total, e_electronic, omega_lo, mstar, alpha, zpr))
      f.write('%-20s \t %-20s \t %-20f \t %-20f \t %-20f \t %-20f \t %-20s\n' %(mpid, formula, eps, omega_lo, mstar, alpha, zpr))
  f.close()



def anaddb_inp():
  for mpid in sys_list:
    an.generate_anaddb_input(mpid)

def anaddb_run():
  for mpid in sys_list:
    an.run_anaddb(mpid)

def comp_eps():
  for mpid in sys_list:
    fn.compare_epsilon(mpid)

def run_gAbinitEM():
    ab.genEffMasses(sys_list)

def frequency_interpolation():
  for mpid in sys_list:
    coeff_mat, error = an.read_omega_q(mpid)
    print(coeff_mat, error)

def displacements_interpolation():
  for mpid in sys_list:
    coeff_mat, coeff_mat_img, error, error_img, modes, nqpts, atoms = an.read_eigendisplacements(mpid)
    print(coeff_mat, coeff_mat_img, error, error_img)

def test_orto():
  for mpid in sys_list:
    an.test_orthogonality(mpid)

def show_disps():
  for mpid in sys_list:
    an.get_q_dot_p_squared(mpid)

# First check options

flag_in = False
if argsdict["ReadFile"] != '':
  in_file = str(args.ReadFile)
  print( "Using %s file"%in_file )
  sys_list = read_file(in_file)
  flag_in = True

if argsdict["MatProjID"] != []:
  if flag_in:
      raise Exception("Only one input option permitted: -rf or -m")
  sys_list= argsdict["MatProjID"]
  for il in sys_list:
    if "mp-" not in il:
        raise Exception("Missing mp- initials")
  flag_in = True

if not flag_in:
    raise Exception("Input option necessary: -rf or -m")

#if argsdisct["Formula"] != []:
#  if flag_in:
#      raise Exception("Only one input option permitted: -rf, -m or -f")
#  formula_list = argsdisct["Formula"]



if args.conduction:
  band  = 1

if args.valence:
  band  = 0

if args.directory:
    target_dir = args.directory

if args.effmassDB:
  opt["eff_masses"] = args.effmassDB

# Then check methods

method = 0

if args.genAbinitEffMasses:
  run_gAbinitEM()
  method = True

if args.oldstandardFroelich:
  run_oldsFr(opt,band)
  method = True

if args.standardFroelich:
  run_sFr(opt,band)
  method = True

if args.standardFeynman:
  run_sfFr(opt,band)
  method = True

if args.generalisedFroehlich:
  run_gFr(opt,band)
  method = True

if args.createAnaddbInp:
  anaddb_inp()
  method = True

if args.runAnaddb:
  anaddb_run()
  method = True

if args.compareEpsilon:
  comp_eps()
  method = True

if args.freqInterpol:
  frequency_interpolation()
  method = True

if args.dispInterpol:
  displacements_interpolation()
  method = True

if args.testOrtho:
  test_orto()
  method = True

if args.showDisp:
  show_disps()
  method = True

if not method:
    print("Only options for materials input were set. Other options should be added to work")
  

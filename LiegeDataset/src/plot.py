#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

import numpy as np
import argparse
import matplotlib.pyplot as plt
from math import ceil


#parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Plot data from ALPS output files')

#parser.add_argument('-b' ,'--Band',  type=str, nargs="?",
#        #choices=["v","valence","c","conduction","vc","valcond"],
#        help="BAND: 'v' or 'valence', 'c' or 'conduction',"+
#        " 'vc' or 'valcond' for both. Ex: -b c. Default=vc", default='vc')
parser.add_argument('-y' ,'--Yaxis',  type=str,  choices=
        ["ZPR","alpha","mstar","omega","eps"],
        help="Default=ZPR."+
        "Where ZPR=Zero-Point Renormalization, alpha=coupling strength, "+
        "mstar=effective mass, omega=phonon frequency, "+
        "eps= effective dielectric constant(Generalized Frohlich only)", default='ZPR')
parser.add_argument('-x' ,'--Xaxis',  type=str, choices=
        ["ZPR","alpha","mstar","omega","eps"],
        help="Default=alpha."+
        "Where ZPR=Zero-Point Renormalization, alpha=coupling strength, "+
        "mstar=effective mass, omega=phonon frequency, "+
        "eps= effective dielectric constant", default='alpha')

parser.add_argument('-sF','--StandardFrohlich',nargs='*', metavar="FILE", 
        help="A list of output files from the standard Frohlich calculations. The file should contain valence or conduction in the name to identify the band."+#
        " Ex: -sF FILE1_valence FILE2valence FILE3-conduction FILE4.conduction", default=[])

parser.add_argument('-gF','--GeneralizedFrohlich',nargs='*', metavar="FILE", 
        help="A list of output files from the generalized Frohlich calculations. The file should contain valence or conduction in the name to identify the band."+#
        " Ex: -gF FILE1_valence FILE2valence FILE3-conduction FILE4.conduction", default=[])


parser.add_argument('-cx','--ColumnXaxis', type=int, nargs='?', metavar="Column number",
        help="Column of the data for the x-axis.")

parser.add_argument('-cy','--ColumnYaxis', type=int, nargs='?', metavar="Column number",
        help="Column of the data for the y-axis.")


args = parser.parse_args()


yaxis = "ZPR"
xaxis = "alpha"

if args.Xaxis:
    xaxis = args.Xaxis

if args.Yaxis:
    yaxis = args.Yaxis



map_axis = {
        "eps" :{
            "keyword":[["EPS"],["AVERAGE_EPSM1"]],
            "column": 3,
            "name" : [r"$\varepsilon^{\infty}$",r"$\left\langle \frac{1}{\varepsilon^*} \right\rangle$"]
            },
        "ZPR" : {
            "keyword":[["ZPR"],["ZPR"]],
            "column": 7,
            "name" : [r"ZPR [meV]","ZPR [meV]"]
            },
        "alpha" : {
            "keyword":[["ALPHA"],["AVERAGE_ALPHA"]],
            "column": 6,
            "name" : [r"$\alpha$",r"$\left\langle \alpha \right\rangle$"]
            },
        "mstar" : {
            "keyword":[["MSTAR"],["AVERAGE_MSTAR"]],
            "column": 5,
            "name" : [r"$\left\langle m^* \right\rangle$"]
            },
        "omega" : {
            "keyword":[["OMEGA_LO"],["AVERAGE_OMEGAM05"]],
            "column": 4,
            "name" : [r"$\omega_{\text{LO}}$ [meV]",r"$\left\langle \frac{1}{\sqrt{\omega_\text{LO}}} \right\rangle$"]
            }
        }

if args.ColumnXaxis:
    map_axis[xaxis]["column"] = args.ColumnXaxis

if args.ColumnYaxis:
    map_axis[yaxis]["column"] = args.ColumnYaxis


#type_list = ["valence","conduction"]


#if args.Band:
#    if args.Band =="v" or args.Band =="valence":
#        type_list=["valence"]
#    elif args.Band =="c" or args.Band =="conduction":
#        type_list=["conduction"]
#    elif args.Band =="vc" or args.Band =="valcond":
#        type_list = ["valence","conduction"]


def plot(x,y,xlabel,ylabel,flag):
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.scatter(x,y)
    if flag:
        ax.plot((-max(x),max(x)),(0,0),"--",color="black")
        ax.plot((0,0),(-max(y),max(y)),"--",color="black")
        ax.set_xlim(-ceil(max(x)),ceil(max(x)))
        #ax.set_ylim(-ceil(max(y)/10)*10,ceil(max(y)/10)*10)
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels([abs(i) for i in ax.get_xticks()])
    plt.show()


def run(xaxis,yaxis,flag,files,option=0):
    x_values = []
    y_values = []
    for f in files:
        present_x = False
        present_y = False
        i_line = -1
        with open(f,"r") as fh:
            for lines in fh:
                i_line+=1
                l_split = lines.split()
                if l_split == []:
                    print("Empty line or end of file")
                    break
                if "#" == lines[0] and i_line == 0:
                    for ma in range(len(map_axis[xaxis]["keyword"][option])):
                        imap_axis = map_axis[xaxis]["keyword"][option][ma]
                        for il in range(len(l_split)):
                            isplit = l_split[il]
                            if imap_axis in isplit:
                                present_x = True
                                column_x = il+1
                    if not present_x:
                        print("Chosen X-axis '"+
                                xaxis+"' not present in the comments line by any of the following keywords:\n",
                                map_axis[xaxis]["keyword"][option],"\nContinuing by reading column "+
                                str( map_axis[xaxis]["column"]) + ". To change use option -cx .")
                    else:
                        if column_x != map_axis[xaxis]["column"]:
                            print("Keyword for x-axis located at column "+
                                    str(column_x)+". Default column for "+
                                    xaxis+" is "+str(map_axis[xaxis]["column"])+". Using default column "+
                                    str(map_axis[xaxis]["column"])+" instead.")
                            #map_axis[xaxis]["column"] = column_x
                    for ma in range(len(map_axis[yaxis]["keyword"][option])):
                        imap_axis = map_axis[yaxis]["keyword"][option][ma]
                        for il in range(len(l_split)):
                            isplit = l_split[il]
                            if imap_axis in isplit:
                                present_y = True
                                column_y = il+1
                    if not present_y:

                               print("Chosen Y-axis '"+
                                       yaxis+"' not present in the comments line by any of the following keywords:\n",
                                       map_axis[yaxis]["keyword"][option],"\nContinuing by reading column "+
                                       str( map_axis[yaxis]["column"]) + ". To change use option -cy .")
                    else:
                        if column_y != map_axis[yaxis]["column"]:

                            print("Keyword for y-axis located at column "+
                                    str(column_y)+". Default column for "+
                                    yaxis+" is "+str(map_axis[yaxis]["column"])+". Using default column "+
                                    str(map_axis[yaxis]["column"])+" instead.")

                            #map_axis[yaxis]["column"] = column_y
                else:
                    xval = float(l_split[map_axis[xaxis]["column"]-1])
                    yval = float(l_split[map_axis[yaxis]["column"]-1])
                    if flag:
                        if "valence" in f:
                            yval = -1*yval
                        if "conduction" in f:
                            xval = -1*xval
                    x_values.append(xval)
                    y_values.append(yval)




    plot(x_values,y_values,map_axis[xaxis]["name"][option], map_axis[yaxis]["name"][option],flag)
                   

if args.StandardFrohlich:
    if args.GeneralizedFrohlich:
        raise Exception("It can only be -sF or -gF")
    files = args.StandardFrohlich
    flag_cond = False
    flag_val = False
    flag = False
    for i in files:
        if "conduction" in i:
            flag_cond = True
        if "valence" in i:
            flag_val = True

    if flag_cond and flag_val:
        flag = True
    run(xaxis,yaxis,flag,files,0)

if args.GeneralizedFrohlich:
    if args.StandardFrohlich:
        raise Exception("It can only be -sF or -gF")
    files = args.GeneralizedFrohlich
    flag_cond = False
    flag_val = False
    flag = False
    for i in files:
        if "conduction" in i:
            flag_cond = True
        if "valence" in i:
            flag_val = True

    if flag_cond and flag_val:
        flag = True

    run(xaxis,yaxis,flag,files,1)
    


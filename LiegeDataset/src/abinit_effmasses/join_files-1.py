import sys
import os
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.flowtk as flowtk
import ast
import numpy as np
import time
import tarfile

pseudo_dir = "/scratch/users/j/a/jabreu/pseudos/PBEsol/"


workdir = "flow"

def is_float(string):
  try:
      float(string) and '.' in string  # True if string is a number contains a dot
      return True
  except ValueError:  # String is not a number
    return False
def is_int(string):
  try:
    int(string) and not "." in string  # True if string is a number contains a dot
    return True
  except ValueError:  # String is not a number
    return False

def startorclean_inp():
    dict_inp = {
        "number" : -1,
        "connect_inp" : -1,
        "path" : "",
        "mpid" : "",
        "info" : {
        "pseudos" : [],
        "pseudos_znucl" : [],
        "natom" : {"type": "i", "data": []},
        "ntypat" : {"type": "i", "data": []},
        "typat" : {"type": "ai", "data": []},
        "znucl" : {"type": "af", "data": []},
        "xred" : {"type": "af", "data": []},
        "acell" : {"type": "af", "data": []},
        "rprim" : {"type": "af", "data": []},
        }}
    readdata = []
    return dict_inp, readdata


def check_finish(path):

    with open(path,"r") as f:
        for lines in f:
            if "Calculation completed" in lines:
                return True
        return False

def read_inp(path):
    path_split = path.split("/")
    dict_inp, readdata = startorclean_inp()
    flag = False
    dict_inp["path"] = path
    dict_inp["number"] = int(path_split[1][1:])
    dict_inp["done"] = False
    if not os.path.exists(path_split[0] + "/" + path_split[1] + "/t1/run.abo"):
        return dict_inp
    if not os.path.exists(path_split[0] + "/" + path_split[1] + "/t0/run.abi"):
        return dict_inp
    with open(path_split[0] + "/" + path_split[1] + "/t1/run.abo", "r") as f:
        for lines in f:
            if "Calculation completed" in lines:
                dict_inp["done"] = True
    with open(path_split[0] + "/" + path_split[1] + "/t0/run.abi", "r") as f:
        for lines in f:
            if "STRUCTURE" in lines:
                flag = True
            if not flag:
                continue
            readdata.append(lines.split())
    idata = 0
    while idata < len(readdata):
        l_split = readdata[idata]
        while l_split == []:
            idata+=1
            l_split = readdata[idata]

        for item in dict_inp["info"]:
            if item == l_split[0]:
                if dict_inp["info"][item]["type"] == "i":
                    dict_inp["info"][item]["data"].append(int(l_split[1]))
                if "a" in dict_inp["info"][item]["type"]:
                    flag2 = False
                    if len(l_split) > 1:
                        for il in l_split[1:]:
                            if "i" in  dict_inp["info"][item]["type"]:
                                dict_inp["info"][item]["data"].append(int(il))
                            elif "f" in  dict_inp["info"][item]["type"]:
                                dict_inp["info"][item]["data"].append(float(il))
                    while is_float(readdata[idata+1][0]) or is_int(readdata[idata+1][0]):
                        flag2 = True
                        idata+=1
                        l_split = readdata[idata]
                        for il in l_split:
                            if "i" in  dict_inp["info"][item]["type"]:
                                dict_inp["info"][item]["data"].append(int(il))
                            elif "f" in  dict_inp["info"][item]["type"]:
                                il=il.replace("d","e")
                                dict_inp["info"][item]["data"].append(float(il))
                        if readdata[idata+1] == []:
                            break


        idata+=1

    with open(path_split[0] + "/" + path_split[1] + "/t0/run.files", "r") as f:
        for lines in f:
            if "psp8" in lines:
                l_split = lines.split("/")
                pseudo = l_split[-1].rstrip()
                dict_inp["info"]["pseudos"].append(pseudo)
                g = open(pseudo_dir + pseudo,"r")
                for pos, line_pseudo in enumerate(g):
                    if pos ==1 :
                        dict_inp["info"]["pseudos_znucl"].append(float(line_pseudo.split()[0]))
                    elif pos > 2:
                        break


    return dict_inp



efmas_inp = []
gsr_inp = []
efmas_tar = tarfile.open("EFMAS.tar.gz", "w:gz")
gsr_tar = tarfile.open("GSR.tar.gz", "w:gz")
output_tar = tarfile.open("OUTPUT.tar.gz", "w:gz")

link_list = []

with open("link_calc_id", "r") as f:
    for lines in f:
        l_split = lines.split()
        if l_split == []:
            continue
        link_list.append([int(l_split[0]),l_split[1]])


link_calc_file = open("link_calc_id_end", "w")

for path, subdirs, files in os.walk(workdir):
        if "out_EFMAS.nc" in files:
            path_file = path + "/out_EFMAS.nc"
            dict_inp = read_inp(path)
            efmas_inp.append(dict_inp)
            #print(dict_inp)
for path, subdirs, files in os.walk(workdir):
    def check_files():
        if "out_GSR.nc" in files:
            path_file = path + "/out_GSR.nc"
            #print(path)
            path_split = path.split("/")
            if path_split[2] == "t1":
                    return
            number = int(path_split[1][1:])
            for ei in efmas_inp:
                    if number == ei["number"]:
                        return

            dict_inp = read_inp(path)
            for ll in link_list:
                    if ll[0] == number:
                        dict_inp["mpid"] = ll[1]
                        break

            for ei in efmas_inp:
                    if dict_inp["info"] == ei["info"] and ei["done"]:
                        dict_inp["connect_inp"] = ei["number"]
                        ei["connect_inp"] = dict_inp["number"]
                        ei["mpid"] = dict_inp["mpid"]
                        efmas_tar.add(ei["path"] + "/out_EFMAS.nc","%s_EFMAS.nc"%ei["mpid"])
                        output_tar.add(ei["path"]+"/../run.abo","%s_output"%ei["mpid"])
                        gsr_tar.add(path_file,"%s_GSR.nc"%dict_inp["mpid"])
                        link_calc_file.write("%d %s %d\n"%(dict_inp["number"], dict_inp["mpid"], ei["number"]))
                        break



    check_files()

efmas_tar.close()
output_tar.close()
gsr_tar.close()
link_calc_file.close()


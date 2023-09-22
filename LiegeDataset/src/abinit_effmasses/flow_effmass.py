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


def read_inp(path):
    path_split = path.split("/")
    dict_inp, readdata = startorclean_inp()
    flag = False
    dict_inp["path"] = path
    dict_inp["number"] = int(path_split[1][1:])
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





def read_input_vars(filename):
    input_vars = {}
    with open(filename) as f:
        for lines in f:
            l1_split = lines.split()
            l_split = lines.split(l1_split[0])
            if "[" in lines:
                input_vars[l1_split[0]] = ast.literal_eval(l_split[1].rstrip().lstrip())
            else:
                input_vars[l1_split[0]] = int(float(l_split[1].rstrip().lstrip()))
    input_vars["nshiftk"] = len(input_vars["shiftk"])
    if input_vars["nshiftk"] == 1:
        input_vars["shiftk"] = input_vars["shiftk"][0]
    for ps in range(len(input_vars["pseudos"])):
        input_vars["pseudos"][ps] = pseudo_dir + input_vars["pseudos"][ps]
    return input_vars


def make_scf_input(structure_cif, input_vars):
    """Returns input for GS-SCF calculation."""

    structure = abilab.Structure.from_file(structure_cif)
    scf_input = abilab.AbinitInput(structure=structure, pseudos=input_vars["pseudos"])

    # Global variables
    scf_input.set_vars(
        ecut=input_vars["ecut"],
        nband=input_vars["nband"]+9,
        paral_kgb=0,
        autoparal=1,
        nbdbuf=4,
        nstep=2000,
        tolvrs=1e-12,
        istwfk="*1",
        ngkpt=input_vars["ngkpt"],
        nshiftk=input_vars["nshiftk"],
        shiftk=input_vars["shiftk"],
    )

    return scf_input

def build_work(scf_input):
    from abipy.flowtk.effmass_works import EffMassDFPTWork, EffMassAutoDFPTWork
    work = EffMassAutoDFPTWork.from_scf_input(scf_input, ndivsm=-30, tolwfr=1e-10)
    return work


if __name__=="__main__":
    data = []
    name = []
    print("hello")

    manager = abilab.TaskManager.from_user_config()
    workdir = "flow"
    flow = flowtk.Flow(workdir=workdir,manager= manager)

    g = open("link_calc_id","w")

    root = os.getcwd()

    #flows = []
    for list_calc in os.listdir("calc_effmasses"):
        print(list_calc)
        lc = list_calc.split(".")
        name_calc= lc[0]
        type_calc = lc[1]

        if not name_calc in name:
            name.append(name_calc)
            data.append([None,None])
        idx_calc = name.index(name_calc)
        if list_calc in data[idx_calc]:
            continue
        else:
            if type_calc == "cif":
                data[idx_calc][0] = list_calc
            elif type_calc == "inp":
                data[idx_calc][1] = list_calc
    for i, n in enumerate(name):
        print(n)
        input_vars = read_input_vars("calc_effmasses/" +data[i][1])
        scf_input = make_scf_input("calc_effmasses/"+data[i][0],input_vars)
        work =  build_work(scf_input)
        flow.register_work(work)
        g.write("%s %s \n"%(i,n))
    g.close()
    flow.get_graphviz()
    flow.show_status()
    if os.path.isdir(workdir):
        import shutil
        shutil.rmtree(workdir)
        
    flow.build_and_pickle_dump()


    print("flow completed")


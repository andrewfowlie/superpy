#! /usr/bin/env python

"""
"""

__author__ = "Kazuki Sakurai <kazuki.sakurai@desy.de>"
__version__ = "0.1"

import os.path
import itertools
from math import *
from basic_func import *
from pyslha import *
import sys


def get_SRs_for_model(procs):
    sr_dict = {}
    accum = 0
    for prc, pData in sorted(procs.items(), key=lambda x:x[1].rate, reverse=True):
        accum += pData.rate                  
        if accum > 0.99: break 
        #print prc, accum
        sr_dict_proc = get_SRs_for_proc(procs, prc)
        for SR, srData in sr_dict_proc.items():
            if not SR in sr_dict:
                Data = Structure()
                Data.xsec = srData.xsec
                Data.MET = srData.xsec * srData.MET 
                Data.proc_dict = {}
                Data.proc_dict[prc] = srData
                sr_dict[SR] = Data
            else:
                sr_dict[SR].xsec += srData.xsec
                sr_dict[SR].MET += srData.xsec * srData.MET
                sr_dict[SR].proc_dict[prc] = srData
    for SR in sr_dict: sr_dict[SR].MET = sr_dict[SR].MET / sr_dict[SR].xsec
    return sr_dict



def get_SRs_for_proc(procs, prc):
    fp_dict = get_fp_dict(procs, prc)
    #print prc
    #print procs[prc].xsec    
    #for key, val in fp_dict.items():
    #    print key, val
    #print " "
    ######  cut #######
    # jet
    pjcut = 20
    try:
        for i in range(len(fp_dict["jet"])): 
            if fp_dict["jet"] < pjcut: fp_dict["jet"].pop([i])
    except: pass
    # b-jet
    pbcut = 20
    try:
        for i in range(len(fp_dict["bjet"])): 
            if fp_dict["bjet"] < pbcut: fp_dict["bjet"].pop([i])
    except: pass
    ################################
    add = []
    ######  z -> (j,j), (l+,l-), (nu, nu) #####
    mz = 91.18
    br_lep = 0.3366
    br_tau = 0.3366
    br_inv = 0.200
    br_had = 1. - br_lep - br_tau - br_inv  
    zjet_beta = 0.9
    try:    
        for p in fp_dict["z"]:
            mode = []
            e = sqrt(p**2 + mz**2)
            v = p / e        
            if v > zjet_beta: obj_dict = { "zjet" : 1 }             
            else: obj_dict = {"jet" : 2}
            mode.append([br_had, obj_dict])        
            obj_dict = { "l+" : 1, "l-" : 1}
            mode.append([br_lep, obj_dict])
            obj_dict = { "tau+" : 1, "tau-" : 1 }
            mode.append([br_tau, obj_dict])
            obj_dict = { "MET" : p }
            mode.append([br_inv, obj_dict])
            add.append(mode)
    except: pass
    ######  h0 -> (b,b), (tau+,tau-), #####
    mh = 125.5
    br_b = 0.569
    br_tau = 0.0624
    br_had = 1 - br_b - br_tau
    higgsjet_beta = 0.9
    try:    
        for p in fp_dict["h0"]:
            mode = []
            e = sqrt(p**2 + mh**2)
            v = p / e
            if v > higgsjet_beta: 
                obj_dict = { "higgsjet" : 1 }             
            else: 
                obj_dict = {"bjet" : 2}
            mode.append([br_b, obj_dict])        
            obj_dict = { "tau+" : 1, "tau-" : 1}
            mode.append([br_tau, obj_dict])
            obj_dict = { "jet" : 4}
            mode.append([br_had, obj_dict])
            add.append(mode)
    except: pass
    ######  w -> (j,j), (l,nu) #####
    mw = 80.4
    br_lep = 0.1080
    br_tau = 0.1125
    br_had = 1. - br_lep - br_tau  
    wjet_beta = 0.9    
    try:
        for p in fp_dict["w"]:
            mode = []
            e = sqrt(p**2 + mw**2)
            v = p / e        
            if v > wjet_beta: obj_dict = { "wjet" : 1 }             
            else: obj_dict = {"jet" : 2 }
            mode.append([br_had, obj_dict])        
            obj_dict = { "l+" : 1, "MET" : mw/2 }
            mode.append([br_lep/2, obj_dict])
            obj_dict = { "tau+" : 1, "MET" : mw/2 }
            mode.append([br_tau/2, obj_dict])
            obj_dict = { "l-" : 1, "MET" : mw/2 }
            mode.append([br_lep/2, obj_dict])
            obj_dict = { "tau-" : 1, "MET" : mw/2 }
            mode.append([br_tau/2, obj_dict])
            add.append(mode)
    except: pass
    try:
        for p in fp_dict["w+"]:
            mode = []
            e = sqrt(p**2 + mw**2)
            v = p / e        
            if v > wjet_beta: obj_dict = { "wjet" : 1 }             
            else: obj_dict = {"jet" : 2 }
            mode.append([br_had, obj_dict])        
            obj_dict = { "l+" : 1, "MET" : mw/2 }
            mode.append([br_lep, obj_dict])
            obj_dict = { "tau+" : 1, "MET" : mw/2 }
            mode.append([br_tau, obj_dict])
            add.append(mode)
    except: pass
    try:
        for p in fp_dict["w-"]:
            mode = []
            e = sqrt(p**2 + mw**2)
            v = p / e        
            if v > wjet_beta: obj_dict = { "wjet" : 1 }             
            else: obj_dict = {"jet" : 2 }
            mode.append([br_had, obj_dict])        
            obj_dict = { "l-" : 1, "MET" : mw/2 }
            mode.append([br_lep, obj_dict])
            obj_dict = { "tau-" : 1, "MET" : mw/2 }
            mode.append([br_tau, obj_dict])
            add.append(mode)    
    except: pass
    ######  top -> b + [(j,j), (l,nu)] #####
    mtop = 172.5
    br_lep = 0.1080
    br_tau = 0.1125
    br_had = 1. - br_lep - br_tau  
    topjet_beta = 0.95    
    try:
        for p in fp_dict["top"]:
            mode = []
            e = sqrt(p**2 + mtop**2)
            v = p / e        
            if v > topjet_beta: obj_dict = { "topjet" : 1 }             
            else: obj_dict = {"bjet" : 1, "jet" : 2}
            mode.append([br_had, obj_dict])        
            obj_dict = { "bjet" : 1, "l+" : 1, "MET" : mw/2 }
            mode.append([br_lep/2, obj_dict])
            obj_dict = { "bjet" : 1, "tau+" : 1, "MET" : mw/2 }
            mode.append([br_tau/2, obj_dict])
            obj_dict = { "bjet" : 1, "l-" : 1, "MET" : mw/2 }
            mode.append([br_lep/2, obj_dict])
            obj_dict = { "bjet" : 1, "tau-" : 1, "MET" : mw/2 }
            mode.append([br_tau/2, obj_dict])
            add.append(mode)
    except: pass
    try:
        for p in fp_dict["top+"]:
            mode = []
            e = sqrt(p**2 + mtop**2)
            v = p / e        
            if v > topjet_beta: obj_dict = { "topjet" : 1 }             
            else: obj_dict = {"bjet" : 1, "jet" : 2}
            mode.append([br_had, obj_dict])        
            obj_dict = { "bjet" : 1, "l+" : 1, "MET" : mw/2 }
            mode.append([br_lep, obj_dict])
            obj_dict = { "bjet" : 1, "tau+" : 1, "MET" : mw/2 }
            mode.append([br_tau, obj_dict])
            add.append(mode)
    except: pass
    try:
        for p in fp_dict["top-"]:
            mode = []
            e = sqrt(p**2 + mtop**2)
            v = p / e        
            if v > topjet_beta: obj_dict = { "topjet" : 1 }             
            else: obj_dict = {"bjet" : 1, "jet" : 2}
            mode.append([br_had, obj_dict])        
            obj_dict = { "bjet" : 1, "l-" : 1, "MET" : mw/2 }
            mode.append([br_lep, obj_dict])
            obj_dict = { "bjet" : 1, "tau-" : 1, "MET" : mw/2 }
            mode.append([br_tau, obj_dict])
            add.append(mode)
    except: pass

    final_list = []
    if len(add) == 0:
        xsec = procs[prc].xsec
        final = {}
        try: final["jet"] = len(fp_dict["jet"])  
        except: pass
        try: final["bjet"] = len(fp_dict["bjet"])  
        except: pass
        try: final["l+"] = len(fp_dict["l+"])  
        except: pass
        try: final["l-"] = len(fp_dict["l-"])  
        except: pass
        try: final["tau+"] = len(fp_dict["tau+"])  
        except: pass
        try: final["tau-"] = len(fp_dict["tau-"])  
        except: pass
        try: final["gam"] = len(fp_dict["gam"])  
        except: pass        
        try: final["MET"] = fp_dict["MET"]  
        except: pass
        final_list.append([xsec, final])
    else:
        #print "# of decaying particles =", len(add)
        #print add
        if len(add) == 1: 
            add.append([[1, {"dmmy": 1}]])
            comb = list(reduce(itertools.product, add))    
        else:
            comb = list(reduce(itertools.product, add))
        #print "# of combinations", len(comb)
        #print "comb", comb
        for nested_elements in comb:        
            elements = flatten_tuple(nested_elements)
            #print elements
            xsec = procs[prc].xsec
            final = {}
            try: final["jet"] = len(fp_dict["jet"])  
            except: pass
            try: final["bjet"] = len(fp_dict["bjet"])  
            except: pass
            try: final["l+"] = len(fp_dict["l+"])  
            except: pass
            try: final["l-"] = len(fp_dict["l-"])  
            except: pass
            try: final["tau+"] = len(fp_dict["tau+"])  
            except: pass
            try: final["tau-"] = len(fp_dict["tau-"])  
            except: pass
            try: final["MET"] = fp_dict["MET"]  
            except: pass
            for mode in elements:
                #print "mode", mode
                br = mode[0]
                xsec = xsec * br
                obj_dict = mode[1]
                if "jet" in obj_dict:
                    try: final["jet"] += obj_dict["jet"]
                    except: final["jet"] = obj_dict["jet"]
                if "bjet" in obj_dict:
                    try: final["bjet"] += obj_dict["bjet"]
                    except: final["bjet"] = obj_dict["bjet"]
                if "wjet" in obj_dict:
                    try: final["wjet"] += obj_dict["wjet"]
                    except: final["wjet"] = obj_dict["wjet"]
                if "zjet" in obj_dict:
                    try: final["zjet"] += obj_dict["zjet"]
                    except: final["zjet"] = obj_dict["zjet"]
                if "topjet" in obj_dict:
                    try: final["topjet"] += obj_dict["topjet"]
                    except: final["topjet"] = obj_dict["topjet"]
                if "higgsjet" in obj_dict:
                    try: final["higgsjet"] += obj_dict["higgsjet"]
                    except: final["higgsjet"] = obj_dict["higgsjet"]
                if "l+" in obj_dict:
                    try: final["l+"] += obj_dict["l+"]
                    except: final["l+"] = obj_dict["l+"]
                if "l-" in obj_dict:
                    try: final["l-"] += obj_dict["l-"]
                    except: final["l-"] = obj_dict["l-"]
                if "tau+" in obj_dict:
                    try: final["tau+"] += obj_dict["tau+"]
                    except: final["tau+"] = obj_dict["tau+"]
                if "tau-" in obj_dict:
                    try: final["tau-"] += obj_dict["tau-"]
                    except: final["tau-"] = obj_dict["tau-"]
                if "MET" in obj_dict:
                    try: final["MET"] = get_quad_ave(final["MET"], obj_dict["MET"])                
                    except: final["MET"] = obj_dict["MET"]                
            final_list.append([xsec, final])

    for i in range(len(final_list)):
        try:
            if final_list[i][1]["l+"] == 1 and final_list[i][1]["l-"] == 1:
                del final_list[i][1]["l+"]
                del final_list[i][1]["l-"]
                final_list[i][1]["OSlep"] = ""
        except: pass
        try:
            if final_list[i][1]["l+"] == 2 and not "l-" in final_list[i][1]:
                del final_list[i][1]["l+"]
                final_list[i][1]["SSlep"] = ""
        except: pass
        try:
            if final_list[i][1]["l-"] == 2 and not "l+" in final_list[i][1]:
                del final_list[i][1]["l-"]
                final_list[i][1]["SSlep"] = ""
        except: pass
        try:
            if final_list[i][1]["tau+"] == 1 and final_list[i][1]["tau-"] == 1:
                del final_list[i][1]["tau+"]
                del final_list[i][1]["tau-"]
                final_list[i][1]["OStau"] = ""
        except: pass
        try:
            if final_list[i][1]["tau+"] == 2 and not "tau-" in final_list[i][1]:
                del final_list[i][1]["tau+"]
                final_list[i][1]["SStau"] = ""
        except: pass
        try:
            if final_list[i][1]["tau-"] == 2 and not "tau+" in final_list[i][1]:
                del final_list[i][1]["tau-"]
                final_list[i][1]["SStau"] = ""
        except: pass

    for i in range(len(final_list)):
        nlep = 0
        try:
            nlep += final_list[i][1]["l+"]
            del final_list[i][1]["l+"]
        except: pass
        try:
            nlep += final_list[i][1]["l-"]
            del final_list[i][1]["l-"]
        except: pass
        if nlep > 0: final_list[i][1]["lep"] = nlep
        ntau = 0
        try:
            ntau += final_list[i][1]["tau+"]
            del final_list[i][1]["tau+"]
        except: pass
        try:
            ntau += final_list[i][1]["tau-"]
            del final_list[i][1]["tau-"]
        except: pass
        if ntau > 0: final_list[i][1]["tau"] = ntau


    obj_name_dict = {"jet":"J", "bjet":"bjet", "wjet":"wjet", "zjet":"zjet", 
                     "topjet":"topjet", "higgsjet":"higgsjet", "lep":"L", "OSlep":"OSL", "SSlep":"SSL", 
                     "tau":"tau", "OStau":"OStau", "SStau":"SStau", "gam":"gam"}
    obj_list = ["jet", "bjet", "wjet", "zjet","topjet", "higgsjet", "lep", "OSlep", "SSlep", "tau", "OStau", "SStau", "gam"]
    final_dict = {}
    for val in final_list:
        xsec = val[0]
        obj_dict = val[1]
        name = ""
        for obj in obj_list:
            if obj in obj_dict: 
                if isinstance(obj_dict[obj], int):
                    name = name + str(obj_dict[obj]) + str(obj_name_dict[obj]) + "_"
                else:
                    name = name + str(obj_name_dict[obj]) + "_"
        MET = 0
        if "MET" in obj_dict:
            name = name + "MET"
            MET = obj_dict["MET"]
        else:
            name = name[:-1]
        Data = Structure()
        Data.xsec = xsec
        Data.MET = MET
        final_dict[name] = Data

    #for name, Data in final_dict.items():
    #    print name.rjust(30), Data.xsec, Data.MET

    return final_dict





def get_fp_dict(procs, prc):
    fp_dict = {}
    pmiss = []
    for pData in procs[prc].plist:
        if pData.pname in ["q", "gl"]: 
            try: fp_dict["jet"].append(pData.p)
            except: fp_dict["jet"] = [pData.p]
        if pData.pname in ["b"]: 
            try: fp_dict["bjet"].append(pData.p)
            except: fp_dict["bjet"] = [pData.p]
        if pData.pname in ["em"]: 
            try: fp_dict["lep-"].append(pData.p)
            except: fp_dict["lep-"] = [pData.p]
        if pData.pname in ["ep"]: 
            try: fp_dict["lep+"].append(pData.p)
            except: fp_dict["lep+"] = [pData.p]
        if pData.pname in ["taum"]: 
            try: fp_dict["tau-"].append(pData.p)
            except: fp_dict["tau-"] = [pData.p]
        if pData.pname in ["taup"]: 
            try: fp_dict["tau+"].append(pData.p)
            except: fp_dict["tau+"] = [pData.p]
        if pData.pname in ["t"]: 
            try: fp_dict["top"].append(pData.p)
            except: fp_dict["top"] = [pData.p]
        if pData.pname in ["tm"]: 
            try: fp_dict["top-"].append(pData.p)
            except: fp_dict["top-"] = [pData.p]
        if pData.pname in ["tp"]: 
            try: fp_dict["top+"].append(pData.p)
            except: fp_dict["top+"] = [pData.p]
        if pData.pname in ["w"]: 
            try: fp_dict["w"].append(pData.p)
            except: fp_dict["w"] = [pData.p]
        if pData.pname in ["wm"]: 
            try: fp_dict["w-"].append(pData.p)
            except: fp_dict["w-"] = [pData.p]
        if pData.pname in ["wp"]: 
            try: fp_dict["w+"].append(pData.p)
            except: fp_dict["w+"] = [pData.p]
        if pData.pname in ["z"]: 
            try: fp_dict["z"].append(pData.p)
            except: fp_dict["z"] = [pData.p]
        if pData.pname in ["h0"]: 
            try: fp_dict["h0"].append(pData.p)
            except: fp_dict["h0"] = [pData.p]
        if pData.pname in ["h2"]: 
            try: fp_dict["h2"].append(pData.p)
            except: fp_dict["h2"] = [pData.p]
        if pData.pname in ["h3"]: 
            try: fp_dict["h3"].append(pData.p)
            except: fp_dict["h3"] = [pData.p]
        if pData.pname in ["hp"]: 
            try: fp_dict["hp"].append(pData.p)
            except: fp_dict["hp"] = [pData.p]
        if pData.pname in ["hm"]: 
            try: fp_dict["hm"].append(pData.p)
            except: fp_dict["hm"] = [pData.p]
        if pData.pname in ["gam"]: 
            try: fp_dict["gam"].append(pData.p)
            except: fp_dict["gam"] = [pData.p]        
        if pData.pname in ["nu"]: pmiss.append(pData.p)
        if pData.pname in ["N1", "R32", "SNU", "NUT"] and pData.pname == procs[prc].plist[len(procs[prc].plist) - 1].pname:
            pmiss.append(pData.p)
    if not len(pmiss) == 0:
        if len(pmiss) == 1: MET = reduce(get_quad_ave, pmiss)
        if len(pmiss) > 1: MET = reduce(get_quad_ave, pmiss)
        fp_dict["MET"] = MET
    return fp_dict

def get_quad_ave(a, b): return sqrt(a**2 + b**2)







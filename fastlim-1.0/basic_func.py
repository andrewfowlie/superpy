#! /usr/bin/env python

from pyslha import *

class Structure: pass

class Paths_and_Data:
    def __init__(self, infile, blocks, decays):

        self.infile = infile
        self.blocks = blocks
        self.decays = decays

        self.Prospino_dir = "Prospino"

        self.xsec_table = {}
        self.xsec_table[("GG", 7)]  = "xsection_tables/7TeV/NLO+NLL/gg_7TeV_NLONLL.xsec"
        self.xsec_table[("GG", 8)]  = "xsection_tables/8TeV/NLO+NLL/gg_8TeV_NLONLL.xsec"
        self.xsec_table[("GG", 14)] = "xsection_tables/Prospino_NLO/gg_pros_14.xsec"

        self.xsec_table[("QG", 7)]  = "xsection_tables/7TeV/NLO+NLL/sg_7TeV_NLONLL.xsec"
        self.xsec_table[("QG", 8)]  = "xsection_tables/8TeV/NLO+NLL/sg_8TeV_NLONLL.xsec"
        self.xsec_table[("QG", 14)] = "xsection_tables/Prospino_NLO/sg_pros_14.xsec"

        self.xsec_table[("QQ", 7)]  = "xsection_tables/7TeV/NLO+NLL/ss_7TeV_NLONLL.xsec"
        self.xsec_table[("QQ", 8)]  = "xsection_tables/8TeV/NLO+NLL/ss_8TeV_NLONLL.xsec"
        self.xsec_table[("QQ", 14)] = "xsection_tables/Prospino_NLO/ss_pros_14.xsec"

        self.xsec_table[("QQbar", 7)]  = "xsection_tables/7TeV/NLO+NLL/sb_7TeV_NLONLL.xsec"
        self.xsec_table[("QQbar", 8)]  = "xsection_tables/8TeV/NLO+NLL/sb_8TeV_NLONLL.xsec"
        self.xsec_table[("QQbar", 14)] = "xsection_tables/Prospino_NLO/sb_pros_14.xsec"

        self.xsec_table[("Q3Q3", 7)]  = "xsection_tables/7TeV/SUSYxsecWG/Q3Q3bar7.xsec"
        self.xsec_table[("Q3Q3", 8)]  = "xsection_tables/8TeV/SUSYxsecWG/Q3Q3bar8.xsec"
        self.xsec_table[("Q3Q3", 13)] = "xsection_tables/13TeV/LHCxsecWG/Q3Q3bar13.xsec"
        self.xsec_table[("Q3Q3", 14)] = "xsection_tables/Prospino_NLO/tb_pros_14.xsec"

        self.xsec_table_dcpl = {}
        self.xsec_table_dcpl[("GG", 7)]     = "xsection_tables/7TeV/NLO+NLL/gdcpl_7TeV_NLONLL.xsec"
        self.xsec_table_dcpl[("GG", 8)]     = "xsection_tables/8TeV/NLO+NLL/gdcpl_8TeV_NLONLL.xsec"
        self.xsec_table_dcpl[("GG", 13)]    = "xsection_tables/LHCxsecWG/GG_Qdcpl13.xsec"
        self.xsec_table_dcpl[("QQ", 7)]     = "xsection_tables/7TeV/NLO+NLL/sdcpl_7TeV_NLONLL.xsec"
        self.xsec_table_dcpl[("QQ", 8)]     = "xsection_tables/8TeV/NLO+NLL/sdcpl_8TeV_NLONLL.xsec"
        self.xsec_table_dcpl[("QQbar", 13)] = "xsection_tables/LHCxsecWG/QQbar_Gdcpl13.xsec"


def get_mass(blocks, pname):
    massesDic = blocks["MASS"].entries
    pid = get_pid(pname)
    pid = abs(pid)
    if pid == 0:
        print "unknow particle", pname
        return
    elif pname == "Q":
        msq = 0
        for pdm in range(1000001, 1000005) + range(2000001, 2000005):
            msq += abs(massesDic[pdm])
        msq = msq/8.
        return msq
    else:
        try:
            return abs(massesDic[pid])
        except:
            return 0
    return

def flatten_tuple(array):
    res = []
    for el in array:
        if isinstance(el, tuple):
            res.extend(flatten_tuple(el))
            continue
        res.append(el)
    return res

def get_quad_ave(a, b): return sqrt(a**2 + b**2)

def is_SUSY(pid):
    if abs(pid) > 1000000:
        return True
    else:
        return False

def is_gamma(pid):
    if abs(pid) in [22]: return True
    return False

def is_neutrino(pid):
    if abs(pid) in [12, 14, 16]: return True
    return False

def is_e_or_mu(pid):
    if abs(pid) in [11, 13]: return True
    return False

def is_jet_or_tau(pid):
    if abs(pid) in [1, 2, 3, 4, 5, 15]: return True
    return False

def get_max_key_len(dic):
    maxlen = -100
    for key in dic:
        if len(key) > maxlen: maxlen = len(key)
    return maxlen

def output_dict(dic):
    space = get_max_key_len(dic)
    for key, var in dic.iteritems():
        print key.rjust(space),":",var

def get_name_simple(pid):
    name = "X"
    if pid == 1000039: name = "R32"
    if pid == 1000021: name = "G"
    if pid == 1000022: name = "N1"
    if pid == 1000023: name = "N2"
    if pid == 1000025: name = "N3"
    if pid == 1000035: name = "N4"
    # Superpy - for fifth neutralino.
    if pid == 1000045: name = "N5"
    if abs(pid) in range(1000001, 1000005): name = "Q"
    if abs(pid) in range(2000001, 2000005): name = "Q"
    if abs(pid) == 1000005: name = "B1"
    if abs(pid) == 2000005: name = "B2"
    if abs(pid) == 1000006: name = "T1"
    if abs(pid) == 2000006: name = "T2"
    if abs(pid) in [1000011, 2000011]: name = "E"
    if abs(pid) in [1000013, 2000013]: name = "M"
    if abs(pid) == 1000015: name = "TAU1"
    if abs(pid) == 2000015: name = "TAU2"
    if abs(pid) == 1000016: name = "NUT"
    if abs(pid) in [1000012, 1000014]: name = "NU"
    if abs(pid) in range(1, 5): name = "q"
    if abs(pid) == 5:  name = "b"
    if abs(pid) in [12, 14, 16]: name = "n"
    if pid == 25: name = "h"
    if pid == 35: name = "h2"
    if pid == 36: name = "h3"
    if pid == 23: name = "z"
    if pid == 22: name = "gam"
    if pid == 21: name = "g"

    if abs(pid) == 1000024: name = "C1"
    if abs(pid) == 1000037: name = "C2"
    if abs(pid) == 37: name = "hp"
    if abs(pid) == 24: name = "w"
    if abs(pid) == 6:  name = "t"
    if abs(pid) in [11]: name = "e"
    if abs(pid) in [13]: name = "m"
    if abs(pid) == 15: name = "ta"

    if name == "X":
        print "ERROR: unknown particle: ", pid
        exit()
    return name

def get_name_charge(pid):
    name = get_name_simple(pid)
    if pid == 37: name = "hp"
    if pid == -37: name = "hm"
    if pid == 24: name = "wp"
    if pid == -24: name = "wm"
    if pid == 6:  name = "tp"
    if pid == -6:  name = "tm"
    if pid in [11, 13]: name = "em"
    if pid in [-11, -13]: name = "ep"
    if pid == 15: name = "tam"
    if pid == -15: name = "tap"
    return name

def all_particles():
    part_list = ["G", "N1", "N2", "N3", "N4", "Q", "B1", "B2", "T1", "T2", "R32", "E", "M", "TAU1", "TAU2"]
    part_list += ["NUT", "NU", "g", "gam", "z", "h3", "h2", "h0", "nu", "b", "q"]
    part_list += ["C1", "C2", "w", "hp", "ta", "e", "m", "t", "C1p", "C1m", "C2p", "C2m", "wp", "wm", "hp", "hm"]
    part_list += ["tap", "tam", "ep", "em", "mm", "mp", "tm", "tp"]
    return part_list

def get_pid(name):
    pid = 0
    if name == "G":  pid = 1000021
    if name == "N1": pid = 1000022
    if name == "N2": pid = 1000023
    if name == "N3": pid = 1000025
    if name == "N4": pid = 1000035
    if name == "Q":  pid = 1000001
    if name == "B1": pid = 1000005
    if name == "B2": pid = 2000005
    if name == "T1": pid = 1000006
    if name == "T2": pid = 2000006
    if name == "R32": pid = 1000039
    if name == "E": pid = 1000011
    if name == "M": pid = 1000013
    if name == "TAU1": pid = 1000015
    if name == "TAU2": pid = 2000015
    if name == "NUT": pid = 1000016
    if name == "NU": pid = 1000012
    if name == "g": pid = 21
    if name == "gam": pid = 22
    if name == "z": pid = 23
    if name == "h3": pid = 36
    if name == "h2": pid = 35
    if name == "h": pid = 25
    if name == "n": pid = 12
    if name == "b": pid = 5
    if name == "q": pid = 1

    if name == "C1": pid = 1000024
    if name == "C2": pid = 1000037
    if name == "w": pid = 24
    if name == "hp": pid = 37
    if name == "ta": pid = 15
    if name == "e": pid = 11
    if name == "m": pid = 13
    if name == "t": pid = 6

    if name == "C1p": pid = 1000024
    if name == "C1m": pid = -1000024
    if name == "C2p": pid = 1000037
    if name == "C2m": pid = -1000037
    if name == "wp": pid = 24
    if name == "wm": pid = -24
    if name == "hm": pid = -37
    if name == "hp": pid = 37
    if name == "tap": pid = -15
    if name == "tam": pid = 15
    if name == "ep": pid = -11
    if name == "em": pid = 11
    if name == "mp": pid = -13
    if name == "mm": pid = 13
    if name == "tm": pid = -6
    if name == "tp": pid = 6

    if pid == 0: print "unknown", name
    return pid


def full_decays(decays):
    dummy_dec = decays.copy()
    for pid in dummy_dec:
        pneut = [1000021, 1000022, 1000023, 1000025, 1000035, 1000039]
        pneut += [1000012, 1000014, 1000016]
        pneut += [21, 22, 23, 25, 35, 36, 39]
        if not pid in pneut:
            decays[-pid] = Particle( -pid, decays[pid].totalwidth, decays[pid].mass )
            for mode in dummy_dec[pid].decays:
                ids_conj = []
                for ids in mode.ids:
                    if ids in pneut:
                        ids_conj.append(ids)
                    else:
                        ids_conj.append(-ids)
                decays[-pid].add_decay(mode.br, mode.nda, ids_conj)
    return decays

def extract_numbers(line):
    vals = []
    for elem in line.split():
        try:
            vals.append(float(elem))
        except: pass
    return vals

def output_warning(err_set):
    worn_list = []
    for mess, var in err_set.iteritems():
        if var == 0: worn_list.append(mess)
    if len(worn_list) == 0: return
    print "#--------------  Warnings  ----------------#"
    for mess in worn_list:
        print mess
    print ""
    #print "#--------------------------#"

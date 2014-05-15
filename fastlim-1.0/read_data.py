#! /usr/bin/env python

# SuperPy - hacked to fix relative filenames.
# I import fast-lim from an external directory.

import os.path
from math import *
from basic_func import *
import sys
import re
sys.path.append('interpolation')
from interpolate import *
sys.path.append('statistics/CLs_limits')
from mycls import *
#from PoissonModelWithBackg import *
from extract_chain import *

def read_2D_xsec(filename, m1, m2):

    # SuperPy hack to fix relative filenames.
    fastlim_dir = os.path.dirname(os.path.realpath(__file__))
    filename = os.path.join(fastlim_dir, filename)

    if not os.path.isfile(filename):
        print ""
        print '!! File "' + str(filename) + '" not found !!'
        print ""
        return 0
    xsec_list = []
    for line in open(filename, 'r'):
        vals = extract_numbers(line)
        try:
            MM1 = float(vals[0])
            MM2 = float(vals[1])
            xsec = float(vals[2])
            logx = log(xsec)
            #print MM1, MM2, logx
            if not xsec == 0: xsec_list.append([MM1, MM2, logx])
        except: pass
    #print m1, m2
    eps = 1.0
    logxsec = Interpolate2D(xsec_list, [m1 + eps, m2 + eps])
    #print logxsec
    #print m1, m2, logxsec
    if isinstance(logxsec, float):
        xsec = exp(logxsec)
    else:
        xsec = 0
        warn = "In " + str(filename) + " xsec is set by 0 because of the point is outside our table."
        print warn
    return xsec

def read_1D_xsec(filename, m1):

    # SuperPy hack to fix relative filenames.
    fastlim_dir = os.path.dirname(os.path.realpath(__file__))
    filename = os.path.join(fastlim_dir, filename)

    if not os.path.isfile(filename):
        print ""
        print '!! File "' + str(filename) + '" not found !!'
        print ""
        return 0
    xsec_list = []
    for line in open(filename, 'r'):
        vals = extract_numbers(line)
        try:
            MM1 = float(vals[0])
            xsec = float(vals[1])
            #print MM1, xsec
            if not xsec == 0: xsec_list.append([MM1, log(xsec)])
        except: pass
    #for x in xsec_list: print x, "a"
    logxsec = Interpolate1D(xsec_list, m1)
    if isinstance(logxsec, float):
        xsec = exp(logxsec)
        #print xsec
    else:
        xsec = 0
        warn = "In " + str(filename) + " xsec is set by 0 because of the point is outside our table."
        print warn
    return xsec


def get_precalc_xsec(mode, root_S, Input):

    root_S = int(root_S)

    if mode in ["T1T1", "T2T2", "B1B1", "B2B2"]:
        if mode == "T1T1": m1 = get_mass(Input.blocks, "T1")
        if mode == "T2T2": m1 = get_mass(Input.blocks, "T2")
        if mode == "B1B1": m1 = get_mass(Input.blocks, "B1")
        if mode == "B2B2": m1 = get_mass(Input.blocks, "B2")
        grid_file = Input.xsec_table[("Q3Q3", root_S)]
        xsec = read_1D_xsec(grid_file, m1)
        return xsec

    mQ = get_mass(Input.blocks, "Q")
    mG = get_mass(Input.blocks, "G")
    #line_mass = "(mQ, mG) = (" + str(mQ) +", "+ str(mG) +")"

    grid_file = Input.xsec_table[(mode, root_S)]
    xsec = read_2D_xsec(grid_file, mQ, mG)
    if xsec == xsec and xsec > 0:
        return xsec

    try:
        grid_file = Input.xsec_table_dcpl[(mode, root_S)]
        xsec = read_1D_xsec(grid_file, mG)
        if xsec == xsec and xsec > 0:
            if mode == "GG": line = "GG xsec: out side of the range => read Q decoupled GG xsec"
            if mode == "QQ": line = "QQ xsec: out side of the range => read G decoupled QQ xsec"
            print line
            return xsec
    except: pass

    print str(mode) + " xsec: The xsec table is not available => use Prospino"

    return 0


def get_xsec(par1, par2, root_S, Input, nlo = False):

    root_S = int(root_S)

    xsec = 0
    if par1 == "G" and par2 == "G": xsec = get_precalc_xsec("GG", root_S, Input)
    if par1 == "Q" and par2 == "Q": xsec = get_precalc_xsec("QQ", root_S, Input)
    if (par1 == "Q" and par2 == "Qbar") or (par1 == "Qbar" and par2 == "Q"):
        xsec = get_precalc_xsec("QQbar", root_S, Input)
    if (par1 == "Q" and par2 == "G") or (par1 == "G" and par2 == "Q"):
        xsec = get_precalc_xsec("QG", root_S, Input)

    if par1 == "T1" and par2 == "T1" and root_S: xsec = get_precalc_xsec("T1T1", root_S, Input)
    if par1 == "T2" and par2 == "T2" and root_S: xsec = get_precalc_xsec("T2T2", root_S, Input)
    if par1 == "B1" and par2 == "B1" and root_S: xsec = get_precalc_xsec("B1B1", root_S, Input)
    if par1 == "B2" and par2 == "B2" and root_S: xsec = get_precalc_xsec("B2B2", root_S, Input)

    if xsec != xsec or xsec == 0:
        mass1 = get_mass(Input.blocks, par1)
        mass2 = get_mass(Input.blocks, par2)
        if max(mass1,mass2) > root_S:
            print "Decoupling " + str(par1) + ": " + str(mass1) + ";" + str(par2)  + ": " + str(mass2)
            xsec = 10**(-13)
        else:
            xsec = prospino_xsec(par1, par2, root_S, Input, nlo)

    return xsec

def get_cls(Data):
    if Data.Nobs == 0:
        cls= "Nobs=0: cannot calculate"
    else:
        cls = mycls(Data.nvis, 0, Data.Nbg, Data.errNbg, Data.Nobs)
    return cls

def get_clsexp(Data):
    cls = mycls(Data.nvis, 0, Data.Nbg, Data.errNbg, Data.Nbg)
    return cls

def get_proc3D_dict(proc_list, procs, m1, m2, m3):
    proc_dict = {}
    for proc in proc_list:
        try:
            Data = Structure()
            Data.xsec = procs[proc].xsec
            Data.m1 = m1
            Data.m2 = m2
            Data.m3 = m3
            proc_dict[proc] = Data
        except: pass
    return proc_dict

def get_proc2D_dict(proc_list, procs, m1, m2):
    proc_dict = {}
    for proc in proc_list:
        try:
            Data = Structure()
            Data.xsec = procs[proc].xsec
            Data.m1 = m1
            Data.m2 = m2
            proc_dict[proc] = Data
        except: pass
    return proc_dict


def get_ana_dict(ana_list, RootS):
    ana_dict = {}
    for ana in ana_list:

        dir_name = "analyses_info/" + str(RootS) + "TeV/" + ana

        # SuperPy hack to fix relative filenames.
        fastlim_dir = os.path.dirname(os.path.realpath(__file__))
        dir_name = os.path.join(fastlim_dir, dir_name)

        if not os.path.exists(dir_name):
            print dir_name, "does not exist!  skip"
            continue
        ana_info = ''
        for line in open(dir_name + "/info.txt", 'r'):
            ana_info += line
        #print ana_info
        sr_dict = {}
        root_S = 0
        lumi = 0
        for line in open(dir_name + "/SR_info.txt", 'r'):
            elem = line.split()
            try:
                root_S = int(elem[0])
                if not int(root_S) == int(RootS):
                    print "Luminosity Error!!"
                    sys.exit()
                lumi   = float(elem[1])
                Data = Structure()
                Data.Nobs   = int(elem[2])
                Data.Nbg    = float(elem[3])
                Data.errNbg = float(elem[4])
                Data.UL_nvis_obs = float(elem[5])
                Data.UL_nvis_exp = float(elem[6])
                Data.UL_xvis_obs = float(elem[7])
                Data.UL_xvis_exp = float(elem[8])
                if Data.UL_nvis_obs < 0 and Data.UL_xvis_obs > 0: Data.UL_nvis_obs = Data.UL_xvis_obs * lumi
                if Data.UL_nvis_exp < 0 and Data.UL_xvis_exp > 0: Data.UL_nvis_exp = Data.UL_xvis_exp * lumi
                sr_name = elem[9]
                Data.sr_info = ''
                for idm in range(10, len(elem)):
                    Data.sr_info += elem[idm]
                    if not idm+1 == len(elem): Data.sr_info +=' '
                if Data.Nbg > 0:  sr_dict[sr_name] = Data
            except: pass
        aData = Structure()
        aData.ana_info = ana_info
        aData.root_S = root_S
        aData.lumi = lumi
        aData.sr_dict = sr_dict
        ana_dict[ana] = aData
    return ana_dict


def get_proc_Data_3D(proc3D_dict, ana, sr_name, root_S, lumi):
    proc_Data = {}
    warn_list = []
    min_eff = 1./100000
    for proc, Data in proc3D_dict.iteritems():
        procname = proc
        if proc in ["GtT2bN1_GtT2bN1", "GtT1bN1_GtT2bN1"]: procname = "GtT1bN1_GtT1bN1"
        if proc in ["GtT2tN1_GtT2tN1", "GtT1tN1_GtT2tN1"]: procname = "GtT1tN1_GtT1tN1"
        if proc in ["GtT2bN1_GtT2tN1", "GtT1bN1_GtT2tN1", "GtT1tN1_GtT2bN1"]: procname = "GtT1bN1_GtT1tN1"
        if proc in ["GbB2bN1_GbB2bN1", "GbB1bN1_GbB2bN1"]: procname = "GbB1bN1_GbB1bN1"
        if proc in ["GbB2tN1_GbB2tN1", "GbB1tN1_GbB2tN1"]: procname = "GbB1tN1_GbB1tN1"
        if proc in ["GbB2bN1_GbB2tN1", "GbB1bN1_GbB2tN1", "GbB1tN1_GbB2bN1"]: procname = "GbB1bN1_GbB1tN1"
        #print proc, ana
        dir_name = "efficiency_tables/" + procname + "/" + str(root_S) + "TeV/" + ana
        eff_file = ""
        for i in range(0,50):
            eff_file = dir_name + "/ana_" + str(i) + "_" + sr_name + ".effi"

            # SuperPy hack to fix relative filenames.
            fastlim_dir = os.path.dirname(os.path.realpath(__file__))
            eff_file = os.path.join(fastlim_dir, eff_file)

            if os.path.exists(eff_file): break
        if eff_file == "":
            print eff_file, "does not exist!!  skep"
            continue
        logefftab = []
        eff = -1
        for line in open(eff_file, 'r'):
            vals = extract_numbers(line)
            try:
                m1 = vals[0]
                m2 = vals[1]
                m3 = vals[2]
                eff = vals[3]
                if eff == 0:
                    logeff = log(min_eff)
                else:
                    logeff = log(eff)
                logefftab.append([m1, m2, m3, logeff])
            except: pass
        #print m1, m2, m3, logefftab
        LogEff = Interpolate3D(logefftab, [Data.m1, Data.m2, Data.m3], err_out = "off")
        #LogEff = 0
        prcData = Structure()
        prcData.eff = 0
        if not isinstance(LogEff, float):
            warn = "Warning in " + str(proc) + ": Efficiency in " + str(ana) + ", " + str(sr_name) + " is outside our table. => [eff is set at 0]"
            warn_list.append(warn)
        else:
            if exp(LogEff) > 1:
                warn = "Warning in " + str(proc) + ", " + str(ana) + ", " + str(sr_name) + ": eff > 1 in interpolation. => [eff is set at 1]"
                warn_list.append(warn)
                prcData.eff = 1
            else:
                prcData.eff = exp(LogEff)
        prcData.xsec = Data.xsec
        prcData.xvis = prcData.eff * Data.xsec
        prcData.nvis = prcData.xvis * lumi
        proc_Data[proc] = prcData
    #print proc_Data.keys()
    return proc_Data, warn_list


def get_proc_Data_2D(proc2D_dict, ana, sr_name, root_S, lumi):
    proc_Data = {}
    warn_list = []
    for proc, Data in proc2D_dict.iteritems():
        procname = proc
        if proc in ["T2bN1_T2bN1", "B1bN1_B1bN1", "B2bN1_B2bN1"]: procname = "T1bN1_T1bN1"
        if proc in ["T2bN1_T2tN1", "B1bN1_B1tN1", "B2bN1_B2tN1"]: procname = "T1bN1_T1tN1"
        if proc in ["T2tN1_T2tN1", "B1tN1_B1tN1", "B2tN1_B2tN1"]: procname = "T1tN1_T1tN1"
        dir_name = "efficiency_tables/" + procname + "/" + str(root_S) + "TeV/" + ana
        eff_file = ""
        for i in range(0,50):
            eff_file = dir_name + "/ana_" + str(i) + "_" + sr_name + ".effi"

            # SuperPy hack to fix relative filenames.
            fastlim_dir = os.path.dirname(os.path.realpath(__file__))
            eff_file = os.path.join(fastlim_dir, eff_file)

            if os.path.exists(eff_file): break
        if eff_file == "":
            print eff_file, "does not exist!!  skep"
            continue
        logefftab = []
        eff = -1
        for line in open(eff_file, 'r'):
            vals = extract_numbers(line)
            try:
                m1 = vals[0]
                m2 = vals[1]
                eff = vals[2]
                if eff > 0:
                    logeff = log(eff)
                else:
                    logeff = -100000000000000.
                logefftab.append([m1, m2, logeff])
            except: pass
        #print logefftab
        Dm1 = Data.m1
        Dm2 = Data.m2
        if Dm2 < 1.: Dm2 = 1.1
        LogEff = Interpolate2D(logefftab, [Dm1, Dm2])
        prcData = Structure()
        prcData.eff = 0
        if not isinstance(LogEff, float):
            warn = "Warning in " + str(proc) + ": Efficiency in " + str(ana) + ", " + str(sr_name) + " is outside our table. => [eff is set at 0]"
            warn_list.append(warn)
        else:
            if exp(LogEff) > 1:
                warn = "Warning in " + str(proc) + ", " + str(ana) + ", " + str(sr_name) + ": eff > 1 in interpolation. => [eff is set at 1]"
                warn_list.append(warn)
                prcData.eff = 1
            else:
                prcData.eff = exp(LogEff)
        prcData.xsec = Data.xsec
        prcData.xvis = prcData.eff * Data.xsec
        prcData.nvis = prcData.xvis * lumi
        proc_Data[proc] = prcData
        #print ana, sr_name, proc, prcData.xsec, prcData.eff
    #print ana, proc_Data.keys()
    return proc_Data, warn_list



def get_proc_Data_3D_scan(proc3D_dict, ana, sr_name, root_S, lumi, eff_list_dict):
    proc_Data = {}
    warn_list = []
    min_eff = 1./100000
    for proc, Data in proc3D_dict.iteritems():
        procname = proc
        if proc in ["GtT2bN1_GtT2bN1", "GtT1bN1_GtT2bN1"]: procname = "GtT1bN1_GtT1bN1"
        if proc in ["GtT2tN1_GtT2tN1", "GtT1tN1_GtT2tN1"]: procname = "GtT1tN1_GtT1tN1"
        if proc in ["GtT2bN1_GtT2tN1", "GtT1bN1_GtT2tN1", "GtT1tN1_GtT2bN1"]: procname = "GtT1bN1_GtT1tN1"
        if proc in ["GbB2bN1_GbB2bN1", "GbB1bN1_GbB2bN1"]: procname = "GbB1bN1_GbB1bN1"
        if proc in ["GbB2tN1_GbB2tN1", "GbB1tN1_GbB2tN1"]: procname = "GbB1tN1_GbB1tN1"
        if proc in ["GbB2bN1_GbB2tN1", "GbB1bN1_GbB2tN1", "GbB1tN1_GbB2bN1"]: procname = "GbB1bN1_GbB1tN1"
        key = (procname, ana, sr_name)
        logefftab = eff_list_dict[key]
        Dm1 = Data.m1
        Dm2 = Data.m2
        Dm3 = Data.m3
        if Dm3 < 1.: Dm3 = 1.1
        LogEff = Interpolate3D(logefftab, [Dm1, Dm2, Dm3], err_out = "off")
        #LogEff = 0
        prcData = Structure()
        prcData.eff = 0
        if not isinstance(LogEff, float):
            warn = "Warning in " + str(proc) + ": Efficiency in " + str(ana) + ", " + str(sr_name) + " is outside our table. => [eff is set at 0]"
            warn_list.append(warn)
        else:
            if exp(LogEff) > 1:
                warn = "Warning in " + str(proc) + ", " + str(ana) + ", " + str(sr_name) + ": eff > 1 in interpolation. => [eff is set at 1]"
                warn_list.append(warn)
                prcData.eff = 1
            else:
                prcData.eff = exp(LogEff)
        prcData.xsec = Data.xsec
        prcData.xvis = prcData.eff * Data.xsec
        prcData.nvis = prcData.xvis * lumi
        proc_Data[proc] = prcData
    #print proc_Data.keys()
    return proc_Data, warn_list


def get_proc_Data_2D_scan(proc2D_dict, ana, sr_name, root_S, lumi, eff_list_dict):
    proc_Data = {}
    warn_list = []
    for proc, Data in proc2D_dict.iteritems():
        procname = proc
        if proc in ["T2bN1_T2bN1", "B1bN1_B1bN1", "B2bN1_B2bN1"]: procname = "T1bN1_T1bN1"
        if proc in ["T2bN1_T2tN1", "B1bN1_B1tN1", "B2bN1_B2tN1"]: procname = "T1bN1_T1tN1"
        if proc in ["T2tN1_T2tN1", "B1tN1_B1tN1", "B2tN1_B2tN1"]: procname = "T1tN1_T1tN1"
        #print logefftab
        key = (procname, ana, sr_name)
        logefftab = eff_list_dict[key]
        Dm1 = Data.m1
        Dm2 = Data.m2
        if Dm2 < 1.: Dm2 = 1.1
        LogEff = Interpolate2D(logefftab, [Dm1, Dm2])
        prcData = Structure()
        prcData.eff = 0
        if not isinstance(LogEff, float):
            warn = "Warning in " + str(proc) + ": Efficiency in " + str(ana) + ", " + str(sr_name) + " is outside our table. => [eff is set at 0]"
            warn_list.append(warn)
        else:
            if exp(LogEff) > 1:
                warn = "Warning in " + str(proc) + ", " + str(ana) + ", " + str(sr_name) + ": eff > 1 in interpolation. => [eff is set at 1]"
                warn_list.append(warn)
                prcData.eff = 1
            else:
                prcData.eff = exp(LogEff)
        prcData.xsec = Data.xsec
        prcData.xvis = prcData.eff * Data.xsec
        prcData.nvis = prcData.xvis * lumi
        proc_Data[proc] = prcData
        #print ana, sr_name, proc, prcData.xsec, prcData.eff
    #print ana, proc_Data.keys()
    return proc_Data, warn_list




def get_results(ana_dict, proc3D_dict, proc2D_dict):
    results = {}
    warning_list = []
    for ana, anaData in ana_dict.iteritems():
        for sr_name, srData in anaData.sr_dict.iteritems():
            #print key
            proc_Data = {}
            pData_dict, warn_list = get_proc_Data_2D(proc2D_dict, ana, sr_name, anaData.root_S, anaData.lumi)
            proc_Data.update(pData_dict)
            warning_list += warn_list
            pData_dict, warn_list = get_proc_Data_3D(proc3D_dict, ana, sr_name, anaData.root_S, anaData.lumi)
            proc_Data.update(pData_dict)
            warning_list += warn_list

            #  make "results" dictionary
            key = (ana, sr_name)
            Data = Structure()
            Data.root_S = anaData.root_S
            Data.lumi = anaData.lumi
            Data.analysis = ana
            Data.SR = sr_name
            Data.ana_info = anaData.ana_info
            Data.SR_info = srData.sr_info
            Data.Nobs = srData.Nobs
            Data.Nbg = srData.Nbg
            Data.errNbg = srData.errNbg
            Data.UL_nvis_obs = srData.UL_nvis_obs
            Data.UL_nvis_exp = srData.UL_nvis_exp
            Data.UL_xvis_obs = srData.UL_xvis_obs
            Data.UL_xvis_exp = srData.UL_xvis_exp
            xvis_tot = 0
            nvis_tot = 0
            for proc, prcData in proc_Data.iteritems():
                if prcData.xvis > 0: xvis_tot += prcData.xvis
                if prcData.nvis > 0: nvis_tot += prcData.nvis
            #print xvis_tot, nvis_tot
            Data.xvis = xvis_tot
            Data.nvis = nvis_tot
            Data.Rvis_obs = -1
            Data.Rvis_exp = -1
            if Data.UL_nvis_obs > 0:  Data.Rvis_obs = nvis_tot / Data.UL_nvis_obs
            if Data.UL_nvis_exp > 0:  Data.Rvis_exp = nvis_tot / Data.UL_nvis_exp
            if Data.UL_xvis_obs > 0:  Data.Rvis_obs = xvis_tot / Data.UL_xvis_obs
            if Data.UL_xvis_exp > 0:  Data.Rvis_exp = xvis_tot / Data.UL_xvis_exp
            Data.proc_Data = proc_Data
            results[key] = Data
    return results, warning_list



def get_results_scan(ana_dict, proc3D_dict, proc2D_dict, eff_list_dict):
    results = {}
    warning_list = []
    for ana, anaData in ana_dict.iteritems():
        for sr_name, srData in anaData.sr_dict.iteritems():
            #print key
            proc_Data = {}
            pData_dict, warn_list = get_proc_Data_2D_scan(proc2D_dict, ana, sr_name, anaData.root_S, anaData.lumi, eff_list_dict)
            proc_Data.update(pData_dict)
            warning_list += warn_list
            pData_dict, warn_list = get_proc_Data_3D_scan(proc3D_dict, ana, sr_name, anaData.root_S, anaData.lumi, eff_list_dict)
            proc_Data.update(pData_dict)
            warning_list += warn_list

            #  make "results" dictionary
            key = (ana, sr_name)
            Data = Structure()
            Data.root_S = anaData.root_S
            Data.lumi = anaData.lumi
            Data.analysis = ana
            Data.SR = sr_name
            Data.ana_info = anaData.ana_info
            Data.SR_info = srData.sr_info
            Data.Nobs = srData.Nobs
            Data.Nbg = srData.Nbg
            Data.errNbg = srData.errNbg
            Data.UL_nvis_obs = srData.UL_nvis_obs
            Data.UL_nvis_exp = srData.UL_nvis_exp
            Data.UL_xvis_obs = srData.UL_xvis_obs
            Data.UL_xvis_exp = srData.UL_xvis_exp
            xvis_tot = 0
            nvis_tot = 0
            for proc, prcData in proc_Data.iteritems():
                if prcData.xvis > 0: xvis_tot += prcData.xvis
                if prcData.nvis > 0: nvis_tot += prcData.nvis
            #print xvis_tot, nvis_tot
            Data.xvis = xvis_tot
            Data.nvis = nvis_tot
            Data.Rvis_obs = -1
            Data.Rvis_exp = -1
            if Data.UL_nvis_obs > 0:  Data.Rvis_obs = nvis_tot / Data.UL_nvis_obs
            if Data.UL_nvis_exp > 0:  Data.Rvis_exp = nvis_tot / Data.UL_nvis_exp
            if Data.UL_xvis_obs > 0:  Data.Rvis_obs = xvis_tot / Data.UL_xvis_obs
            if Data.UL_xvis_exp > 0:  Data.Rvis_exp = xvis_tot / Data.UL_xvis_exp
            Data.proc_Data = proc_Data
            results[key] = Data
    return results, warning_list


def get_eff_list2D(proclist_2D, ana_dict):

    eff_list_dict = {}
    min_eff = 1./50000/3
    for proc in proclist_2D:
        for ana in ana_dict:
            for sr_name in ana_dict[ana].sr_dict:
                key = (proc, ana, sr_name)
                root_S = ana_dict[ana].root_S
                dir_name = "efficiency_tables/" + proc + "/" + str(root_S) + "TeV/" + ana
                eff_file = ""
                for i in range(0,50):
                    eff_file = dir_name + "/ana_" + str(i) + "_" + sr_name + ".effi"
                    if os.path.exists(eff_file): break
                #print eff_file
                logefftab = []
                eff = -1
                for line in open(eff_file, 'r'):
                    vals = extract_numbers(line)
                    try:
                        m1 = float(vals[0])
                        m2 = float(vals[1])
                        eff = float(vals[2])
                        if eff == 0:
                            logeff = log(min_eff)
                        else:
                            logeff = log(eff)
                        logefftab.append([m1, m2, logeff])
                    except: pass
                eff_list_dict[key] = logefftab
    return eff_list_dict



def get_eff_list3D(proclist_3D, ana_dict):

    eff_list_dict = {}
    min_eff = 1./50000/3
    for proc in proclist_3D:
        for ana in ana_dict:
            for sr_name in ana_dict[ana].sr_dict:
                key = (proc, ana, sr_name)
                root_S = ana_dict[ana].root_S
                dir_name = "efficiency_tables/" + proc + "/" + str(root_S) + "TeV/" + ana
                eff_file = ""
                for i in range(0,50):
                    eff_file = dir_name + "/ana_" + str(i) + "_" + sr_name + ".effi"
                    if os.path.exists(eff_file): break
                #print eff_file
                logefftab = []
                eff = -1
                for line in open(eff_file, 'r'):
                    vals = extract_numbers(line)
                    try:
                        m1 = vals[0]
                        m2 = vals[1]
                        m3 = vals[2]
                        eff = vals[3]
                        if eff == 0:
                            logeff = log(min_eff)
                        else:
                            logeff = log(eff)
                        logefftab.append([m1, m2, m3, logeff])
                    except: pass
                eff_list_dict[key] = logefftab
    return eff_list_dict



def get_filelist(barelist):
    barelist.pop(0)
    outputdir = barelist.pop(-1)
    inputdir = ''
    pointlist = []
    for line in barelist:
        point = line.split("/")[-1]
        pointlist.append(point)
        if len(line.split("/")) > 1: inputdir = line.split("/" + point)[0]
    if not os.path.exists(inputdir):
        print inputdir, 'does not exist'
        print './ScanPoints.py [inputs] [output directory]'
        exit()
    if not os.path.exists(outputdir):
        print outputdir, 'does not exist'
        print './ScanPoints.py [inputs] [output directory]'
        exit()
    return inputdir, outputdir, pointlist

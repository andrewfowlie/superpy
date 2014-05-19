#! /usr/bin/env python

"""
"""

__author__ = "Kazuki Sakurai <kazuki.sakurai@desy.de>"
__version__ = "0.1"

import os.path
from math import *
from basic_func import *
from pyslha import *
import sys
#sys.path.append('statistics/CLs_limits')
#from mycls import *
#from PoissonModelWithBackg import *

###########################################

def soft_removal(procs, Ejet_thres = 20., Elep_thres = 10.):
    procs_org = procs
    #i = 0
    new_procs = {}
    for proc, data_org in procs_org.iteritems():
        #if i > 10: break
        #i += 1
        dlist = data_org.plist
        #print proc, len(dlist)
        pstr, dstr, chains, newplist0, newplists = "", "", [], [], []
        visible = False
        for pData in dlist:
            if pData.pname == "-":
                if dstr in all_particles() and is_SUSY(get_pid(dstr)): # dstr = LSP
                    pstr += dstr
                    newplist0 += declist
                    newplist0.append(pData)
                    newplists.append(newplist0)
                    newplist0 = []
                if visible:
                    pstr += dstr   # case for RPV
                    newplist0 += declist
                    newplist0.append(pData)
                    newplists.append(newplist0)
                    newplist0 = []
                chains.append(pstr)
                pstr = ""
            elif is_SUSY(pData.pid):
                if visible:
                    pstr += dstr
                    newplist0 += declist
                dstr = pData.pname
                declist = [pData]
                visible = False
            elif not is_SUSY(pData.pid):
                cut = 0
                e = sqrt(pData.mass**2 + pData.p**2)
                if is_gamma(pData.pid)      and e < Elep_thres: cut = 1
                if is_e_or_mu(pData.pid)    and e < Elep_thres: cut = 1
                if is_neutrino(pData.pid)   and e < Elep_thres: cut = 1
                if is_jet_or_tau(pData.pid) and e < Ejet_thres: cut = 1
                if cut == 0:  # visible particle is found
                    visible = True
                    dstr += pData.pname
                    declist.append(pData)
            #print str(pData.pid).rjust(15), str(pData.pname).rjust(5), str(visible).rjust(7),
            #print str(pData.p).rjust(15), dstr.rjust(10), pstr, "|",
            #if len(newplist0) > 0:
            #    for qData in newplist0: print qData.pname,
            #print ""
        if chains[0] < chains[1]: icmin, icmax = 0, 1
        if chains[1] < chains[0]: icmin, icmax = 1, 0
        process = chains[icmin] + "_" + chains[icmax]
        plist = newplists[icmin] + newplists[icmax]
        #print process, "##",
        #for pData in plist: print pData.pname,
        #print ""
        if process in new_procs:
            new_procs[process].xsec += data_org.xsec
        else:
            Data = Structure()
            Data.xsec = data_org.xsec
            Data.plist = plist
            new_procs[process] = Data
    set_rates(new_procs)
    return new_procs

def Replace(procs, target, replace):
    proc_data = procs
    dummy_data = procs.copy()
    ic = 0
    for proc in dummy_data:
        ic += 1
        #if ic > 1000: exit()
        proc_old = proc
        proc_new = proc
        found = False
        string_list = ("0" + proc_old + "0").split( target )
        if len(string_list) < 2: continue
        found = True
        #print string_list, len(string_list)
        if len(string_list) == 2: proc_new = string_list[0] + replace + string_list[1]
        if len(string_list) == 3: proc_new = string_list[0] + replace + string_list[1] + replace + string_list[2]
        if len(string_list) > 3:
            print "ERROR in sms_reassign: too many targets are found"
            exit()
        proc_new = proc_new[1:-1]
        if found:
            double_proc = (proc_new).split("_")
            #print proc_old, double_proc, proc_new
            proc_new = min(double_proc[0], double_proc[1]) +"_"+ max(double_proc[0], double_proc[1])
            org_data = proc_data.pop(proc_old)
            if proc_new in proc_data:
                proc_data[proc_new].xsec += org_data.xsec
            else:
                proc_data[proc_new] = org_data
    procs = proc_data
    set_rates(procs)
    return procs

class TheProdMode(object):

    def __init__(self):
        self.modes = {}
        self.processes = {}
        self.prod_xsec = {}
        self.xsec_tot = 0

    def scale(self, prod, prod_xsec):
        for proc, sig in self.modes[prod].iteritems():
            self.modes[prod][proc].xsec *= prod_xsec
        self.processes.update(self.modes[prod])

    def add_prod(self, part_dict, nameA, nameB, prod_xsec):
        pb2fb = 1000.
        prod_xsec = pb2fb * prod_xsec
        name1 = min(nameA, nameB)
        name2 = max(nameA, nameB)
        Part1 = part_dict[name1]
        Part2 = part_dict[name2]
        prod = name1 + "_" + name2
        proc_dict = {}
        #----------------------------
        for ch1, Data1 in Part1.chain.iteritems():
            for ch2, Data2 in Part2.chain.iteritems():
                proc = min(ch1, ch2) + "_" + max(ch1, ch2)
                sig = Data1.Br * Data2.Br
                if proc in proc_dict:
                    proc_dict[proc].xsec += sig
                else:
                    Data = Structure()
                    Data.xsec = sig
                    Data.plist = Data1.plist + get_plist_empty() + Data2.plist + get_plist_empty()
                    proc_dict[proc] = Data

        #----------------------------
        self.modes[prod] = proc_dict
        self.scale(prod, prod_xsec)
        #self.processes.update(proc_dict)
        self.prod_xsec[prod] = prod_xsec
        self.xsec_tot += prod_xsec
        for proc in self.modes[prod]:
            self.modes[prod][proc].rate = self.modes[prod][proc].xsec/prod_xsec

###########################################

def set_rates(procs, ignore = []):
    xvis = 0
    for proc in procs:
        if proc not in ignore:
            xvis += procs[proc].xsec
    for proc in procs:
        procs[proc].rate = 0
        if proc not in ignore: procs[proc].rate = procs[proc].xsec/xvis

def get_plist_empty():
    plist = []
    pData = Structure()
    pData.pname = "-"
    pData.mass = 0
    pData.pid = 0
    pData.p, pData.v = 0, 0
    plist.append(pData)
    return plist

def get_lsp(decays):
    lsp_id = 0
    n_stables = 0
    ip_stables = []
    worn = ""
    mmin = 1000000000000000000000.
    for spar in decays:
        if not is_SUSY(spar): continue
        if abs(decays[spar].mass) < mmin:
            mmin = abs(decays[spar].mass)
            lsp_id = abs(spar)
    #print lsp_id, mmin, decays[lsp_id].totalwidth
    if decays[lsp_id].totalwidth > 0: lsp_id = 0
    return lsp_id, worn


def get_lsp_old(decays):
    lsp_id = 0
    n_stables = 0
    ip_stables = []
    worn = ""
    for spar in decays:
        if is_SUSY(spar) and len(decays[spar].decays) == 0:
            n_stables += 1
            ip_stables.append(spar)
            lsp_id = spar
    if n_stables > 1:
        worn = "There are " + str(n_stables) + " stable SUSY particles !!"
        print worn
        for ip in ip_stables: print ip
    return lsp_id, worn

def extract_chain(initial_part_list, Input, charge_flag = False):

    cf = charge_flag
    blocks = Input.blocks
    decays = Input.decays
    decays = full_decays(decays)

    if charge_flag:
        get_name = get_name_charge
    else:
        get_name = get_name_simple

    part_dict = {}
    err_dict = {}

    mass_dic = blocks["MASS"].entries

    lsp_id, warning = get_lsp(decays)
    if len(warning) > 0: err_dict[warning] = 0

    if "Q" in initial_part_list:
        ###########################
        ##  for squark
        ###########################
        sq_list = []

        msq = 0
        for ipart in range(1000001, 1000005) + range(2000001, 2000005):
            msq += mass_dic[ipart]
        msq = msq/8.

        for ipart in range(1000001, 1000005) + range(2000001, 2000005):

            #print ipart
            Sq = TheParticle(ipart, lsp_id, msq, cf)
            Sq_Data = Sq.chain.pop(get_name(ipart))
            sigtot = Sq_Data.Br
            oldname = get_name(ipart)
            #############
            decay_dic = {}
            for dec in decays[ipart].decays:
                psusy = 0
                newname = oldname
                dlist = [["Q", msq]]
                for dauid in sorted(dec.ids):
                    if not is_SUSY(dauid):
                        newname += get_name(dauid)
                        dlist.append([get_name(dauid), get_mass(blocks, get_name(dauid)), dauid])
                    else:
                        psusy = abs(dauid)
                if not psusy == 0:
                    newname += get_name(psusy)
                    dlist.append([get_name(psusy), get_mass(blocks, get_name(psusy)), psusy])
                oldsig = 0
                if newname in Sq.chain: oldsig = Sq.chain[newname].Br
                newsig = oldsig + sigtot * dec.br
                newData = Structure()
                newData.lastid = psusy
                newData.Br = newsig
                newData.plist = Sq_Data.plist + get_plist(dlist)
                Sq.chain[newname] = newData
            #############
            sq_list.append(Sq)

        #==============================================
        #  Now we have a list of light flavour squarks
        #-------------------------------------------
        #for sq_dm in part_sq:
        #    for dm_mode, dm_sig in sq_dm.chain.iteritems(): print dm_mode, dm_sig
        #-------------------------------------------
        #  Add the squarks and normalise
        #==============================================
        Squark = reduce(add_particles, sq_list)
        #Squark.output()
        #print Squark.get_sig()
        Squark.scale_to_one()
        #Squark.output()
        #print Squark.get_sig()

        #######################################
        ##  get the particle list
        #######################################
        #---  for the light flovour squarks
        part_dict["Q"] = Squark

    #"---  Now all the other particles"

    for npart in initial_part_list:
        ipart = get_pid(npart)
        part_dict[npart] = TheParticle(ipart, lsp_id, mass_dic[abs(ipart)], cf)
    #for key, pdm in part_dict.iteritems(): print key, pdm.chain
    #print ""
    #########################################
    ##   get all processes
    #########################################
    dummy_list = part_dict.copy()
    for npart, Part in dummy_list.iteritems():
        #print npart, Part.chain.keys()

        #======================================
        while Part.can_decay():
            Part.get_status()
            #print ""
            #inpdm = raw_input(' press return ')
            #print "---    ", Part.status, "   ---"
            Part.get_status()
            dummy_dict = Part.chain.copy()
            #=========================
            for proc, chain_Data in dummy_dict.iteritems():
                idecpart = chain_Data.lastid
                sigtot = chain_Data.Br
                old_plist = chain_Data.plist
                #print proc, idecpart, sigtot
                if idecpart == lsp_id or not is_SUSY(idecpart): continue
                Part.chain.pop(proc)
                oldname = proc
                #===========
                for dec in decays[idecpart].decays:
                    psusy = 0
                    br = dec.br
                    newname = oldname
                    dlist = [[get_name(idecpart), get_mass(blocks, get_name(idecpart)), idecpart]]
                    newname_list = []
                    for dauid in sorted(dec.ids):
                        if not is_SUSY(dauid):
                            newname_list.append(get_name(dauid))
                            dlist.append([get_name(dauid), get_mass(blocks, get_name(dauid)), dauid])
                        else:
                            psusy = abs(dauid)
                    newname_list.sort()
                    for dauname in newname_list: newname += dauname
                    if not psusy == 0:
                        newname += get_name(psusy)
                        dlist.append([get_name(psusy), get_mass(blocks, get_name(psusy)), psusy])
                    if br < 0:
                        mess = "Br<0 for " + str(idecpart) + "->" + str(dec.ids) \
                            + " " + str(dec.br) + " :Flip the sign and include"
                        err_dict[mess] = 0
                        br = abs(br)
                    oldsig = 0
                    if newname in Part.chain: oldsig = Part.chain[newname].Br
                    newsig = oldsig + sigtot * br
                    newData = Structure()
                    #print newname, psusy
                    newData.lastid = psusy
                    newData.Br = newsig
                    newData.plist = old_plist + get_plist(dlist)
                    if newsig < 10**-6 or isnan(newsig):
                        #print 'skip', newname, newsig
                        continue
                    Part.chain[newname] = newData
                #===========
            #===========================
        #=====================================================
        #Part.output()
        #print Part.get_sig()
        stot = Part.get_sig()
        if abs(1 - stot) > 0.1:
           print "!!!  ERROR  !!!"
           print "[Br(" + npart + "->all) = ", stot, "]:  !!! The total Br is deviated from 1 !!!"
           print ""
           exit
        Part.scale_to_one()
        part_dict[npart] = Part
    #=====================================================
    return part_dict, err_dict


def get_plist(dlist):
    plist = []
    for i in range(1,len(dlist)):
        pname = dlist[i][0]
        pData = Structure()
        pData.pname = pname
        pData.mass = dlist[i][1]
        pData.pid = dlist[i][2]
        nbody = len(dlist) - 1
        pData.p, pData.v = get_p_and_v(dlist, i)
        plist.append(pData)
    return plist

def get_p_and_v(dlist, i):
    if len(dlist) == 3:  # 2-body decay
        m = [dlist[0][1], dlist[1][1], dlist[2][1]]
        p = get_p_2body(m)
        e = sqrt(p**2 + m[i]**2)
        v = 0
        if e > 0: v = p/e
        return p, v
    elif len(dlist) == 4:  # 3-body decay
        m = [dlist[0][1], dlist[1][1], dlist[2][1], dlist[3][1]]
        p = get_p_3body(m)
        e = sqrt(p**2 + m[i]**2)
        v = 0
        if e > 0: v = p/e
        return p, v
    else:
        m = [dlist[0][1], dlist[1][1], dlist[2][1], dlist[3][1], dlist[4][1]]
        p = get_p_4body(m)
        e = sqrt(p**2 + m[i]**2)
        v = 0
        if e > 0: v = p/e
        return p, v

def get_p_2body(mass):
    m0 = mass[0]
    m1 = mass[1]
    m2 = mass[2]
    if m0 == m2: return 0
    aa = (m0**2 - (m1 + m2)**2) * (m0**2 - (m1 - m2)**2)
    try:
        p = sqrt( aa ) / (2 * m0)
    except:
        p = 0
    return p

def get_p_3body(mass):
    if (mass[1] - mass[2])/mass[0] < 0.0001:
        m0 = mass[0]
        m1 = mass[1]
        m2 = mass[3]
        if m0 == m2: return 0
        aa = m0**2 * (m0**2 + 3*m1**2 - 3*m2**2)
        try:
            bb = 5*m0**2 + 3*m1**2 - 12*m2**2 - 4*sqrt(aa)
        except:
            bb = -1.
        try:
            p = (1./3) * sqrt( bb )
        except:
            p = (mass[0] - mass[1] - mass[2] - mass[3]) / 3
        return p
    elif min(mass[1], mass[2])/mass[0] < 0.0001:
        m0 = mass[0]
        m1 = max(mass[1], mass[2])
        m2 = mass[3]
        r2 = sqrt(2)
        r3 = sqrt(3)
        p = (m0 - m1 - m2) / 3
        return p
    else:
        p = (mass[0] - mass[1] - mass[2] - mass[3]) / 3
        return p

def get_p_4body(mass):
    p = (mass[0] - mass[1] - mass[2] - mass[3] - mass[4]) / 4
    return p



#########################################################

class TheParticle(object):

    def __init__(self, ipd, lsp, mass, charge_flag):
        self.status = "can decay"
        if charge_flag:
            self.name = get_name_charge(ipd)
        else:
            self.name = get_name_simple(ipd)
        self.ipd = ipd
        Data = Structure()
        Data.Br = 1.
        Data.lastid = ipd
        pData = Structure()
        pData.pname = self.name
        pData.pid = ipd
        pData.mass = mass
        pData.p = 0.
        pData.v = 0.
        Data.plist = [pData]
        #self.chain = {get_name(ipd):[abs(ipd), 1]}
        self.chain = {self.name : Data}
        self.lsp_id = lsp
        self.mass = mass

    def get_sig(self, lowlim = -10000):
        sigtot = 0
        for mode, chain_Data in self.chain.iteritems():
          if chain_Data.Br > lowlim:
             sigtot += chain_Data.Br
        return sigtot

    def scale(self, fac):
        for key in self.chain:
          self.chain[key].Br *= fac

    def scale_to_one(self):
        sigtot = self.get_sig()
        self.scale(1/sigtot)

    def get_status(self):
        self.status = "done"
        for proc, Data in self.chain.iteritems():
            spart = Data.lastid
            #print proc, spart, self.lsp_id
            if is_SUSY(spart) and not spart == self.lsp_id:
                self.status = "can decay"

    def decayed(self):
        self.get_status()
        if self.status == "done":
            return True
        else:
            return False

    def can_decay(self):
        self.get_status()
        if self.status == "done":
            return False
        else:
            return True


#########################################################
def add_particles(Par1, Par2):
    Added = Par1
    for name2, Data2 in Par2.chain.iteritems():
        if name2 in Added.chain:
            Added.chain[name2].Br += Data2.Br
        else:
            Added.chain[name2] = Data2
    return Added



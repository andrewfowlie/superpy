#! /usr/bin/env python

"""
"""

__author__ = "M. Papucci, K. Sakurai, A. Weiler, L. Zeune"
__version__ = "1.0"

import os
import sys
from pyslha import * 
from math import *
from basic_func import *
sys.path.append('interpolation')
from interpolate import *
from interpolate import *
sys.path.append('statistics/CLs_limits')
from mycls import *
from extract_chain import *
from read_data import *
from display import *
from finalstate import *


###############################################################################

if __name__ == "__main__":

    try:
        inputdir, outputdir, pointlist = get_filelist(sys.argv)
        print "inputdir =", inputdir
        print "outputdir =", outputdir
        print "pointlist =", pointlist        
    except:
        print "Input Error!"
        print "Use correct command: (e.g.) ./ScanPoints.py slha_files/* ScanOutput"
        exit()

    print "Reading efficiency tables..."

    ana_list_8 = [
                "ATLAS_CONF_2013_024",    #0 lepton + 6 (2 b-)jets + Etmiss [Heavy stop] at 8TeV with $20.5fb^{-1}$
                "ATLAS_CONF_2013_035",    #3 leptons + Etmiss [EW production] at 8 TeV with $20.7fb^{-1}$.
                "ATLAS_CONF_2013_037",    #1 lepton + 4(1 b-)jets + Etmiss [Medium / heavy stop] at 8TeV with $20.7fb^{-1}$
                "ATLAS_CONF_2013_047",    #0 leptons + 2-6 jets + Etmiss [Incl. squarks & gluinos] at 8TeV with $20.3fb^{-1}$
                "ATLAS_CONF_2013_048",    #2 leptons (+ jets) + Etmiss [Medium stop] at 8 TeV with $20.3fb^{-1}$.
                "ATLAS_CONF_2013_049",    #2 leptons + Etmiss [EW production] at 8 TeV with $20.3fb^{-1}$.
                "ATLAS_CONF_2013_053",    #0 leptons + 2 b-jets + Etmiss [Sbottom/stop] at 8TeV with $20.1fb^{-1}$                
                "ATLAS_CONF_2013_054",    #0 leptons + >=7-10 jets + Etmiss [Incl. squarks & gluinos] at 8TeV with $20.3fb^{-1}$
                "ATLAS_CONF_2013_061",    #0-1 leptons + >=3 b-jets + Etmiss [3rd gen. squarks] at 8TeV with $20.1fb^{-1}$
                "ATLAS_CONF_2013_062",    #1-2 leptons + 3-6 jets + Etmiss [Incl. squarks & gluinos, mUED] at 8TeV with $20.3fb^{-1}$
                "ATLAS_CONF_2013_093"    #1-2 leptons + 3-6 jets + Etmiss [Incl. squarks & gluinos, mUED] at 8TeV with $20.3fb^{-1}$                
                ]

    ana_dict = get_ana_dict(ana_list_8, 8)

    proc_GN1_tab  = [
                    "GbbN1_GbbN1","GbbN1_GbtN1","GbbN1_GqqN1","GbbN1_GttN1","GbtN1_GbtN1",
                    "GbtN1_GqqN1","GbtN1_GttN1","GqqN1_GqqN1","GqqN1_GttN1","GttN1_GttN1",
                    "GgN1_GgN1","GgN1_GqqN1","GgN1_GttN1","GbbN1_GgN1","GbtN1_GgN1"
                    ]
    proc_T1N1_tab  = ["T1bN1_T1bN1","T1bN1_T1tN1","T1tN1_T1tN1"]
    proclist_2D = proc_GN1_tab + proc_T1N1_tab 

    proc_GT1N1_tab = ["GtT1bN1_GtT1bN1", "GtT1bN1_GtT1tN1", "GtT1tN1_GtT1tN1"]
    proc_GB1N1_tab = ["GbB1bN1_GbB1bN1", "GbB1bN1_GbB1tN1", "GbB1tN1_GbB1tN1"]
    proclist_3D = proc_GT1N1_tab + proc_GB1N1_tab

    eff_list_dict = {}
    eff_list_dict.update( get_eff_list2D(proclist_2D, ana_dict) )
    eff_list_dict.update( get_eff_list3D(proclist_3D, ana_dict) )


    ##############################################
    ###          Scanning model points         ### 
    ##############################################
    for point in pointlist:

        print "-----", point, "-----"
        outputlist = []

        infile = inputdir + '/' + point
        outfile = outputdir + '/' + point + '.FLout'

        warning_list = []
        blocks, decays = readSLHAFile(infile)

        Input = Paths_and_Data(infile, blocks, decays)

        mG = abs(get_mass(blocks,"G"))
        mC1 = abs(get_mass(blocks,"C1"))
        mN2 = abs(get_mass(blocks,"N2"))
        mN1 = abs(get_mass(blocks,"N1"))
        mT1 = abs(get_mass(blocks,"T1"))
        mB1 = abs(get_mass(blocks,"B1"))    
        mT2 = abs(get_mass(blocks,"T2"))
        mB2 = abs(get_mass(blocks,"B2"))    

        outputlist.append( ('mG', mG) )
        outputlist.append( ('mT1', mT1) )
        outputlist.append( ('mB1', mB1) )
        outputlist.append( ('mT2', mT2) )
        outputlist.append( ('mC1', mC1) )
        outputlist.append( ('mN2', mN2) )
        outputlist.append( ('mN1', mN1) )
        for out in outputlist: print out[0], out[1]            

        initial_part_list = []
        if mG < 2000: initial_part_list.append('G') 
        if mT1 < 1500: initial_part_list.append('T1') 
        if mB1 < 1500: initial_part_list.append('B1') 
        if mT2 < 1500: initial_part_list.append('T2') 
        if mB2 < 1500: initial_part_list.append('B2') 

        if len(initial_part_list) == 0:
            print 'no particles to be considered'
            exit()
        print 'Reading decays...'
        part_dict, err_dict = extract_chain(initial_part_list, Input, charge_flag = False)


        root_S = 8
        Prod_8 = TheProdMode()
        print 'Calculating sigma * BR...'    
        if 'G' in initial_part_list: Prod_8.add_prod(part_dict, "G", "G", get_xsec("G", "G", root_S, Input)) 
        if 'T1' in initial_part_list: Prod_8.add_prod(part_dict, "T1", "T1", get_xsec("T1", "T1", root_S, Input)) 
        if 'B1' in initial_part_list: Prod_8.add_prod(part_dict, "B1", "B1", get_xsec("B1", "B1", root_S, Input)) 
        if 'T2' in initial_part_list: Prod_8.add_prod(part_dict, "T2", "T2", get_xsec("T2", "T2", root_S, Input)) 
        if 'B2' in initial_part_list: Prod_8.add_prod(part_dict, "B2", "B2", get_xsec("B2", "B2", root_S, Input)) 

        procs_8 = Prod_8.processes  ### Do not comment out
        if 0 < mN2 - mC1 < 10:
            Replace(procs_8, "N2qqC1", "C1")
            Replace(procs_8, "N2enC1", "C1")
            Replace(procs_8, "N2mnC1", "C1")        
            Replace(procs_8, "N2ntaC1", "C1")        
        if 0 < mC1 - mN2 < 10:
            Replace(procs_8, "C1qqN2", "N2")
            Replace(procs_8, "C1enN2", "N2")
            Replace(procs_8, "C1mnN2", "N2")        
            Replace(procs_8, "C1ntaN2", "N2")        
        if 0 < mC1 - mN1 < 10:
            Replace(procs_8, "C1qqN1", "N1")
            Replace(procs_8, "C1enN1", "N1")
            Replace(procs_8, "C1mnN1", "N1")        
            Replace(procs_8, "C1ntaN1", "N1")        
        if 0 < mN2 - mN1 < 10:   
            Replace(procs_8, "N2qqN1", "N1")
            Replace(procs_8, "N2eeN1", "N1")
            Replace(procs_8, "N2mmN1", "N1")        
            Replace(procs_8, "N2tataN1", "N1")        
            Replace(procs_8, "N2nnN1", "N1")    
            Replace(procs_8, "N2gamN1", "N1")        
        set_rates(procs_8)          ### Do not comment out     
        #show_Prod(Prod_8, "C1m_C1p", upto = 0.0000000001)


        print 'Calculating results...'    

        procs_in_model = [proc for proc in procs_8]

        proc2D = {}      
        proc_GN1_all = proc_GN1_tab            
        proc_GN1 = list(set(procs_in_model) & set(proc_GN1_all))
        proc2D.update(get_proc2D_dict(proc_GN1, procs_8, mG, mN1))

        proc_T1N1_all = proc_T1N1_tab            
        proc_T1N1 = list(set(procs_in_model) & set(proc_T1N1_all))
        proc2D.update(get_proc2D_dict(proc_T1N1, procs_8, mT1, mN1))

        proc_T2N1_all  = ["T2bN1_T2bN1","T2bN1_T2tN1","T2tN1_T2tN1"]        
        proc_T2N1 = list(set(procs_in_model) & set(proc_T2N1_all))
        proc2D.update(get_proc2D_dict(proc_T2N1, procs_8, mT2, mN1))

        proc_B1N1_all  = ["B1bN1_B1bN1","B1bN1_B1tN1","B1tN1_B1tN1"]
        proc_B1N1 = list(set(procs_in_model) & set(proc_B1N1_all))
        proc2D.update(get_proc2D_dict(proc_B1N1, procs_8, mB1, mN1))

        proc_B2N1_all  = ["B2bN1_B2bN1","B2bN1_B2tN1","B2tN1_B2tN1"]
        proc_B2N1 = list(set(procs_in_model) & set(proc_B2N1_all))
        proc2D.update(get_proc2D_dict(proc_B2N1, procs_8, mB2, mN1))

        proc3D = {}
        proc_GT1N1_all = proc_GT1N1_tab
        proc_GT1N1 = list(set(procs_in_model) & set(proc_GT1N1_all))
        proc3D.update(get_proc3D_dict(proc_GT1N1, procs_8, mG, mT1, mN1))

        proc_GT2N1_all = ["GtT2bN1_GtT2bN1", "GtT2bN1_GtT2tN1", "GtT2tN1_GtT2tN1"]
        proc_GT2N1 = list(set(procs_in_model) & set(proc_GT2N1_all))
        proc3D.update(get_proc3D_dict(proc_GT2N1, procs_8, mG, mT2, mN1))

        if (mT2 - mT1)/mT2 < 10.: 
            proc_GT1_GT2_all = set(["GtT1bN1_GtT2bN1", "GtT1tN1_GtT2bN1", "GtT1bN1_GtT2tN1", "GtT1tN1_GtT2tN1"])
            proc_GT1_GT2 = list(set(procs_in_model) & set(proc_GT1_GT2_all)) 
            proc3D.update(get_proc3D_dict(proc_GT1_GT2, procs_8, mG, (mT2 + mT1)/2., mN1))
        else: proc3D_GT1_GT2 = []    

        proc_GB1N1_all = proc_GB1N1_tab
        proc_GB1N1 = list(set(procs_in_model) & set(proc_GB1N1_all))
        proc3D.update(get_proc3D_dict(proc_GB1N1, procs_8, mG, mB1, mN1))

        proc_GB2N1_all = ["GbB2bN1_GbB2bN1", "GbB2bN1_GbB2tN1", "GbB2tN1_GbB2tN1"]        
        proc_GB2N1 = list(set(procs_in_model) & set(proc_GB2N1_all))
        proc3D.update(get_proc3D_dict(proc_GT2N1, procs_8, mG, mB2, mN1))

        if (mB2 - mB1)/mB2 < 10.: 
            proc_GB1_GB2_all = set(["GtB1bN1_GtB2bN1", "GtB1tN1_GtB2bN1", "GtB1bN1_GtB2tN1", "GtB1tN1_GtB2tN1"])
            proc_GB1_GB2 = list(set(procs_in_model) & set(proc_GB1_GB2_all)) 
            proc3D.update(get_proc3D_dict(proc_GB1_GB2, procs_8, mG, (mB2 + mB1)/2., mN1))
        else: proc3D_GB1_GB2 = []    

        proc3D_list = proc_GT1N1 + proc_GB1N1 + proc_GT2N1 + proc_GB2N1 + proc_GT1_GT2 + proc_GB1_GB2
        proc2D_list = proc_GN1 + proc_T1N1 + proc_B1N1 + proc_T2N1 + proc_B2N1
        proc_list = proc3D_list + proc2D_list    

        results = {}
        res, warn = get_results_scan(ana_dict, proc3D, proc2D, eff_list_dict)
        results.update(res)

        ############################
        ###       Output         ### 
        ############################

        fout = open(outfile, "w")
        fout.write('\n')
        fout.close()

        show_outputlist(outputlist, outfile)
        #show_summary2(Prod_8, proc_list, results, __author__, __version__)        
        show_summary_out(Prod_8, proc_list, results, __author__, __version__, outfile)
        show_processes(Prod_8, proc_list, 8, outfile, upto = 0.005)
        show_analysis(ana_dict, results, outfile)


    exit()


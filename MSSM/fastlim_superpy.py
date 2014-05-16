#########################################################################
#                                                                       #
#    F a s t L i m - S u p e r P y                                      #
#                                                                       #
#########################################################################

# A python interface to FastLim. Provide a function that when imported,
# takes as an arguements an SLHA file with decay tables, and returns whether it
# it excluded.

import os
import sys

sys.path.append('../fastlim-1.0/')
sys.path.append('../fastlim-1.0/interpolation')
sys.path.append('../fastlim-1.0/statistics/CLs_limits')

from pyslha import *
from math import *
from basic_func import *
from interpolate import *
from mycls import *
from extract_chain import *
from read_data import *
from display import *
from finalstate import *
from update import *

__author__ = "M.Papucci, K.Sakurai, A.Weiler, L.Zeune"
__version__ = "1.0"

################################################################################

def excluded(infile):
    """" Find whether a spectrum in an SLHA file is excluded.
    Arguments:
    infile - name of SLHA file.

    Returns:
    Excluded - logical, whether excluded or not.
    Name -- name of first analysis that excluded it.
    """

    ############################################################################
    # New SuperPy code:
    # Check that infile exists.
    if not os.path.isfile(infile):
      print "Error: could not find SLHA file - excluded."
      return True

    ############################################################################
    # Original fast_lim.py code.

    logging.basicConfig(format='%(levelname)s:  %(message)s', level=logging.INFO)

    #try:
        #infile = sys.argv[1]
    #except:
        #logging.error("Invalid input!!")
        #logging.info("Use correct command: (e.g.) ./fastlim.py slha_files/testspectrum.slha")
        #exit()

    options = {'fastlimdir' : os.path.join(os.path.dirname(os.path.realpath(__file__)),"fastlim-1.0"), \
               'version_major' : ((__version__).split("."))[0],
               'version_minor' : ((__version__).split("."))[1],
               'authors' : __author__,
               'update_interval' : 7} # checking updates every 7 days

    print options["fastlimdir"]
    print_logo(options)

    #check_install_update(options)

    warning_list = []
    blocks, decays = readSLHAFile(infile)

    Input = Paths_and_Data(infile, blocks, decays)


    # Here the masses are stored. For example, get_mass(blocks,"T1") gives the T1 mass.
    mG = abs(get_mass(blocks,"G"))
    mC1 = abs(get_mass(blocks,"C1"))
    mN2 = abs(get_mass(blocks,"N2"))
    mN1 = abs(get_mass(blocks,"N1"))
    mT1 = abs(get_mass(blocks,"T1"))
    mB1 = abs(get_mass(blocks,"B1"))
    mT2 = abs(get_mass(blocks,"T2"))
    mB2 = abs(get_mass(blocks,"B2"))

    logging.info('mG  = %s' % mG)
    logging.info('mT1 = %s' % mT1)
    logging.info('mB1 = %s' % mB1)
    logging.info('mT2 = %s' % mT2)
    logging.info('mB2 = %s' % mB2)
    logging.info('mC1 = %s' % mC1)
    logging.info('mN2 = %s' % mN2)
    logging.info('mN1 = %s' % mN1)

    # Making a list of initial particles for which the branching ratios for all the decay chains are calculated.
    initial_part_list = []
    if mG < 2000: initial_part_list.append('G')
    if mT1 < 1500: initial_part_list.append('T1')
    if mB1 < 1500: initial_part_list.append('B1')
    if mT2 < 1500: initial_part_list.append('T2')
    if mB2 < 1500: initial_part_list.append('B2')

    if len(initial_part_list) == 0:
        logging.error('no particles to be considered')
        # exit()
        return False, ""

    logging.info('Reading decays...')
    # The branching ratops are calculated here.
    part_dict, err_dict = extract_chain(initial_part_list, Input, charge_flag = False)
    output_warning(err_dict)

    # To find the Brs, use the following function.
    #show_chain(part_dict, "T1", upto = 0.001)

    # Here, the cross section of all the topologies are calculated combining the produciton cross section and the branching ratios.  The 8TeV T1T1 production cross section can be obtained by get_xsec("T1", "T1", 8, Input)
    root_S = 8
    Prod_8 = TheProdMode()
    logging.info('Calculating sigma * BR...')
    if 'G' in initial_part_list: Prod_8.add_prod(part_dict, "G", "G", get_xsec("G", "G", root_S, Input))
    if 'T1' in initial_part_list: Prod_8.add_prod(part_dict, "T1", "T1", get_xsec("T1", "T1", root_S, Input))
    if 'B1' in initial_part_list: Prod_8.add_prod(part_dict, "B1", "B1", get_xsec("B1", "B1", root_S, Input))
    if 'T2' in initial_part_list: Prod_8.add_prod(part_dict, "T2", "T2", get_xsec("T2", "T2", root_S, Input))
    if 'B2' in initial_part_list: Prod_8.add_prod(part_dict, "B2", "B2", get_xsec("B2", "B2", root_S, Input))

    # The topologies are renamed if there is mass degeneracies satisfying the conditions below.
    procs_8 = Prod_8.processes  ### Do not comment out
    if 0 < mN2 - mC1 < 10:
        Replace(procs_8, "N2qqC1", "C1")
        Replace(procs_8, "N2enC1", "C1")
        Replace(procs_8, "N2mnC1", "C1")
        Replace(procs_8, "N2e3nC1", "C1")
    if 0 < mC1 - mN2 < 10:
        Replace(procs_8, "C1qqN2", "N2")
        Replace(procs_8, "C1enN2", "N2")
        Replace(procs_8, "C1mnN2", "N2")
        Replace(procs_8, "C1e3nN2", "N2")
    if 0 < mC1 - mN1 < 10:
        Replace(procs_8, "C1qqN1", "N1")
        Replace(procs_8, "C1enN1", "N1")
        Replace(procs_8, "C1mnN1", "N1")
        Replace(procs_8, "C1e3nN1", "N1")
    if 0 < mN2 - mN1 < 10:
        Replace(procs_8, "N2qqN1", "N1")
        Replace(procs_8, "N2eeN1", "N1")
        Replace(procs_8, "N2mmN1", "N1")
        Replace(procs_8, "N2e3e3N1", "N1")
        Replace(procs_8, "N2nnN1", "N1")
        Replace(procs_8, "N2gamN1", "N1")
    set_rates(procs_8)          ### Do not comment out

    # List of all the topologies with the T1_T1 production is displayed by the following funciton.
    #show_Prod(Prod_8, "T1_T1", upto = 0.0000000001)

    ################################################

    procs_in_model = set([proc for proc in procs_8])

    # Here, efficiency tables and the corresponding masses are set.
    proc3D = {}
    proc3D_GT1N1_all = set(["GtT1bN1_GtT1bN1", "GtT1bN1_GtT1tN1", "GtT1tN1_GtT1tN1"])
    proc3D_GT1N1 = list(procs_in_model & proc3D_GT1N1_all)
    proc3D.update(get_proc3D_dict(proc3D_GT1N1, procs_8, mG, mT1, mN1))

    proc3D_GT2N1_all = set(["GtT2bN1_GtT2bN1", "GtT2bN1_GtT2tN1", "GtT2tN1_GtT2tN1"])
    proc3D_GT2N1 = list(procs_in_model & proc3D_GT2N1_all)
    proc3D.update(get_proc3D_dict(proc3D_GT2N1, procs_8, mG, mT2, mN1))

    if (mT2 - mT1)/mT2 < 0.1:
        proc3D_GT1_GT2_all = set(["GtT1bN1_GtT2bN1", "GtT1tN1_GtT2bN1", "GtT1bN1_GtT2tN1", "GtT1tN1_GtT2tN1"])
        proc3D_GT1_GT2 = list(procs_in_model & proc3D_GT1_GT2_all)
        proc3D.update(get_proc3D_dict(proc3D_GT1_GT2, procs_8, mG, (mT2 + mT1)/2., mN1))
    else: proc3D_GT1_GT2 = []

    proc3D_GB1N1_all = set(["GbB1bN1_GbB1bN1", "GbB1bN1_GbB1tN1", "GbB1tN1_GbB1tN1"])
    proc3D_GB1N1 = list(procs_in_model & proc3D_GB1N1_all)
    proc3D.update(get_proc3D_dict(proc3D_GB1N1, procs_8, mG, mB1, mN1))

    proc3D_GB2N1_all = set(["GbB2bN1_GbB2bN1", "GbB2bN1_GbB2tN1", "GbB2tN1_GbB2tN1"])
    proc3D_GB2N1 = list(procs_in_model & proc3D_GB2N1_all)
    proc3D.update(get_proc3D_dict(proc3D_GB2N1, procs_8, mG, mB2, mN1))

    if (mB2 - mB1)/mB2 < 0.1:
        proc3D_GB1_GB2_all = set(["GtB1bN1_GtB2bN1", "GtB1tN1_GtB2bN1", "GtB1bN1_GtB2tN1", "GtB1tN1_GtB2tN1"])
        proc3D_GB1_GB2 = list(procs_in_model & proc3D_GB1_GB2_all)
        proc3D.update(get_proc3D_dict(proc3D_GB1_GB2, procs_8, mG, (mB2 + mB1)/2., mN1))
    else: proc3D_GB1_GB2 = []

    proc3D_list = proc3D_GT1N1 + proc3D_GB1N1 + proc3D_GT2N1 + proc3D_GB2N1 + proc3D_GT1_GT2 + proc3D_GB1_GB2


    proc2D = {}
    proc_GN1_all  = set([
                    "GbbN1_GbbN1","GbbN1_GbtN1","GbbN1_GqqN1","GbbN1_GttN1","GbtN1_GbtN1",
                    "GbtN1_GqqN1","GbtN1_GttN1","GqqN1_GqqN1","GqqN1_GttN1","GttN1_GttN1",
                    "GgN1_GgN1","GgN1_GqqN1","GgN1_GttN1","GbbN1_GgN1","GbtN1_GgN1"
                    ])
    proc_GN1 = list(procs_in_model & proc_GN1_all)
    proc2D.update(get_proc2D_dict(proc_GN1, procs_8, mG, mN1))

    proc_T1N1_all  = set(["T1bN1_T1bN1","T1bN1_T1tN1","T1tN1_T1tN1"])
    proc_T1N1 = list(procs_in_model & proc_T1N1_all)
    proc2D.update(get_proc2D_dict(proc_T1N1, procs_8, mT1, mN1))

    proc_T2N1_all  = set(["T2bN1_T2bN1","T2bN1_T2tN1","T2tN1_T2tN1"])
    proc_T2N1 = list(procs_in_model & proc_T2N1_all)
    proc2D.update(get_proc2D_dict(proc_T2N1, procs_8, mT2, mN1))

    proc_B1N1_all  = set(["B1bN1_B1bN1","B1bN1_B1tN1","B1tN1_B1tN1"])
    proc_B1N1 = list(procs_in_model & proc_B1N1_all)
    proc2D.update(get_proc2D_dict(proc_B1N1, procs_8, mB1, mN1))

    proc_B2N1_all  = set(["B2bN1_B2bN1","B2bN1_B2tN1","B2tN1_B2tN1"])
    proc_B2N1 = list(procs_in_model & proc_B2N1_all)
    proc2D.update(get_proc2D_dict(proc_B2N1, procs_8, mB2, mN1))

    proc2D_list = proc_GN1 + proc_T1N1 + proc_B1N1 + proc_T2N1 + proc_B2N1


    # Defining a list of the analyses that you want to consider
    ana_list_8 = os.listdir(os.path.join(options['fastlimdir'],'analyses_info','8TeV'))
    # ana_list_8 = [
    #             "ATLAS_CONF_2013_024",    #0 lepton + 6 (2 b-)jets + Etmiss [Heavy stop] at 8TeV with $20.5fb^{-1}$
    #             "ATLAS_CONF_2013_035",    #3 leptons + Etmiss [EW production] at 8 TeV with $20.7fb^{-1}$.
    #             "ATLAS_CONF_2013_037",    #1 lepton + 4(1 b-)jets + Etmiss [Medium / heavy stop] at 8TeV with $20.7fb^{-1}$
    #             "ATLAS_CONF_2013_047",    #0 leptons + 2-6 jets + Etmiss [squarks & gluinos] at 8TeV with $20.3fb^{-1}$
    #             "ATLAS_CONF_2013_048",    #2 leptons (+ jets) + Etmiss [Medium stop] at 8 TeV with $20.3fb^{-1}$.
    #             "ATLAS_CONF_2013_049",    #2 leptons + Etmiss [EW production] at 8 TeV with $20.3fb^{-1}$.
    #             "ATLAS_CONF_2013_053",    #0 leptons + 2 b-jets + Etmiss [Sbottom/stop] at 8TeV with $20.1fb^{-1}$
    #             "ATLAS_CONF_2013_054",    #0 leptons + >=7-10 jets + Etmiss [squarks & gluinos] at 8TeV with $20.3fb^{-1}$
    #             "ATLAS_CONF_2013_061",    #0-1 leptons + >=3 b-jets + Etmiss [3rd gen. squarks] at 8TeV with $20.1fb^{-1}$
    #             "ATLAS_CONF_2013_062",    #1-2 leptons + 3-6 jets + Etmiss [squarks & gluinos, mUED] at 8TeV with $20.3fb^{-1}$
    #             "ATLAS_CONF_2013_093"     #1 lepton + bb(H) + Etmiss [EW production] at 8TeV with $20.3fb^{-1}$
    #             ]

    # Information about the analyses is obtained from the files in analyses_info.
    ana_dict_8 = get_ana_dict(ana_list_8, 8)
    ana_dict = {}
    ana_dict.update(ana_dict_8)

    ###############################################################

    logging.info('Reading efficiency tables...')
    # Here calculates the visible cross sections and R measures. These are stored in results dictionary
    results = {}
    res, warn = get_results(ana_dict_8, proc3D, proc2D)
    results.update(res)
    warning_list += warn
    for w in warning_list: print w

    # The efficiency in T1tN1_T1tN1 topology for ATLAS_CONF_2013_047, cut_0 is obtained by the following command.
    #print "eff(ATLAS_CONF_2013_047, cut_0) =", results["ATLAS_CONF_2013_047", "cut_0"].SR_info
    #print "eff(ATLAS_CONF_2013_047, cut_0) =", results["ATLAS_CONF_2013_047", "cut_0"].proc_Data["T1tN1_T1tN1"].eff

    ###############################################################

    proc_list = proc3D_list + proc2D_list

    outputfile="fastlim.out"
    fout = open(outputfile, "w")
    fout.write('\n')
    fout.close()

    # Creating the output file
    show_processes(Prod_8, proc_list, 8, outputfile, upto = 0.005)
    show_analysis(ana_dict, results, outputfile)

    ####################################################################

    # Showing display output
    show_summary(Prod_8, proc_list, results, __author__, __version__)

    ############################################################################

    # End of original fast_lim.py code.
    # New SuperPy code:
    for key in sorted(results.keys()):
      if results[key].Rvis_obs > 1:
        print "Point is excluded."
        return True, key[0]
    return False, ""




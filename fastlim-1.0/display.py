#! /usr/bin/env python
import sys
from basic_func import *
sys.path.append('statistics/CLs_limits')
from mycls import *
from read_data import *
from finalstate import *
#from PoissonModelWithBackg import *


def print_logo(options):
    print " "
    print "#" * 65
    print "#", " " * 61, "#"
    print "#", " " * 61, "#"
    print "#", " " *17, "This is Fastlim  ver. %(version_major)s.%(version_minor)s"  % options, " "* 17,    "#"
    print "#", " " * 61, "#"
    print "#", " " *8, "A fast LHC limit calculator for BSM models.", " "*8,"#"
    print "#", " " * 61, "#"
    print "#"," "*18,"http://cern.ch/fastlim/"," "*18,"#"
    print "#", " " * 61, "#"
    print "#    for bugs/feature requests: fastlim.developers@gmail.com    #"
    print "#", " " * 61, "#"
    print "#    please refer to:", " " * 41, "#"
    print "#    %(authors)s arXiv:1402.0492    #" % options
    print "#", " " * 61, "#"
    print "#" * 65
    print " "

##########################################################################
def show_relevantSRs(procs, upto = 10, upto_proc = 0.01):
    print "#" * 60
    print " " * 12, "Summary of Relevant Signal Regions"
    print "#" * 60
    SR_dict = get_SRs_for_model(procs)
    SR_short = {} 
    i = 0
    for SR, srData in sorted(SR_dict.items(), key=lambda x:x[1].xsec, reverse=True):
        i += 1
        SR_short[SR] = srData
        if i > upto: break
    for SR, srData in sorted(SR_short.items(), key=lambda x:x[1].xsec, reverse=True):
        print "----  [" + str(SR) + "]  ----"
        print "Cross Section: ", str(round(SR_short[SR].xsec, 5)), "fb"
        print "Typical MET: ", str(round(SR_short[SR].MET, 1))
        sp_prc = get_max_key_len(SR_short[SR].proc_dict)
        print "Process".rjust(sp_prc), "Xsec/fb".rjust(9), "<MET>".rjust(8)
        for prc, prcData in sorted(SR_short[SR].proc_dict.items(), key=lambda x:x[1].xsec, reverse=True):
            print prc.rjust(sp_prc), str(round(prcData.xsec, 5)).rjust(9), str(round(prcData.MET, 1)).rjust(8)        
            if prcData.xsec / SR_short[SR].xsec < upto_proc: break

##########################################################################
def show_relevantSRs_short(procs, upto = 10):
    print "----  Output relevant Signal Regions  ----"
    SR_dict = get_SRs_for_model(procs)
    SR_short = {}
    i = 0
    for SR, srData in sorted(SR_dict.items(), key=lambda x:x[1].xsec, reverse=True):
        i += 1
        SR_short[SR] = srData
        if i > upto: break
    space = get_max_key_len(SR_short)
    print "Signal Region".rjust(space), "Xsec/fb".rjust(9), "<MET>".rjust(8)
    for SR, srData in sorted(SR_short.items(), key=lambda x:x[1].xsec, reverse=True):
        print SR.rjust(space), str(round(srData.xsec, 5)).rjust(9), str(round(srData.MET, 1)).rjust(8)

##################################
def show_SR(procs, SR, upto = 0.01):
    SR_dict = get_SRs_for_model(procs)
    if not SR in SR_dict:
        print SR, "is not in the list"
        return
    print "----  [" + str(SR) + "]  ----"
    print "Cross Section: ", str(round(SR_dict[SR].xsec, 5)), "fb"
    print "Typical MET: ", str(round(SR_dict[SR].MET, 1))
    sp_prc = get_max_key_len(SR_dict[SR].proc_dict)
    print "Process".rjust(sp_prc), "Xsec/fb".rjust(9), "<MET>".rjust(8)
    for prc, prcData in sorted(SR_dict[SR].proc_dict.items(), key=lambda x:x[1].xsec, reverse=True):
        print prc.rjust(sp_prc), str(round(prcData.xsec, 5)).rjust(9), str(round(prcData.MET, 1)).rjust(8)        
        if prcData.xsec / SR_dict[SR].xsec < upto: break

##########################################################################
def show_finalstate(procs, prc):
    if not prc in procs:
        print prc, "is not in the processes"
        return
    print "-" * 4, " " + str(prc) + ": Final state information ", "-" * 4
    print "Cross Section:", str(round(procs[prc].xsec, 5)) + " fb"
    print "particle".rjust(8), "mass".rjust(9), "<momentum>".rjust(11), "<velocity>".rjust(11)
    for pData in procs[prc].plist:        
        print str(pData.pname).rjust(8), str(round(pData.mass, 2)).rjust(9), 
        print str(round(pData.p, 2)).rjust(11), str(round(pData.v, 2)).rjust(11)
    print "[Relevant signal regions]"
    relevant_SR_dict = get_SRs_for_proc(procs, prc)
    space = get_max_key_len(relevant_SR_dict) 
    print "Final state".rjust(space), "Xsec/fb".rjust(9), "<MET>".rjust(8)
    for sr, srData in sorted(relevant_SR_dict.items(), key=lambda x:x[1].xsec, reverse=True):
        print sr.rjust(space), str(round(srData.xsec, 5)).rjust(9), str(round(srData.MET, 1)).rjust(8)


##########################################################################
def show_chain(part_dict, pname, upto = -100):
    print "-" * 60 
    if not pname in part_dict:
        print pname, "is not found in part_dict"
        print "Choose from", part_dict.keys()
        return
    if upto > 0: 
        print "[" + str(pname) + "] decays with Br > " + str(upto * 100) + "%"
    else:
        print "[" + str(pname) + "] decays"
    chain_short = {}
    for chain, chain_Data in sorted(part_dict[pname].chain.items(), key=lambda x:x[1].Br, reverse=True):
        chain_short[chain] = chain_Data
        if chain_Data.Br < upto: break
    space = get_max_key_len(chain_short) + 1
    print "Chain".rjust(space), ":", "Brancing ratio"
    for chain, short_Data in sorted(chain_short.items(), key=lambda x:x[1].Br, reverse=True):
        print chain.rjust(space), ":",
        print ('%.2f' % float(100 * short_Data.Br)) + "%"  
    return


##########################################################################
def show_Prod_modes(Prod, upto, prod):
    print "-" * 20, " ", prod, " ", "-" * 20 
    print "The production xsec (" + str(prod) + ") =", 
    print ('%.2f' % float( Prod.prod_xsec[prod] ) + " [fb]")
    if upto > 0:
        print "Output upto " + str(upto * 100) + "%"
    else:
        print "Output all"
    modes_short = {}
    for proc, modes_Data in sorted(Prod.modes[prod].items(), key=lambda x:x[1].rate, reverse=True):
        modes_short[proc] = modes_Data
        if modes_Data.rate < upto: break
    space = get_max_key_len(modes_short)
    print "Modes".rjust(space),":","xsec/fb   (rate)"
    for proc, short_Data in sorted(modes_short.items(), key=lambda x:x[1].rate, reverse=True):
        print proc.rjust(space),":",
        print ('%.7f' % float( short_Data.xsec) ),
        print ('(%.2f' % float(100 * short_Data.rate)) + "%)"  
    return

def show_Prod_procs(Prod, upto):
    print "-" * 15, " Output inclusive processes ", "-" * 15 
    print "The total xsec =", 
    print ('%.2f' % float( Prod.xsec_tot ) + " [fb]")
    if upto > 0:
        print "Output upto " + str(upto * 100) + "%"
    else:
        print "Output all"
    proc_short = {}
    for proc, Data in sorted(Prod.processes.items(), key=lambda x:x[1].rate, reverse=True):
        proc_short[proc] = Data
        if Data.rate < upto: break
    space = get_max_key_len(proc_short)
    print "Modes".rjust(space),":","xsec/fb   (rate)"
    for proc, short_Data in sorted(proc_short.items(), key=lambda x:x[1].rate, reverse=True):
        print proc.rjust(space),":",
        print ('%.7f' % float( short_Data.xsec) ),
        print ('(%.2f' % float(100 * short_Data.rate)) + "%)"  
    return

def show_Prod(Prod, _prod_, upto = -100):	
	if _prod_ in Prod.modes:
		show_Prod_modes(Prod, upto, _prod_)
		return
	elif _prod_ == "inc":
		show_Prod_procs(Prod, upto)
		return
	else:
		print "production", _prod_, "not found"
        print "Choose from:",
        for prod_mode in Prod.modes: print "[" + str(prod_mode) + "]",
        print ""
        return

#################################################################################



def show_summary(Prod8, model_list, results, __author__, __version__):
    # print "_" * 90
    # print " "
    # print "  " * 18, "FastLim " + __version__
    # print " "
    # print "  " * 10, "(" + __author__ + ")"
    # print " "
    print "_" * 90
    print " "
    xtot_8 = Prod8.xsec_tot
    xsec_imp8 = get_implemented_xsec(Prod8, model_list)
    print "----------    Cross Section    ----------"
    print "Ecm ".rjust(4), "Total".rjust(11), "  Implemented".rjust(12), "  Coverage".rjust(01)
    print "8TeV".rjust(4), str(round(xtot_8, 3)).rjust(9) + "fb", str(round(xsec_imp8, 3)).rjust(11) + "fb", str(round(100. * xsec_imp8/xtot_8, 2)).rjust(9) + "%" 
    print " "
    #print ""    
    print "-" * 90
    print "Analysis".rjust(19), "E/TeV".rjust(6), "L*fb".rjust(5), "Signal Region:".rjust(25), "Nev/N_UL".rjust(9), "CLs".rjust(7)
    print "-" * 90
    i_line = 0
    for key in sorted(results.keys()):
        vals = results[key]
        i_line += 1        
        cls = "--"
        if 0.9 < float(vals.Rvis_obs) < 1.1: 
            try:
                cls = get_cls(vals)
            except:
                pass        
        print vals.analysis.rjust(19), str(vals.root_S).rjust(6), str(vals.lumi).rjust(5), (vals.SR_info + ':').rjust(25),
        #print ":",
        print ( '%.4f' % float(vals.Rvis_obs) ).rjust(9),
        try:
            print ( '%.4f' % float(cls) ).rjust(7),
        except:
            print cls.rjust(7),        
        excl = ""
        if vals.Rvis_obs > 1: excl = "  <== Exclude"
        print excl
    print "-" * 90
    print " "
    print " "    




def show_summary_out(Prod8, model_list, results, __author__, __version__, outputfile = 'fastlim.out'):

    fout = open(outputfile, "a")

    lines = []
    lines.append( "_" * 90 )
    lines.append( " " )
    lines.append( "  " * 18 + "FastLim " + __version__ )
    lines.append( " " )
    lines.append( "  " * 10 + "(" + __author__ + ")" )
    lines.append( " " )
    lines.append( "_" * 90 )
    lines.append( " " )
    xtot_8 = Prod8.xsec_tot
    xsec_imp8 = get_implemented_xsec(Prod8, model_list)
    lines.append( "-----    Cross Section    -----" )
    lines.append( " Ecm ".rjust(5) + "Total".rjust(12) + "Implemented".rjust(12) )
    lines.append( "8 TeV".rjust(5) + str(round(xtot_8, 4)).rjust(9) + " fb" + str(round(xsec_imp8, 4)).rjust(9) + " fb" )
    lines.append( " " )
    lines.append( "-" * 90 )
    lines.append( "Analysis".rjust(19) + "Ecm".rjust(5) + "L/fb".rjust(6) + "Signal Region".rjust(26) + ":" +
                 "Nev/Nev_UL".rjust(11) + "CLs".rjust(8) )
    lines.append( "-" * 90 )
    i_line = 0
    for key in sorted(results.keys()):
        vals = results[key]
        i_line += 1        
        cls = "--"
        if 0.9 < float(vals.Rvis_obs) < 1.1: 
            try:
                cls = get_cls(vals)
            except:
                pass        
        L0 = vals.analysis.rjust(19) + str(vals.root_S).rjust(5) + str(vals.lumi).rjust(6) + vals.SR_info.rjust(26) + ":" + str('%.4f' % float(vals.Rvis_obs) ).rjust(11)
        try:
            L0 += str( '%.4f' % float(cls) ).rjust(8)
        except:
            L0 += cls.rjust(8)        
        excl = ""
        if vals.Rvis_obs > 1: excl = "  <== Exclude"
        L0 += excl
        lines.append(L0)
    lines.append( "-" * 90 )
    lines.append( " " )
    lines.append( " " )    
    for line in lines: fout.write( line + '\n' )


def show_outputlist(outputlist, outputfile = "fastlim.out"):

    fout = open(outputfile, "a")
    fout.write("#" * 60+"\n")
    fout.write("                 Outputlist " +"\n")
    fout.write("#" * 60+"\n")
    for out in outputlist: fout.write( str(out[0]) +'  '+ str(out[1]) + '\n' )


def show_analysis(ana_dict, results, outputfile = "fastlim.out"):

    #outputfile="fastlim.out"
    fout = open(outputfile, "a")
    fout.write("#" * 60+"\n")
    fout.write("                 Analyses Details " +"\n")
    fout.write("#" * 60+"\n")
    ana_prev = ""
    sp = " " * 3
    for ana, aData in ana_dict.iteritems(): 
        fout.write("-" * 60+"\n")
        fout.write("[" + str(ana) + "]" +"\n")
        #fout.write("(" + str(aData.ana_info) + ")"+"\n")
        fout.write( str(aData.ana_info) +"\n" )
        fout.write("Ecm/TeV = "+ str(aData.root_S)+"\n")
        fout.write("lumi*fb = "+ str(aData.lumi) +"\n")
        for sr, srData in aData.sr_dict.iteritems():
	    	vals = results[ana, sr]
	        fout.write(sp+ "#----  " + str(vals.SR_info) + "  ----#" +"\n")
	        fout.write(sp+ "Nobs:              "+ str(vals.Nobs)  +"\n")
	        fout.write(sp+ "Nbg:               "+ str(vals.Nbg) + "(" + str(vals.errNbg) + ")"  +"\n")
	        if vals.UL_nvis_obs > 0: fout.write(sp+ "Nvis_UL[observed]: "+ str(vals.UL_nvis_obs)  +"\n")
	        if vals.UL_nvis_exp > 0: fout.write(sp+ "Nvis_UL[expected]: "+ str(vals.UL_nvis_exp)  +"\n")
	        #for prc, Data in vals.proc_Data.iteritems():
	        pspace = get_max_key_len(vals.proc_Data) 
	        fout.write(sp+ "Process".rjust(pspace)+"Nev".rjust(10)+"R[obs]".rjust(10))
	        if vals.Rvis_exp > 0: 
	            fout.write("R[exp]".rjust(10)  +"\n")
	        else:
	            fout.write(""  +"\n")

	        fout.write(sp+ "Total".rjust(pspace)) 
	        fout.write(( '%.4f' % float(vals.nvis) ).rjust(10))
	        try:
	            fout.write(( '%.4f' % float(vals.Rvis_obs) ).rjust(10))
	        except:
	            fout.write( "".rjust(10))
	        if vals.Rvis_exp > 0: 
	            fout.write(('%.4f' % float(vals.Rvis_exp) ).rjust(10))
	        else:
	            fout.write("")
	        note = ""
	        if vals.Rvis_obs > 1: note = " <== Exclude"        
	        fout.write(note+"\n")

	        for prc, Data in sorted(vals.proc_Data.items(), key=lambda x:x[1].xvis, reverse=True):
	            Robs = ""
	            Rexp = ""
	            if vals.UL_nvis_obs > 0: Robs = Data.nvis/vals.UL_nvis_obs
	            if vals.UL_nvis_exp > 0: Rexp = Data.nvis/vals.UL_nvis_exp
	            if vals.UL_xvis_obs > 0: Robs = Data.xvis/vals.UL_xvis_obs
	            if vals.UL_xvis_exp > 0: Rexp = Data.xvis/vals.UL_xvis_exp
	            fout.write(sp+ prc.rjust(pspace))
	            fout.write(( '%.4f' % float(Data.nvis) ).rjust(10))
	            try:
	                fout.write(( '%.4f' % float(Robs) ).rjust(10))
	            except:
	                fout.write( "".rjust(10))
	            try:
	                fout.write( ( '%.4f' % float(Rexp) ).rjust(10))
	            except:
	                fout.write( "".rjust(10))
		    fout.write("\n")
        fout.write("\n")	    
    fout.close()

def show_proc_dict(proc_dict, upto = 0.99):
    proc_short = {}
    Accum = 0
    for proc, Data in sorted(proc_dict.items(), key=lambda x:x[1].rate, reverse=True):
        Accum += Data.rate
        proc_short[proc] = Data
        if Accum > upto: break
    space = get_max_key_len(proc_short)    
    print "Modes".rjust(space) + " :  xsec/fb   (Rate)   (Accum)" 
    Accum = 0
    count = 0
    for proc, short_Data in sorted(proc_short.items(), key=lambda x:x[1].rate, reverse=True):
        Accum += short_Data.rate
        count += 1
        line = proc.ljust(space) 
        line += " : " + ('%.7f' % float( short_Data.xsec) ) 
        line += "  " + ('(%.2f' % float(100 * short_Data.rate)) + "%)" 
        line += "  " + ('(%.2f' % float(100 * Accum)) + "%)"
        line += "  " + str(count)
        print line 

def show_processes(Prod, model_list, root_S, outputfile = "fastlim.out", upto = 0.001):
    fout = open(outputfile, "a")
    xtot = Prod.xsec_tot 
    fout.write("#"*50 +"\n")
    fout.write("    Branching Ratio x Cross Section @ " + str(root_S) + " TeV  " +"\n")
    fout.write("#" * 50 +"\n")
    fout.write("-" * 50 +"\n")
    fout.write("Production:   Xsec/fb       Rate"  +"\n")  
    fout.write("     Total:".rjust(11))
    fout.write(( '%.3f' % float(xtot) ).rjust(10) + ( '%.2f' % float( 100. ) ).rjust(10) + "%" +"\n")
    for prod in Prod.prod_xsec: 
        fout.write((str(prod) + ":").rjust(11))
        fout.write( ( '%.3f' % float(Prod.prod_xsec[prod]) ).rjust(10) + ( '%.2f' % float( 100. * Prod.prod_xsec[prod]/xtot ) ).rjust(10) + "%" +"\n")
    fout.write("-" * 50 +"\n")
    Racum = 0 
    proc_show = {}
    for proc, pData in sorted(Prod.processes.items(), key=lambda x:x[1].rate, reverse=True):
        proc_show[proc] = pData
        if pData.rate < upto: break
    pspace = get_max_key_len(proc_show) + 1
    fout.write("Output processes upto " + str(upto * 100) + "%" + "\n")
    fout.write("Process".rjust(pspace)+ ":"+ "Br*Xsec/fb".rjust(12) + "Rate".rjust(9) + "Accum".rjust(9)  +"\n")
	#,"individual".rjust(12), "Accumulated".rjust(12)
    Racum = 0
    for proc, sData in sorted(proc_show.items(), key=lambda x:x[1].rate, reverse=True):
        note = ""
        if proc in model_list: note = "  <== Implemented"
        Racum += sData.rate
        fout.write(str(proc).rjust(pspace) + ":")
        fout.write(( '%.5f' % float(sData.xsec) ).rjust(12))          
        fout.write(( '%.2f' % float(sData.rate * 100) + "%").rjust(9)),          
        fout.write(( '%.2f' % float(Racum * 100) + "%").rjust(9)),
        fout.write(note +"\n")
    fout.write("..." +"\n")
    fout.write("\n")
    fout.write("\n")
    fout.close()

def get_implemented_xsec(Prod, model_list):
    xsec_imp = 0
    for proc in model_list:
        try:
            xsec_imp += Prod.processes[proc].xsec
        except:
            pass
    return xsec_imp


#! /usr/bin/env python

from math import *
from basic_func import *
#import numpy as np
import sys
try:
    import numpy as np
except Exception:
    print ''' WARNING:
    Could not load numpy. CLs value cannot be calculated!
    '''


def mycls(Nsig, errNsig, Nbkg, errNbkg, Nobs):
    distListSig = []
    distList = []
    nb = 0
    nsb = 0
    ntot = 100000
    bgdlist = np.random.normal(Nbkg,errNbkg,ntot)
    for i in range(0,ntot):
        distListSig.append(np.random.poisson(max(0.000001,Nsig+bgdlist[i])))
        distList.append(np.random.poisson(max(0.000001,bgdlist[i])))
        if distList[i]<Nobs:
            nb=nb+1
        if distListSig[i]<Nobs:
            nsb=nsb+1
    #print "nsb, nb", nsb, nb, Nsig, errNsig, Nbkg, errNbkg, Nobs
    if nsb>nb:
        cls=1
    else:
        cls=float(nsb)/float(nb)
    return cls


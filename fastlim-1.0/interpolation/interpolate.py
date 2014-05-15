
__author__ = "Andreas Weiler <andreas.weiler@desy.de>"
__version__ = "0.3"

from convexHull import *

_scipyloaded = False

try:
    import numpy as np
    from scipy.interpolate import LinearNDInterpolator
    _scipyloaded = True
except Exception:
    print ''' WARNING:
    Could not load numpy/scipy, please install! Will use cruder nearest-neighbor interpolation, might not always check if points are within interpolating
    boundaries.
    '''
    _scipyloaded = False

# see interpolate.nb for tests and caluclations

def Interpolate1D(listIn, x):
    '''linear interpolated function y = f(x)

    Expects a (not necessarily sorted) list (xi ,yi)

    listIn = [ [1.,2.] , [2,-3.], [4.,3.] ]

    and calculates linearly interpolated y at x.
    '''

    listIn = [[float(xi), float(yi)] for xi, yi in listIn]

    listIn.sort()   # sort in place
    minX = listIn[0][0]
    maxX = listIn[-1][0]    # last element
    if minX > maxX:
        errmsg = "No valid list given: " + listIn + " !"
        raise Exception(errmsg)
    elif not minX <= x <= maxX:
        errmsg = "x-value (" + str(x) + ") outside boundaries (" + \
        str(minX) + "," + str(maxX) + ")! \n"
        print errmsg
        #raise Exception(errmsg)
        return 'nan'

    distlist = [[abs(xel - x), listIn.index([xel, yel])] for xel, yel in listIn]
    distlist.sort()     # list of distance and index number
    x1Index = distlist[0][1]   # index of closest element
    x2Index = distlist[1][1]

    x1, y1 = listIn[x1Index]
    x2, y2 = listIn[x2Index]

    if x1 == x2:
        errmsg = "two entries for same x-value : " + str(x1)
        raise Exception(errmsg)

    yinterpolate = (y1 - y2) / (x1 - x2) * x - (x2 * y1 - x1 * y2) / (x1 - x2)

    return yinterpolate

def Interpolate2D(listIn, x):

    listCoordinates = [[float(xi), float(yi)] for xi, yi, zi in listIn]
    listIn = [[float(xi), float(yi), float(zi)] for xi, yi, zi in listIn]
    listValues = [ float(zi) for xi, yi, zi in listIn]

    p = (float(x[0]), float(x[1]))

    if _scipyloaded == False:       # do nearest neighbor

        boundaryPolygon = convexHull(listCoordinates)

        if not isPointInPolygon(p, boundaryPolygon):
            errmsg = "Point" + str(x) + "is not inside interpolated region!"
            print errmsg
            #raise Exception(errmsg)
            return "nan"

        distlist = [[abs(xel - p[0]) ** 2 + abs(yel - p[1]) ** 2, \
                    listIn.index([xel, yel, zel])] for xel, yel, zel in listIn]
        distlist.sort()     # list of distance and index number
        x1Index = distlist[0][1]
        x1, y1, z1 = listIn[x1Index]
        return z1

    evalPoint = np.array([ p ])
    lV = np.array(listValues)
    lC = np.array(listCoordinates)

    li=LinearNDInterpolator(lC, lV, 'nan')  #last argument is fill_value if outside boundaries!
    #licubic = CloughTocher2DInterpolator(lC, lV, 0. ) #  last one is fill_value, wich is returned is outside of boundaries

    result = li(evalPoint).item(0)

    #print result, licubic(evalPoint)
    return result




def Interpolate3D(listIn, x, err_out = "off"):
    '''linear interpolated function fint = f(x,y,z)


    listIn = [ [1, 1, 1, 2.123]  , [2, 1,1 ,-3.], [4,2, 3.,1.23] ,, [2.,2,2, `3.] ]

    and calculates linearly interpolated f at (x,y,z). Checks if interpolated
    point is inside of supplied list.
    '''

    listCoordinates = [[float(xi), float(yi),  float(zi)] for xi, yi, zi, fi in listIn]
    listIn = [[float(xi), float(yi), float(zi), float(fi)] for xi, yi, zi, fi in listIn]

    p = (float(x[0]), float(x[1]), float(x[2]))


    if _scipyloaded == False:

        for pis in listCoordinates:
            if p[0] == pis[0] and p[1] == pis[1] and p[2] == pis[2]:
                print "no interpolation needed"
                return listIn[listCoordinates.index(pis)][3]

        # use nearest neighbor instead

        distlist = [[abs(xel - p[0]) ** 2 + abs(yel - p[1]) ** 2 + abs(zel - p[2]) ** 2, \
                    listIn.index([xel, yel, zel, fel])] for xel, yel, zel, fel in listIn]

        distlist.sort()     # list of distance and index number
        x1Index = distlist[0][1]
        x1, y1, z1, f1 = listIn[x1Index]
        return f1


    # numpy is there

    evalPoint = np.array([ p ] )

    listValues = [float(fi) for xi, yi, zi, fi in listIn]

    lV = np.array(listValues)
    lC = np.array(listCoordinates)

    li=LinearNDInterpolator(lC, lV, 'nan')  #last argument is fill_value if outside boundaries!
    #licubic = CloughTocher2DInterpolator(lC, lV, 0. ) #  last one is fill_value!

    result = li(evalPoint).item(0)

    #print result, licubic(evalPoint)
    return result


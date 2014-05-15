#!/usr/bin/env python

__author__ = "Andreas Weiler <andreas.weiler@desy.de>"
__version__ = "0.1"


"""convexHull.py

Calculate the convex hull of a set of n 2D-points in O(n log n) time.
Taken from Berg et al., Computational Geometry, Springer-Verlag, 1997.

When run from the command line it generates a random set of points
inside a square of given length and finds the convex hull for those

Usage:

sam = [(121, 107),
 (14, 88),
 (123, 166),
 (57, 118),
 (85, 181),
 (82, 29),
 (143, 85),
 (98, 73),
 (86, 11),
 (110, 195),
 (47, 147),
 (121, 95),
 (2, 64),
 (142, 141),
 (4, 39),
 (22, 59),
 (160, 106),
 (158, 21),
 (47, 9)]

find convex Hull:

c = convexHull(sam)

see if point is in Poly

_isPointInPolygon( (x,y) , c)



"""


import sys, string, random


######################################################################
# Helpers
######################################################################

def _myDet(p, q, r):
    """Calc. determinant of a special matrix with three 2D points.

    The sign, "-" or "+", determines the side, right or left,
    respectivly, on which the point r lies, when measured against
    a directed vector from p to q.
    """

    # We use Sarrus' Rule to calculate the determinant.
    # (could also use the Numeric package...)
    sum1 = q[0]*r[1] + p[0]*q[1] + r[0]*p[1]
    sum2 = q[0]*p[1] + r[0]*q[1] + p[0]*r[1]

    return sum1 - sum2


def _isRightTurn((p, q, r)):
    "Do the vectors pq:qr form a right turn, or not?"

    assert p != q and q != r and p != r

    if _myDet(p, q, r) < 0:
	return 1
    else:
        return 0


######################################################################
# Public interface
######################################################################

def convexHull(P):
    "Calculate the convex hull of a set of points."

    # Get a local list copy of the points and sort them lexically.
    points = map(None, P)
    points.sort()

    # Build upper half of the hull.
    upper = [points[0], points[1]]
    for p in points[2:]:
	upper.append(p)
	while len(upper) > 2 and not _isRightTurn(upper[-3:]):
	    del upper[-2]

    # Build lower half of the hull.
    points.reverse()
    lower = [points[0], points[1]]
    for p in points[2:]:
	lower.append(p)
	while len(lower) > 2 and not _isRightTurn(lower[-3:]):
	    del lower[-2]

    # Remove duplicates.
    del lower[0]
    del lower[-1]

    # Concatenate both halfs and return.
    return tuple(upper + lower)



def isPointInPolygon(r, P):
    "Is point r inside a given polygon P?"

    # We assume the polygon is a list of points, listed clockwise!
    for i in xrange(len(P[:-1])):
        p, q = P[i], P[i+1]
        if not _isRightTurn((p, q, r)):
            return 0 # Out!

    return 1 # It's within!


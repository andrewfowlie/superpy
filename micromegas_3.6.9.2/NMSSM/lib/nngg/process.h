* process.h
* defines all process-dependent parameters
* this file is part of FormCalc
* last modified 15 Jan 03 th

* Definition of the external particles.
* The TYPEn may be one of SCALAR, FERMION, PHOTON, or VECTOR.
* (PHOTON is equivalent to VECTOR, except that longitudinal
* modes are not allowed)

* Note: The initial definitions for particles 2...5 are of course
* sample entries for demonstration purposes.

#define TYPE1 FERMION
#define MASS1 MNE1
#define CHARGE1 0

#define TYPE2 FERMION
#define MASS2 MNE1
#define CHARGE2 0

#define TYPE3 PHOTON
#define MASS3 0
#define CHARGE3 0

#define TYPE4 PHOTON
#define MASS4 0
#define CHARGE4 0


#define TYPE5 PHOTON
#define MASS5 0
#define CHARGE5 0

* The combinatorial factor for identical particles in the final state:
* 1/n! for n identical particles, 1 otherwise

#define IDENTICALFACTOR 5D-1

* Possibly a colour factor, e.g.
* - an additional averaging factor if any of the incoming particles
*   carry colour,
* - the overall colour factor resulting from the external particles
*   if that cannot computed by FormCalc (e.g. if the model has no
*   colour indices, as SM.mod).

#define COLOURFACTOR 1

* Whether to include soft-photon bremsstrahlung.
* ESOFTMAX is the maximum energy a soft photon may have and may be
* defined in terms of sqrtS, the CMS energy.

c#define BREMSSTRAHLUNG
#define ESOFTMAX .1D0*sqrtS

* Possibly some wave-function renormalization
* (e.g. if calculating in the background-field method)

c#define WF_RENORMALIZATION (nW*dWFW1 + nZ*dWFZ1)

#define NCOMP 2

* Include the kinematics-dependent part of the code

#include "2to2.F"


*     LanHEP output produced at Tue Apr 23 15:14:21 2013
*     Model named 'Z3 Inert Dublet'

      double precision Sqrt2, pi, degree, hbar_c2,bogus
      parameter (Sqrt2=1.41421356237309504880168872421D0)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      parameter (degree = pi/180D0)
      parameter (hbar_c2 = 3.8937966D8)
      parameter (bogus = -1D123)
      double complex cI
      parameter (cI = (0D0, 1D0))

      double precision Divergence
      common /renorm/ Divergence

      double precision EE, GG, SW, CW, MZ, MW, Q, Mb, V, wZ, wW
      double precision Mm, Mt, Mu, Md, Mc, Ms, Mtop, wtop, Mh
      double precision la3, la2, la4, laS, laS1, laS2, laS21, Mdm1
      double precision Mdm2, sinDm, muS, cosDm, muSq, mush, mu2q
      double precision MHC

      double precision AAABR(99)

      common /mdl_para/
     &    EE, GG, SW, CW, MZ, MW, Q, Mb, V, wZ, wW, Mm, Mt, Mu,
     &    Md, Mc, Ms, Mtop, wtop, Mh, la3, la2, la4, laS, laS1,
     &    laS2, laS21, Mdm1, Mdm2, sinDm, muS, cosDm, muSq, mush,
     &    mu2q, MHC, AAABR
      double precision slhaValFormat
      external slhaValFormat


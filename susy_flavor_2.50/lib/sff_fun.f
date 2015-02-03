c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: SFF_FUN.F
c     Released: 08: 6: 2013(J.R.)
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains standard functions for the                   c 
c     scalar(pseudoscalar)-fermion-fermion vertex formfactor          c
c                                                                     c
c                                                                     c
c                            V1______/_____ f_j                       c
c                            /|      \  q                             c
c                   k+p1+p2 / |                                       c 
c                         |/  |                                       c
c                         /~  |                                       c
c                     L2 /    |                                       c
c             S   /     /    /|\ k                                    c
c             ~~~~~~~~~~\V2   |                                       c
c             p2  \      \    | L1                                    c
c                     L3  \|  |                                       c 
c                         ~\  |                                       c
c                      k+p1 \ |                                       c
c                            \|______\______ f_i                      c
c                             V3     /  p1                            c
c                                                                     c
c     General form of the vertex                                      c
c                                                                     c


c     V = V_tree + i(V_L P_L + V_R P_R)                               c
c                                                                     c
c     Momentum arguments in formfactors: p1^2=fmi^2                   c
c                                        q^2=fmj^2                    c
c                                        p2^2=sm^2                    c
c                                                                     c
c     Other arguments:                                                c
c     dm1,dm2,dm3:     masses of particles circulating in loop        c 
c                      on lines L1,L2,L3 respectively                 c
c     vl_i,vr_i: LR couplings in the vertices (complex in general)    c
c     form:      complex output array containing formfactor values    c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine sff_svert(dm1,dm2,dm3,vl2,vr2,vl1,vr1,vl3,vr3,form)
c     Fermion-fermion-scalar loop
c     Vertices in loop:
c       V1: i (vl1 PL + vr1 PR)
c       V2: i (vl2 PL + vr2 PR)
c       V3: i (vl3 PL + vr3 PR)
      implicit double precision (a-h,o-z)
      double complex vl1,vr1,vl2,vr2,vl3,vr3
      double complex form(2)
      double complex c0,c11,c12,b0
      common/sff_args/sm,fmj,fmi
      p2 = fmi**2
      q2 = sm**2
      pq = (fmj**2 - sm**2 - fmi**2)/2
c     Left formfactor
      form(1) = (dm3*fmj*vl2*vl3*vr1 + dm1**2*vl1*vl3*vr2 + fmi*fmj*vl2
     $     *vr1*vr3 + dm2*vl1*vl2*(dm3*vl3 + fmi*vr3))*c0(p2,q2,pq,dm1
     $     ,dm3,dm2) + fmi*(fmi*vl1*vl3*vr2 + (dm2*vl1*vl2 + fmj*vl2
     $     *vr1+ dm3*vl1*vr2)*vr3)*c11(p2,q2,pq,dm1,dm3,dm2) + (fmj*vl3
     $     *(fmj*vl1 + dm2*vr1)*vr2 - fmi**2*vl1*vl3*vr2 - dm2*fmi*vl1
     $     *vl2*vr3 + dm3*(fmj*vl2*vl3*vr1 - fmi*vl1*vr2*vr3))*c12(p2,q2
     $     ,pq,dm1,dm3,dm2) - vr2*vl1*vl3*b0(q2,dm3,dm2)
c     Right formfactor
      form(2) = (fmi*vl3*(fmj*vl1 + dm2*vr1)*vr2 + (dm3*(fmj*vl1 + dm2
     $     *vr1)*vr2 + dm1**2*vl2*vr1)*vr3)*c0(p2,q2,pq,dm1,dm3,dm2) +
     $     fmi*(dm3*vl2*vl3*vr1 + fmj*vl1*vl3*vr2 + dm2*vl3*vr1*vr2 +
     $     fmi*vl2*vr1*vr3)*c11(p2,q2,pq,dm1,dm3,dm2) + ((fmj**2 - fmi
     $     **2)*vl2*vr1*vr3 + dm2*(fmj*vl1*vl2*vr3 - fmi*vl3*vr1*vr2) +
     $     dm3*(fmj*vl1*vr2*vr3 - fmi*vl2*vl3*vr1))*c12(p2,q2,pq,dm1
     $     ,dm3,dm2) - vl2*vr1*vr3*b0(q2,dm3,dm2) 
      return  
      end

      subroutine fss_svert(dm1,dm2,dm3,vl1,vr1,vl3,vr3,form)
c     Scalar-scalar-fermion loop
c     Vertices in loop:
c       V1: i (vl1 PL + vr1 PR)
c       V2: i a
c       V3: i (vl3 PL + vr3 PR)
c     a  - factorized constant
      implicit double precision (a-h,o-z)
      double complex vl1,vr1,vl3,vr3
      double complex form(2)
      double complex c0,c11,c12
      common/sff_args/sm,fmj,fmi
      p2 = fmi**2
      q2 = sm**2
      pq = (fmj**2 - fmi**2 - sm**2)/2
c     Left formfactor
      form(1) = dm1*vl1*vl3*c0(p2,q2,pq,dm1,dm3,dm2) - fmi*vl1*vr3
     $     *c11(p2,q2,pq,dm1,dm3,dm2) + (fmi*vl1*vr3 - fmj*vl3*vr1)
     $     *c12(p2,q2,pq,dm1,dm3,dm2)
c     Right formfactor
      form(2) = dm1*vr1*vr3*c0(p2,q2,pq,dm1,dm3,dm2) - fmi*vl3*vr1
     $     *c11(p2,q2,pq,dm1,dm3,dm2) + (fmi*vl3*vr1 - fmj*vl1*vr3)
     $     *c12(p2,q2,pq,dm1,dm3,dm2)
      return
      end

      subroutine vff_svert(dm1,dm2,dm3,vl2,vr2,form)
c     Fermion-fermion-vector boson loop
c     Vertices in loop:
c       V1: i A G(al) P_L  
c       V2: i (vl2 PL + vr2 PR)
c       V3: i B G(be) P_L
c     A, B - factorized constants
      implicit double precision (a-h,o-z)
      double complex vl2,vr2
      double complex form(2)
      double complex c0,c11,c12
      common/sff_args/sm,fmj,fmi
      p2 = fmi**2
      q2 = sm**2
      pq = (fmj**2 - sm**2 - fmi**2)/2
c     Left formfactor
      form(1) =  2*(dm2*fmi*vr2*c0(p2,q2,pq,dm1,dm3,dm2) + fmi*(dm3*vl2
     $     + dm2*vr2)*c11(p2,q2,pq,dm1,dm3,dm2) - fmi*(dm3*vl2 + dm2
     $     *vr2)*c12(p2,q2,pq,dm1,dm3,dm2))
c     Right formfactor
      form(2) = 2*(dm3*fmj*vl2*c0(p2,q2,pq,dm1,dm3,dm2) + fmj*(dm3*vl2 +
     $     dm2*vr2)*c12(p2,q2,pq,dm1,dm3,dm2))
      return
      end

      subroutine fvv_svert(dm1,dm2,dm3,form)
c     Vector boson-vector boson-fermion loop
c     Vertices in loop:
c       V2: i A G(al) P_L
c       V1: i S 
c       V3: i B G(be) P_L
c     A,B,S factorized
      implicit double precision (a-h,o-z)
      double complex form(2)
      double complex c11,c12
      common/sff_args/sm,fmj,fmi
      p2 = fmi**2
      q2 = sm**2
      pq = (fmj**2 - sm**2 - fmi**2)/2
c     Left formfactor
      form(1) = 2*fmi*(c11(p2,q2,pq,dm1,dm3,dm2) 
     $     - c12(p2,q2,pq,dm1,dm3,dm2))
c     Right formfactor
      form(2) = 2*fmj*c12(p2,q2,pq,dm1,dm3,dm2)
      return
      end

      subroutine fsv_svert(dm1,dm2,dm3,vl1,vr1,form)
c     Vertices in loop:
c       V1: i (vl1 PL + vr1 PR)
c       V2: i cv (p_S + P_V)
c       V3: i A G(al) P_L
c     A,cv - factorized constants
      implicit double precision (a-h,o-z)
      double complex vl1,vr1
      double complex form(2)
      double complex c0,c11,c12,b0
      common/sff_args/sm,fmj,fmi
      p2 = fmi**2
      q2 = sm**2
      pq = (fmj**2 - sm**2 - fmi**2)/2
c     Left formfactor
      form(1) = dm1*fmi*vl1*c0(p2,q2,pq,dm1,dm3,dm2) - fmi*(dm1*vl1 + 2
     $     *fmj*vr1)*c11(p2,q2,pq,dm1,dm3,dm2) + fmi*(dm1*vl1 + fmj*vr1)
     $     *c12(p2,q2,pq,dm1,dm3,dm2)
c     Right formfactor
      form(2) = dm1*(dm1*vr1 - 2*fmj*vl1)*c0(p2,q2,pq,dm1,dm3,dm2) +
     $     (fmi**2 + 2*fmj**2 - 2*sm**2)*vr1*c11(p2,q2,pq,dm1,dm3,dm2) -
     $     (dm1*fmj*vl1 + (fmi**2 - 2*sm**2)*vr1)*c12(p2,q2,pq,dm1,dm3
     $     ,dm2) - vr1*b0(sm**2,dm3,dm2)
      return
      end

      subroutine fvs_svert(dm1,dm2,dm3,vl3,vr3,form)
c     Vector boson-scalar-fermion loop
c     Vertices in loop:
c       V1: i A G(al) P_L
c       V2: i cv (p_S + P_V)
c       V3: i (vl3 PL + vr3 PR)
      implicit double precision (a-h,o-z)
      double complex vl3,vr3
      double complex form(2)
      double complex c0,c11,c12,b0
      common/sff_args/sm,fmj,fmi
      p2 = fmi**2
      q2 = sm**2
      pq = (fmj**2 - sm**2 - fmi**2)/2
c     Left formfactor
      form(1) = dm1*(2*fmi*vr3 - dm1*vl3)*c0(p2,q2,pq,dm1,dm3,dm2) + fmi
     $     *(dm1*vr3 - 2*fmi*vl3)*c11(p2,q2,pq,dm1,dm3,dm2) - (fmj**2
     $     *vl3 - 2*sm**2*vl3 + dm1*fmi*vr3)*c12(p2,q2,pq,dm1,dm3,dm2)
     $     + vl3*b0(sm**2,dm3,dm2)
c     Right formfactor
      form(2) =  - dm1*fmj*vr3*c0(p2,q2,pq,dm1,dm3,dm2) + fmi*fmj*vl3
     $     *c11(p2,q2,pq,dm1,dm3,dm2) + fmj*(fmi*vl3 + dm1*vr3)*c12(p2
     $     ,q2,pq,dm1,dm3,dm2)
      return
      end


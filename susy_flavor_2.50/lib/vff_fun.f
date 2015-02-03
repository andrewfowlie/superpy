c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: VFF_FUN.F
c     Released: 30: 4: 2014(J.R.)
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains standard functions for the                   c 
c     gluon(photon)-fermion-fermion vertex formfactor                 c
c     on-shell version only (dipole moments)                          c
c                                                                     c
c                            V1______/_____ f_i                       c
c                            /|      \ pi                             c
c                        k  / |                                       c 
c                         |/  |                                       c
c                         /~  |                                       c
c                     L2 /    |                                       c
c            V^mu /     /    /|\ k-pi                                 c
c             ~~~~~~~~~~\V3   |                                       c
c             q   \      \    | L1                                    c
c                     L3  \|  |                                       c 
c                         ~\  |                                       c
c                      k-q  \ |                                       c
c                            \|______\______ f_j                      c
c                             V2     /  pj                            c
c                                                                     c
c     General form of the vertex                                      c
c                                                                     c
c     V = V_tree - sigma^munu (V_L P_L + V_R P_R) q_nu                c
c                                                                     c
c     Momentum arguments in formfactors: pi^2=fmi^2                   c
c                                        pj^2=fmj^2                   c
c                                        q^2=0                        c
c                                                                     c
c     Other arguments:                                                c
c     dm1,dm2:   masses of particles circulating in loop on lines     c
c     L1,L2,L3 respectively                                           c
c     vl_i,vr_i: LR couplings in the vertices (complex in general)    c
c     form:      complex output array containing formfactor values    c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine vsf_vert(fm,sm,vf,vs,vl1,vr1,vl2,vr2,form)
c     Fermion-fermion-scalar + scalar-scalar-fermion + self-energy loops
c     Vertices in loop:
c       V1: i (vl1 PL + vr1 PR)
c       V2: i (vl2 PL + vr2 PR)
c       vf,vs: effective scalar(fermion)-photon(gluon) coupling
c       in case of gluon contain resummed color factors
      implicit double precision (a-h,o-z)
      double complex vl1,vr1,vl2,vr2
      double complex form(2)
      double complex c0,f
      common/vff_args/fmi
      pq = fmi**2/2
c     Left tensor formfactor
      form(1) = (vf - vs)*vl2*vr1/2 + (vf + vs)*vl2*(fm*fmi*vl1 + (fm**2
     $     - sm**2)*vr1)*f(fmi**2,sm,fm)/fmi**2 + vs*vl2*vr1*sm**2
     $     *c0(0.d0,0.d0,pq,sm,sm,fm) - vf*vl2*(fmi*vl1 + fm*vr1)*fm
     $     *c0(0.d0,0.d0,pq,fm,fm,sm) 
c     Right tensor formfactor
      form(2) = (vf - vs)*vl1*vr2/2 + (vf + vs)*vr2*(fm*fmi*vr1 + (fm**2
     $     - sm**2)*vl1)*f(fmi**2,sm,fm)/fmi**2 + vs*vl1*vr2*sm**2
     $     *c0(0.d0,0.d0,pq,sm,sm,fm) - vf*vr2*(fm*vl1 + fmi*vr1)*fm
     $     *c0(0.d0,0.d0,pq,fm,fm,sm) 
      do kk=1,2
         form(kk) = form(kk)/fmi
      end do
      return  
      end

      subroutine vvf_vert(fm,vm,vf,vl1,vl2,form)
c     Fermion-fermion-W + W-fermion self-energy loops
c     Vertices in loop:
c       V1: i vl1 Gamma^al PL
c       V2: i vl2 Gamma^be PL
c       vf: fermion-photon(gluon) coupling
      implicit double precision (a-h,o-z)
      double complex vl1,vl2
      double complex form(2)
      double complex c0,f
      common/vff_args/fmi
      pq = fmi**2/2
c     Left tensor formfactor
      form(1) = 0.d0
c     Right tensor formfactor
      form(2) = vf*vl1*vl2*(1 - 2*fm**2*c0(0.d0,0.d0,pq,fm,fm,vm) 
     $     + 2*(fm**2 + fmi**2 - vm**2)*f(fmi**2,vm,fm)/fmi**2)/fmi
      return  
      end

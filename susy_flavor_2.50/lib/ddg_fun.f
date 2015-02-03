c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: DDG_FUN.F
c     Released: 02:07:1996 (J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for the coefficients          c
c     of the expansion of the generic diagrams with scalar or      c
c     vector exchange into the gauge-invariant operators           c
c     H^1_L,R - H^5_L,R for the d^J -> d^I + photon/gluon decay    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dd_ffs(ai,aj,bi,bj,qf,fm,sm,cfl,cfr)
c     Scalar-fermion-fermion loop (i/(4pi)^2 factorized)
c     fbar-d^J-S   vertex: i(aj P_L + bj P_R)
c     dbar^I-f-S^* vertex: i(ai P_R + bi P_L)
c     fbar-f-gamma vertex: ie qf gamma(mu)
c     fm,sm: fermion and scalar mass in loop
c     cfl(1..5),cfr(1..5): coefficients of H^i_L and H^i_R
      implicit double precision (a-h,o-z)
      double complex ai,aj,bi,bj
      double complex cfl(5),cfr(5)
      cfl(1) = cfl(1) + qf*aj*bi*fm*cp11(fm,sm)
      cfl(2) = cfl(2) + qf*aj*bi*fm*cp12(fm,sm)
      cfl(3) = cfl(3) + qf*ai*aj*(2*cp11(fm,sm) - 2*cp21(fm,sm)
     $                - cp23(fm,sm))
      cfl(4) = cfl(4) - qf*ai*aj*cp23(fm,sm)/2
      cfl(5) = cfl(5) + qf*ai*aj*cp22(fm,sm)
      cfr(1) = cfr(1) + qf*bj*ai*fm*cp11(fm,sm)
      cfr(2) = cfr(2) + qf*bj*ai*fm*cp12(fm,sm)
      cfr(3) = cfr(3) + qf*bi*bj*(2*cp11(fm,sm) - 2*cp21(fm,sm)
     $                - cp23(fm,sm))
      cfr(4) = cfr(4) - qf*bi*bj*cp23(fm,sm)/2
      cfr(5) = cfr(5) + qf*bi*bj*cp22(fm,sm)
      return
      end

      subroutine dd_ssf(ai,aj,bi,bj,qs,fm,sm,cfl,cfr)
c     Fermion-scalar-scalar in loop (i/(4pi)^2 factorized)
c     fbar-d^J-S   vertex: i(aj P_L + bj P_R)
c     dbar^I-f-S^* vertex: i(ai P_R + bi P_L)
c     S*-S-gamma   vertex: ie qs (p_S + p_S*)
c     fm,sm: fermion and scalar mass in loop
c     cfl(1..5),cfr(1..5): coefficients of H^i_L and H^i_R
      implicit double precision (a-h,o-z)
      double complex ai,aj,bi,bj
      double complex cfl(5),cfr(5)
      cfl(1) = cfl(1) + qs*aj*bi*fm*cp12(sm,fm)/2
      cfl(2) = cfl(2) - qs*aj*bi*fm*cp12(sm,fm)
      cfl(3) = cfl(3) + qs*ai*aj*(cp11(sm,fm) - 2*cp21(sm,fm)
     $                - cp23(sm,fm))
      cfl(4) = cfl(4) - qs*ai*aj*cp23(sm,fm)/2
      cfl(5) = cfl(5) - 2*qs*ai*aj*cp23(sm,fm)
      cfr(1) = cfr(1) + qs*bj*ai*fm*cp12(sm,fm)/2
      cfr(2) = cfr(2) - qs*bj*ai*fm*cp12(sm,fm)
      cfr(3) = cfr(3) + qs*bi*bj*(cp11(sm,fm) - 2*cp21(sm,fm)
     $                - cp23(sm,fm))
      cfr(4) = cfr(4) - qs*bi*bj*cp23(sm,fm)/2
      cfr(5) = cfr(5) + 2*qs*bi*bj*cp23(sm,fm)
      return
      end

      subroutine dd_ffv(ai,aj,qf,fm,vm,cfl,cfr)
c     Vector-fermion-fermion in loop (i/(4pi)^2 factorized)
c     fbar-d^J-V   vertex: -i aj gamma(al) P_L
c     dbar^I-f-V^* vertex: -i ai gamma(be) P_L
c     fbar-f-gamma vertex:  i e qf gamma(mu)
c     fm,vm: fermion and vector mass in loop
c     cfl(1..5),cfr(1..5): coefficients of H^i_L and H^i_R
      implicit double precision (a-h,o-z)
      double complex ai,aj
      double complex cfl(5),cfr(5)
      cfl(3) = cfl(3) + 2*qf*ai*aj*(cp0_1(fm,vm) - 2*cp21(fm,vm)
     $                - cp23(fm,vm))
      cfl(4) = cfl(4) + qf*ai*aj*(cp12(fm,vm) - cp23(fm,vm))
      cfl(5) = cfl(5) + 2*qf*ai*aj*cp22(fm,vm)
      return
      end

      subroutine dd_vvf(ai,aj,qv,fm,vm,cfl,cfr)
c     Fermion-vector-vector in loop (i/(4pi)^2 factorized)
c     fbar-d^J-V   vertex: -i aj gamma(al) P_L
c     dbar^I-f-V^* vertex: -i ai gamma(be) P_L
c     V*-V-gamma   vertex:  background field gauge
c       ie qv (2p(be)g(mu,al) - 2p(al)g(mu,be) + (p(mu)-2k(mu))g(al,be))
c     fm,vm: fermion and vector mass in loop
c     cfl(1..5),cfr(1..5): coefficients of H^i_L and H^i_R
      implicit double precision (a-h,o-z)
      double complex ai,aj
      double complex cfl(5),cfr(5)
      cfl(3) = cfl(3) + 2*qv*ai*aj*(cp23(vm,fm) + 2*cp21(vm,fm)
     $                - 3*cp11(vm,fm))
      cfl(4) = cfl(4) + qv*ai*aj*(cp23(vm,fm) - 2*cp11(vm,fm))
      cfl(5) = cfl(5) + 4*qv*ai*aj*cp23(vm,fm)
      return
      end



